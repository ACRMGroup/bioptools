/*************************************************************************

   Program:    distancematric
   File:       distancematrix.c
   
   Version:    V1.0
   Date:       01.12.16
   Function:   Calculate inter-CA distances on a set of common-labelled
               PDB files
   
   Copyright:  (c) UCL, Dr. Andrew C. R. Martin 2016
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   This is a very crude and simple program for calculating means and 
   standard deviations for inter residue distances from a set of PDB files
   having common numbering (e.g. antibodies). The program is *slow* as it
   doesn't do anything clever with the list of residue pair labels - it
   simply stores them in an array and searches that array every time to
   find a label. This really should be replaced by a hash function!

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   01.12.16  Original - Based loosely on distmat.c

*************************************************************************/
/* #define DEBUG 1 */

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "bioplib/MathUtil.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF          160
#define DEF_MAXRES       300     /* Max residues pairs to handle        */
#define MAXLABEL         16      /* ResidueLabel_ResidueLabel size      */
#define DEF_NCHAINS      2       /* Number of chains to keep            */

/************************************************************************/
/* Globals
*/
REAL *gSx       = NULL, 
     *gSxSq     = NULL;
int  *gNVal     = NULL,
     gNStrings  = 0;
char **gStrings = NULL;

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *nchains, int *maxrespairs, BOOL *singleFile);
BOOL HandleInput(FILE *in, FILE *out, int nrespairs, int maxchain,
                 BOOL singleFile);
PDB *GetPDB(char *filename, int maxchain);
PDB *GetPDBFp(FILE *fp, int maxchain);
BOOL AllocateMemory(int nrespairs);
void FreeMemory(int nrespairs);
BOOL ProcessPDB(PDB *pdb, int nrespairs);
void DisplayResults(FILE *out);
void Usage(void);
BOOL StoreData(char chain1, int res1, char ins1,
               char chain2, int res2, char ins2,
               REAL d, int nrespairs);
int  LookUpString(char *string, int tablesize);

#define ATOMS_CA  0
#define ATOMS_ALL 1
#define ATOMS_SC  2

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for distance analysis

   01.04.09 Original   By: ACRM
   06.04.09 Added -n and -m parameters
*/
int main(int argc, char **argv)
{
   char  infile[MAXBUFF],
         outfile[MAXBUFF];
   FILE  *in  = stdin,
         *out = stdout;
   BOOL  singleFile  = FALSE;
   ULONG hashSize    = DEF_MAXRES * DEF_MAXRES;
   int   atomTypes   = ATOMS_CA;
   

   if(ParseCmdLine(argc, argv, infile, outfile, &singleFile, &atomTypes))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if(hashTable = blInitializeHash(hashSize))
         {
            if(HandleInput(in, out, singleFile, hashTable, atomTypes))
            {
               DisplayResults(out, hashTable);
            }
         }
         else
         {
            fprintf(stderr,"ERROR: Unable to initialize hash table.\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"ERROR: Unable to open files.\n");
         return(1);
      }
   }
   else
   {
      Usage();
   }

   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *singleFile, int *atomTypes)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            BOOL   *singleFile  Input is a single PDB file instead of
                                a list
   Returns: BOOL                Success?

   Parse the command line

   01.04.09 Original    By: ACRM
   06.04.09 Added -n and -m and their parameters
   30.11.16 Added -p
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *singleFile, int *atomTypes)
{
   argc--;
   argv++;

   infile[0]    = outfile[0] = '\0';
   *singleFile  = FALSE;
   
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'p':
            *singleFile = TRUE;
            break;
         case 'c':
            *atomTypes = ATOMTYPES_CA;
            break;
         case 'a':
            *atomTypes = ATOMTYPES_ALL;
            break;
         case 's':
            *atomTypes = ATOMTYPES_SC;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are only 1 or 2 arguments left             */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         /* If there's another, copy it to outfile                      */
         argc--;
         argv++;
         if(argc)
            strcpy(outfile, argv[0]);
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL HandleInput(FILE *in, FILE *out, BOOL singleFile, HASHTABLE *hashTable)
   ------------------------------------------------------------------
   Handle the input file - extract the PDB filenames and process each 
   in turn

   01.04.09 Original   By: ACRM
   06.04.09 Handles maxchain
   30.11.16 Added singleFile
*/
BOOL HandleInput(FILE *in, FILE *out, BOOL singleFile, HASHTABLE *hashTable, int atomTypes)
{
   char filename[MAXBUFF];
   
   PDB  *pdb;
   
   if(singleFile)
   {
      ProcessFile(in, hashTable, atomTypes);
   }
   else
   {
      while(fgets(filename,MAXBUFF,in))
      {
         FILE *fp;
         
         TERMINATE(filename);
         
         if(fp = fopen(filename,"r"))
         {
            fprintf(stderr,"INFO: Processing file: %s\n",filename);
            ProcessFile(fp, hashTable);
            fclose(fp);
         }
         else
         {
            fprintf(stderr,"WARNING: Unable to read file: %s\n", filename);
         }
      }
   }
   
   return(TRUE);
}

void ProcessFile(FILE *fp, HASHTABLE *hashTable, int atomTypes)
{
   PDB *pdb;
   int natoms;
   
   if((pdb = blReadPDBAtoms(fp, &natoms))!=NULL)
   {
      pdb = ReduceAtomList(pdb, atomTypes);
      ProcessPDB(pdb, hashTable);
      FREELIST(pdb, PDB);
   }
   else
   {
      fprintf(stderr,"Error: No memory to read PDB file\n");
      exit(1);
   }
}

void ProcessPDB(PDB *pdb, HASHTABLE hashTable)
{
   PDB *p, *q;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(q=pdb; q!=NULL; NEXT(q))
      {
         
      }
   }
   
}



PDB *ReduceAtomList(PDB *pdb, int atomTypes)
{
   PDB  *reduced;
   char *sel[40];
   
   switch(atomTypes)
   {
   case ATOMTYPES_ALL:
      return(pdb);
      break;
   case ATOMTYPES_CA:
      SELECT(sel[0],"CA  ");
      if(sel[0] == NULL) return(NULL);
      reduced = SelectAtomsPDB(pdb, 1, sel, &natoms);
      FREELIST(pdb, PDB);
      return(reduced);
      break;
   case ATOMTYPES_SC:
      SELECT(sel[0], "CB  ");
      SELECT(sel[1], "CD  ");
      SELECT(sel[2], "CD1 ");
      SELECT(sel[3], "CD2 ");
      SELECT(sel[4], "CE  ");
      SELECT(sel[5], "CE1 ");
      SELECT(sel[6], "CE2 ");
      SELECT(sel[7], "CE3 ");
      SELECT(sel[8], "CG  ");
      SELECT(sel[9], "CG1 ");
      SELECT(sel[10],"CG2 ");
      SELECT(sel[11],"CH2 ");
      SELECT(sel[12],"CZ  ");
      SELECT(sel[13],"CZ2 ");
      SELECT(sel[14],"CZ3 ");
      SELECT(sel[15],"ND1 ");
      SELECT(sel[16],"ND2 ");
      SELECT(sel[17],"NE  ");
      SELECT(sel[18],"NE1 ");
      SELECT(sel[19],"NE2 ");
      SELECT(sel[20],"NH1 ");
      SELECT(sel[21],"NH2 ");
      SELECT(sel[22],"NZ  ");
      SELECT(sel[23],"OD1 ");
      SELECT(sel[24],"OD2 ");
      SELECT(sel[25],"OE1 ");
      SELECT(sel[26],"OE2 ");
      SELECT(sel[27],"OG  ");
      SELECT(sel[28],"OG1 ");
      SELECT(sel[29],"OH  ");
      SELECT(sel[30],"SD  ");
      SELECT(sel[31],"SG  ");
      if(sel[31] == NULL) return(NULL);
      reduced = SelectAtomsPDB(pdb, 32, sel, &natoms);
      FREELIST(pdb, PDB);
      return(reduced);
      break;
   }
   return(pdb);
}







/************************************************************************/
/*>BOOL ProcessPDB(PDB *pdb, int nrespairs)
   ----------------------------------------
   Does the main work of handling a PDB linked list and getting the
   distances then calling the code ro store them.

   01.04.09 Original   By: ACRM
*/
BOOL ProcessPDB(PDB *pdb, int nrespairs)
{
   REAL d;
   PDB *p, *q;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(q=pdb; q!=NULL; NEXT(q))
      {
         
         if(p==q)
         {
            d = (REAL)0.0;
         }
         else
         {
            d = DIST(p, q);
         }
         if(!StoreData(p->chain[0], p->resnum, p->insert[0],
                       q->chain[0], q->resnum, q->insert[0],
                       d, nrespairs))
         {
            return(FALSE);
         }
      }
   }
   return(TRUE);
}


/************************************************************************/
/*>BOOL StoreData(char chain1, int res1, char ins1,
                  char chain2, int res2, char ins2,
                  REAL d, int nrespairs)
   ------------------------------------------------
   Takes two residue specifications. Calls the routine to find an index
   by which to store the data for mean and SD calculation and then stores 
   it.

   01.04.09 Original   By: ACRM
*/
BOOL StoreData(char chain1, int res1, char ins1,
               char chain2, int res2, char ins2,
               REAL d, int nrespairs)
{
   char resids[32];
   int  i;
   
   sprintf(resids, "%c%d%c %c%d%c", 
           chain1, res1, ins1,
           chain2, res2, ins2);
   i = LookUpString(resids, nrespairs);
   if(i==(-1))
   {
      return(FALSE);
   }
   
   CalcExtSD(d, 0, &gSx[i], &gSxSq[i], &gNVal[i], NULL, NULL);
   return(TRUE);
}


/************************************************************************/
/*>int LookUpString(char *string, int tablesize)
   ---------------------------------------------
   A very-crude string indexing function. Simply looks in the table of
   strings to see if it is there already - if so returns the index for
   the string. If not, then it adds it onto the end of the table and
   returns its index.

   Really this should be replaced by some sensible and clever hashing
   function!

   01.04.09 Original   By: ACRM
*/
int LookUpString(char *string, int tablesize)
{
   int i;
   int l;
   l = strlen(string);
   
   /* See if string is already in the table                             */
   for(i=0; i<gNStrings; i++)
   {
      if(!strncmp(gStrings[i], string, l))
      {
#ifdef DEBUG
         fprintf(stderr,"DEBUG: Found in table: '%s' (%d)\n", string, i);
#endif
         return(i);
      }
   }
   /* If not, then add it to the table                                  */
   if(gNStrings < tablesize)
   {
      i=gNStrings++;
#ifdef DEBUG
      fprintf(stderr,"DEBUG: Added to table: '%s' (%d)\n", string, i);
#endif
      strncpy(gStrings[i], string, l+1);
      return(i);
   }
   
   return(-1);
}


/************************************************************************/
/*>void DisplayResults(FILE *out)
   ------------------------------
   Display the results. Run through each indexed location and calculate 
   the mean and standard deviation then print the residue IDs with these
   values.

   01.04.09 Original   By: ACRM
*/
void DisplayResults(FILE *out)
{
   int  i;
   REAL mean, 
        sd;

   for(i=0; i<gNStrings; i++)
   {
      CalcExtSD((REAL)0.0, 1, 
                &gSx[i], &gSxSq[i], &gNVal[i], &mean, &sd);
      fprintf(out,"%s %6.3f %6.3f\n", 
              gStrings[i], mean, sd);
   }            
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   01.04.09 Original   By: ACRM
   06.04.09 V1.1
   30.11.16 V1.2
*/
void Usage(void)
{
   fprintf(stderr,"\nDistMat V1.2 (c) 2009, Dr. Andrew C.R. Martin, \
UCL\n");

   fprintf(stderr,"\nUsage: distmat [-n numchain] [-m numres] \
[input [output]]\n");
   fprintf(stderr,"       -n Specify max number of chains to keep \
(default: %d)\n", DEF_NCHAINS);
   fprintf(stderr,"       -m Specify max number of residues in PDB file \
(default: %d)\n", DEF_MAXRES);
   fprintf(stderr,"I/O Through stdin/stdout if not specified\n");

   fprintf(stderr,"\nDistMat analyses inter-CA distances in the first \
'numchain' chains\n");
   fprintf(stderr,"(default %d) of a PDB file containing at most \
'numres' residues\n", DEF_MAXRES);
   fprintf(stderr,"(default %d) in total across the chains. Defaults \
are designed for\n", DEF_NCHAINS);
   fprintf(stderr,"analyzing antibodies.\n");

   fprintf(stderr,"\nThe input file simply contains a list of the PDB \
files to be processed.\n\n");
}
