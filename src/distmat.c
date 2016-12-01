/*************************************************************************

   Program:    distmat
   File:       distmat.c
   
   Version:    V2.0
   Date:       01.12.16
   Function:   Calculate inter-CA distances on a set of common-labelled
               PDB files
   
   Copyright:  (c) UCL, Dr. Andrew C. R. Martin 2009-2016
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

   This is a simple program for calculating means and standard
   deviations for inter residue distances from a set of PDB files
   having common numbering (e.g. antibodies). Input is either a single
   PDB file or a 'file of files' - i.e. a file containing a list of
   PDB files to be processed.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   01.04.09  Original
   V1.1   06.04.09  Added -n and -m options
   V1.2   30.11.16  Minor cleanup for bioptools - Added option to take
                    a single PDB file as input rather than a set (-p)
   V2.0   01.12.16  Rewrite to use hashes. Removed -n and -m since these
                    are't needed any more - everything is dynamically
                    allocated.

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
#include "bioplib/hash.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF     160
#define DEF_MAXRES  300     /* Approximate number of residues in a file.
                               This is simply used to initialize the hash
                               table. Increasing it will increase 
                               efficiency for large structures, but waste
                               memory for small ones
                            */
#define MAXLABEL    16      /* ResidueLabel size                        */
#define ATOMS_CA    0       /* Selection types                          */
#define ATOMS_ALL   1
#define ATOMS_SC    2

typedef struct respair
{
   REAL sx;
   REAL sxsq;
   int  nval;
}  RESPAIR;
   

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *singleFile, int *atomTypes);
BOOL HandleInput(FILE *in, FILE *out, BOOL singleFile, 
                 HASHTABLE *hashTable, int atomTypes);
void ProcessFile(FILE *fp, HASHTABLE *hashTable, int atomTypes);
void ProcessPDB(PDB *pdb, HASHTABLE *hashTable);
void StoreData(HASHTABLE *hashTable, PDB *res1, PDB *res2, REAL dist);
PDB *ReduceAtomList(PDB *pdb, int atomTypes);
void DisplayResults(FILE *out, HASHTABLE *hashTable);
void Usage(void);

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
   HASHTABLE *hashTable = NULL;
   

   if(ParseCmdLine(argc, argv, infile, outfile, &singleFile, &atomTypes))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((hashTable = blInitializeHash(hashSize))!=NULL)
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
            *atomTypes = ATOMS_CA;
            break;
         case 'a':
            *atomTypes = ATOMS_ALL;
            break;
         case 's':
            *atomTypes = ATOMS_SC;
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
/*>BOOL HandleInput(FILE *in, FILE *out, BOOL singleFile, 
                    HASHTABLE *hashTable, int atomTypes)
   ------------------------------------------------------
   Handle the input file - extract the PDB filenames and process each 
   in turn

   01.04.09 Original   By: ACRM
   06.04.09 Handles maxchain
   30.11.16 Added singleFile
*/
BOOL HandleInput(FILE *in, FILE *out, BOOL singleFile, 
                 HASHTABLE *hashTable, int atomTypes)
{
   char filename[MAXBUFF];
   
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
         
         if((fp = fopen(filename,"r"))!=NULL)
         {
            fprintf(stderr,"INFO: Processing file: %s\n",filename);
            ProcessFile(fp, hashTable, atomTypes);
            fclose(fp);
         }
         else
         {
            fprintf(stderr,"WARNING: Unable to read file: %s\n", 
                    filename);
         }
      }
   }
   
   return(TRUE);
}


/************************************************************************/
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


/************************************************************************/
void ProcessPDB(PDB *pdb, HASHTABLE *hashTable)
{
   PDB *atom1, *atom2, 
       *res1,  *res1Next, 
       *res2,  *res2Next;
   REAL minDistSq = (REAL)0.0;

   /* Step through a residue at a time                                  */
   for(res1=pdb; res1!=NULL; res1=res1Next)
   {
      /* Find the start of the next residue                             */
      res1Next = blFindNextResidue(res1);
      
      /* Step through again a residue at a time                         */
      for(res2=pdb; res2!=NULL; res2=res2Next)
      {
         /* Find the start of the next residue                          */
         res2Next = blFindNextResidue(res2);

         /* Initialize minimum distance between the residues            */
         minDistSq = DISTSQ(res1, res2);

         /* Step through atoms in first residue                         */
         for(atom1=res1; atom1!=res1Next; NEXT(atom1))
         {
            /* Step through atoms in second residue to find 
               the minimum distance between the two residues
            */
            for(atom2=res2; atom2!=res2Next; NEXT(atom2))
            {
               REAL dSq = DISTSQ(atom1, atom2);
               if(dSq < minDistSq)
               {
                  minDistSq = dSq;
               }
            }
         }

         StoreData(hashTable, res1, res2, sqrt(minDistSq));
      }
   }
}


/************************************************************************/
void StoreData(HASHTABLE *hashTable, PDB *res1, PDB *res2, REAL dist)
{
   RESPAIR *rp = NULL;
   char resID1[MAXLABEL], resID2[MAXLABEL], resPair[MAXLABEL*2];
   
   MAKERESID(resID1, res1);
   MAKERESID(resID2, res2);
   sprintf(resPair, "%s-%s", resID1, resID2);
   
   /* This key not in the hash so create it */
   if(!blHashKeyDefined(hashTable, resPair))
   {
      if((rp = (RESPAIR *)malloc(sizeof(RESPAIR))) == NULL)
      {
         fprintf(stderr, "Error: no memory for RESPAIR data\n");
         exit(1);
      }

      rp->sx   = 0.0;
      rp->sxsq = 0.0;
      rp->nval = 0;

      blSetHashValuePointer(hashTable, resPair, (BPTR)rp);
   }
   else
   {
      rp = (RESPAIR *)blGetHashValuePointer(hashTable, resPair);
      if(rp == NULL)
      {
         fprintf(stderr, "Error: internal Hash confused!\n");
         exit(1);
      }
   }
   
   blCalcExtSD(dist, 0, &(rp->sx), &(rp->sxsq), &(rp->nval), NULL, NULL);
}


/************************************************************************/
PDB *ReduceAtomList(PDB *pdb, int atomTypes)
{
   PDB  *reduced;
   char *sel[40];
   int  natoms;
   
   switch(atomTypes)
   {
   case ATOMS_ALL:
      return(pdb);
      break;
   case ATOMS_CA:
      SELECT(sel[0],"CA  ");
      if(sel[0] == NULL) return(NULL);
      reduced = blSelectAtomsPDBAsCopy(pdb, 1, sel, &natoms);
      FREELIST(pdb, PDB);
      return(reduced);
      break;
   case ATOMS_SC:
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
      reduced = blSelectAtomsPDBAsCopy(pdb, 32, sel, &natoms);
      FREELIST(pdb, PDB);
      return(reduced);
      break;
   }
   return(pdb);
}






/************************************************************************/
/*>void DisplayResults(FILE *out, HASHTABLE *hashTable)
   -----------------------------------------------------
   Display the results. Run through each indexed location and calculate 
   the mean and standard deviation then print the residue IDs with these
   values.

   01.04.09 Original   By: ACRM
*/
void DisplayResults(FILE *out, HASHTABLE *hashTable)
{
   char **keys = NULL;
   int  i;
   REAL mean, 
        sd;


   keys = blGetHashKeyList(hashTable);
   
   for(i=0; keys[i] != NULL; i++)
   {
      RESPAIR *rp;
      char    res1[MAXLABEL],
              res2[MAXLABEL],
              *chp;
      
      rp = (RESPAIR *)blGetHashValuePointer(hashTable, keys[i]);

      strncpy(res1, keys[i], MAXLABEL);
      TERMAT(res1, '-');
      chp = strchr(keys[i], '-');
      strncpy(res2, chp+1, MAXLABEL);
      
      blCalcExtSD((REAL)0.0, 1, 
                  &(rp->sx), &(rp->sxsq), &(rp->nval), &mean, &sd);
      fprintf(out,"%s %s %6.3f %6.3f\n", 
              res1, res2, mean, sd);
   }            
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   01.04.09 Original   By: ACRM
   06.04.09 V1.1
   30.11.16 V1.2
   01.12.16 V2.0
*/
void Usage(void)
{
   fprintf(stderr,"\nDistMat V2.0 (c) 2009-2016, Dr. Andrew C.R. Martin, \
UCL\n");

   fprintf(stderr,"\nUsage: distmat [-p][-c][-a][-s] [input [output]]\n");
   fprintf(stderr,"       -p Input is a single PDB file instead \
of a file of files\n");
   fprintf(stderr,"       -c Only look at CAs (default)\n");
   fprintf(stderr,"       -a Look at all atoms\n");
   fprintf(stderr,"       -s Look at sidechain atoms\n");
   fprintf(stderr,"\nI/O Through stdin/stdout if not specified\n");

   fprintf(stderr,"\nDistMat analyses inter-CA distances in one or \
more PDB files\n");

   fprintf(stderr,"\nThe default input file simply contains a list of \
the PDB files to be processed.\n");
   fprintf(stderr,"If -p is specified the input is a single PDB \
file.\n\n");
}
