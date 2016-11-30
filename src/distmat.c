/*************************************************************************

   Program:    DistMat
   File:       distmat.c
   
   Version:    V1.1
   Date:       06.04.09
   Function:   Calculate inter-CA distances on a set of common-labelled
               antibody PDB files
   
   Copyright:  (c) Dr. Andrew C. R. Martin 2009, UCL
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
   V1.0   01.04.09  Original
   V1.1   06.04.09  Added -n and -m options

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
                  int *nchains, int *maxrespairs);
BOOL HandleInput(FILE *in, FILE *out, int nrespairs, int maxchain);
PDB *GetPDB(char *filename, int maxchain);
BOOL AllocateMemory(int nrespairs);
void FreeMemory(int nrespairs);
BOOL ProcessPDB(PDB *pdb, int nrespairs);
void DisplayResults(FILE *out);
void Usage(void);
BOOL StoreData(char chain1, int res1, char ins1,
               char chain2, int res2, char ins2,
               REAL d, int nrespairs);
int  LookUpString(char *string, int tablesize);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for distance analysis

   01.04.09 Original   By: ACRM
   06.04.09 Added -n and -m parameters
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   FILE *in  = stdin,
        *out = stdout;
   int  nchains = DEF_NCHAINS;
   int  maxrespairs = DEF_MAXRES * DEF_MAXRES;


   if(ParseCmdLine(argc, argv, infile, outfile, &nchains, &maxrespairs))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if(AllocateMemory(maxrespairs))
         {
            if(HandleInput(in, out, maxrespairs, nchains))
            {
               DisplayResults(out);
            }
            FreeMemory(maxrespairs);
         }
         else
         {
            /* Since Allocate may have partially succeed                */
            FreeMemory(maxrespairs);
            fprintf(stderr,"ERROR: Unable to allocate memory.\n");
            return(1);
         }
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
                     int *nchains, int *maxrespairs)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            int    *nchains     Number of chains to keep
            int    *maxrespairs Max number of res pairs to handle
   Returns: BOOL                Success?

   Parse the command line

   01.04.09 Original    By: ACRM
   06.04.09 Added -n and -m and their parameters
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *nchains, int *maxrespairs)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *nchains = DEF_NCHAINS;
   *maxrespairs = DEF_MAXRES * DEF_MAXRES;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'n':
            argv++;
            argc--;
            if((argc < 0) || (!sscanf(argv[0], "%d", nchains)))
            {
               return(FALSE);
            }
            break;
         case 'm':
            argv++;
            argc--;
            if((argc < 0) || (!sscanf(argv[0], "%d", maxrespairs)))
            {
               return(FALSE);
            }
            *maxrespairs = (*maxrespairs)*(*maxrespairs);
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
/*>BOOL HandleInput(FILE *in, FILE *out, int nrespairs, int maxchain)
   ------------------------------------------------------------------
   Handle the input file - extract the PDB filenames and process each 
   in turn

   01.04.09 Original   By: ACRM
   06.04.09 Handles maxchain
*/
BOOL HandleInput(FILE *in, FILE *out, int nrespairs, int maxchain)
{
   char filename[MAXBUFF];
   
   PDB  *pdb;
   BOOL InLoop;
   
   InLoop = TRUE;
   
   while(fgets(filename,MAXBUFF,in))
   {
      TERMINATE(filename);

      fprintf(stderr,"INFO: Processing file: %s\n",filename);
      
      if((pdb = GetPDB(filename, maxchain))==NULL)
      {
         fprintf(stderr,"WARNING: Unable to read structure\n");
      }
      else
      {
         if(!ProcessPDB(pdb, nrespairs))
         {
            fprintf(stderr,"ERROR: Unable to store all residue label pairs. \
Use -m to increase maxres\n");
            return(FALSE);
         }
         
         FREELIST(pdb, PDB);
      }
   }
   return(TRUE);
}


/************************************************************************/
/*>PDB *GetPDB(char *filename, int maxchain)
   -----------------------------------------
   Returns a PDB linked list for a PDB file of specified name. The 
   resulting linked list contains only the CAs of the first maxchain 
   chains

   01.04.09 Original    By: ACRM
   03.04.09 Rejects light or heavy chain dimers
   06.04.09 maxchain now a parameter
*/
PDB *GetPDB(char *filename, int maxchain)
{
   PDB  *pdb, *capdb, *prev, *p;
   FILE *fp;
   int  natoms;
   char *sel[2];
   int  nchain, lastresnum;
   char lastchain;
   
   if((fp=fopen(filename,"r"))==NULL)
   {
      return(NULL);
   }
   
   if((pdb=ReadPDBAtoms(fp,&natoms))==NULL)
   {
      fclose(fp);
      return(NULL);
   }

   fclose(fp);

   /* Select out only the CAs                                           */
   SELECT(sel[0],"CA  ");
   if(sel[0] == NULL)
      return(NULL);

   /* capdb will be NULL if routine fails                               */
   capdb = SelectAtomsPDB(pdb, 1, sel, &natoms);
   FREELIST(pdb, PDB);
   if(capdb == NULL)
   {
      return(NULL);
   }

   /* See if there are two identically labelled chains following eachother
      as in a light chain dimer
   */
   lastresnum = capdb->resnum;
   lastchain  = capdb->chain[0];
   for(p=capdb; p!=NULL; NEXT(p))
   {
      if((p->resnum < lastresnum) && (p->chain[0] == lastchain))
      {
         fprintf(stderr,"WARNING: File discarded as it has two \
identically labelled chains\n");
         FREELIST(capdb, PDB);
         return(NULL);
      }
      
      lastresnum = p->resnum;
      lastchain  = p->chain[0];
   }
   

   /* If the file contains more than one pair of light and heavy chains,
      terminate after the first pair
   */
   nchain = 1;
   lastchain = capdb->chain[0];
   prev = NULL;
   for(p=capdb; p!=NULL; NEXT(p))
   {
      if(p->chain[0] != lastchain)
      {
         if(++nchain > maxchain)
         {
            fprintf(stderr, "WARNING: Discarded additional chains\n");
            if(prev!=NULL)
            {
               prev->next = NULL;
            }
            FREELIST(p, PDB);
            return(capdb);
         }
         lastchain = p->chain[0];
      }
      prev = p;
   }

   return(capdb);
}


/************************************************************************/
/*>BOOL AllocateMemory(int nrespairs)
   ----------------------------------
   Allocates global memory depending on the number of residue pairs which
   have been requested.

   01.04.09 Original    By: ACRM
*/
BOOL AllocateMemory(int nrespairs)
{
   int i;
   
   /* Allocate space for storing the `constraints' and pointers         */
   if((gSx=(REAL *)malloc(nrespairs * sizeof(REAL)))==NULL)
      return(FALSE);
   for(i=0; i<nrespairs; i++)
      gSx[i] = (REAL)0.0;
   
   if((gSxSq=(REAL *)malloc(nrespairs * sizeof(REAL)))==NULL)
      return(FALSE);
   for(i=0; i<nrespairs; i++)
      gSxSq[i] = (REAL)0.0;

   if((gNVal=(int *)malloc(nrespairs * sizeof(int)))==NULL)
      return(FALSE);
   for(i=0; i<nrespairs; i++)
      gNVal[i] = (REAL)0.0;

   if((gStrings=(char **)malloc(nrespairs * sizeof(char *)))==NULL)
      return(FALSE);
   for(i=0; i<nrespairs; i++)
   {
      gStrings[i] = NULL;
   }
   for(i=0; i<nrespairs; i++)
   {
      /* Allocate strings for each                                      */
      if((gStrings[i]=(char *)malloc(MAXLABEL * sizeof(char)))==NULL)
         return(FALSE);
   }

   return(TRUE);
}


/************************************************************************/
/*>void FreeMemory(void)
   ---------------------
   Free all globally allocated memory

   01.04.09 Original    By: ACRM
*/
void FreeMemory(int nrespairs)
{
   int i;
   
   if(gSx!=NULL)   free(gSx);
   if(gSxSq!=NULL) free(gSxSq);
   if(gNVal!=NULL) free(gNVal);
   if(gStrings!=NULL)
   {
      for(i=0; i<nrespairs; i++)
      {
         if(gStrings[i]!=NULL) free(gStrings[i]);
         gStrings[i] = NULL;
      }
      free(gStrings);
   }
   gSx=NULL;
   gSxSq=NULL;
   gNVal=NULL;
   gStrings=NULL;
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
*/
void Usage(void)
{
   fprintf(stderr,"\nDistMat V1.1 (c) 2009, Dr. Andrew C.R. Martin, \
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
