/************************************************************************/
/**

   \file       checkpdb.c
   
   \version    V1.0
   \date       28.01.18
   \brief      Check a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin 1996-2018
   \author     Dr. Andrew C. R. Martin
   \par
               Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This code is NOT IN THE PUBLIC DOMAIN, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC.

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified.

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   Note:

   Need to deal with 3zeu where residue A41 and A194 is both MET and MSE!


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

-  V1.0  16.08.18 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>

#include "bioplib/pdb.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF       512
#define MAXRESSPEC     16
#define MAXBONDDISTSQ   4.0    /* Distance = 2.0A - typical is 1.32A    */

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL CheckBackboneAtoms(FILE *out, PDBSTRUCT *pdbs, BOOL verbose);
BOOL CheckBackboneContinuity(FILE *out, PDBSTRUCT *pdbs, BOOL verbose);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *verbose);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program for checking a PDB file - currently just does backbone
   continuity

-  16.08.18 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE    *in  = stdin,
           *out = stdout;
   char    inFile[MAXBUFF],
           outFile[MAXBUFF];
   int     natoms  = 0,
           retval  = 0;
   PDB     *pdb    = NULL;
   BOOL    verbose = FALSE;

   if(ParseCmdLine(argc, argv, inFile, outFile, &verbose))
   {
      if(blOpenStdFiles(inFile, outFile, &in, &out))
      {
         if((pdb=blReadPDBAtoms(in, &natoms))!=NULL)
         {
            PDBSTRUCT *pdbs = NULL;
            if((pdbs = blAllocPDBStructure(pdb))!=NULL)
            {
               if(!CheckBackboneAtoms(out, pdbs, verbose))
               {
                  fprintf(out, "BAD\n");
                  retval = 1;
               }
               else
               {
                  if(CheckBackboneContinuity(out, pdbs, verbose))
                  {
                     fprintf(out, "OK\n");
                     retval = 0;
                  }
                  else
                  {
                     fprintf(out, "BAD\n");
                     retval = 1;
                  }
               }
            }
            else
            {
               fprintf(stderr,"checkpdb: Error - unable to allocate \
PDB structure\n");
               return(1);
            }
         }
         else
         {
            fprintf(stderr,"checkpdb: Error - no atoms read from PDB \
file\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"checkpdb: Error - unable to open input or \
output file\n");
         return(1);
      }
   }
   else
   {
      Usage();
   }
   
   return(retval);
}


/************************************************************************/
/*>BOOL CheckBackboneAtoms(FILE *out, PDBSTRUCT *pdbs, BOOL verbose)
   -----------------------------------------------------------------
*//**
   \param[in]   *out    Output file pointer
   \param[in]   *pdbs   Structured PDB data
   \param[in]   verbose Whether to print information
   \return              PDB OK?

   Checks whether all backbone atoms are present thoughout the PDB file

-  16.08.18 Original   By: ACRM
*/
BOOL CheckBackboneAtoms(FILE *out, PDBSTRUCT *pdbs, BOOL verbose)
{
   PDBCHAIN *pdbchain = NULL;
   BOOL     retval    = TRUE;

   /* For each chain                                                    */
   for(pdbchain=pdbs->chains; pdbchain!=NULL; NEXT(pdbchain))
   {
      PDBRESIDUE *pdbres = NULL;

      /* For each residue in this chain                                 */
      for(pdbres=pdbchain->residues; pdbres!=NULL; NEXT(pdbres))
      {
         PDB  *start = pdbres->start,
              *stop  = pdbres->stop,
              *p;
         BOOL gotN   = FALSE,
              gotCA  = FALSE,
              gotC   = FALSE,
              gotO   = FALSE;
         
         /* For each atom in this residue, check we have the backbone
          atoms
         */
         for(p=start; p!=stop; NEXT(p))
         {
            if(!strncmp(p->atnam, "N   ", 4))
            {
               gotN  = TRUE;
            }
            else if(!strncmp(p->atnam, "CA  ", 4))
            {
               gotCA = TRUE;
            }
            else if(!strncmp(p->atnam, "C   ", 4))
            {
               gotC  = TRUE;
            }
            else if(!strncmp(p->atnam, "O   ", 4))
            {
               gotO  = TRUE;
            }
            else if(!strncmp(p->atnam, "OXT ", 4))
            {
               gotO  = TRUE;
            }
         }

         /* If any backbone atom is missing, set the return status and
            print an error message
         */
         if(!(gotN && gotCA && gotC && gotO))
         {
            retval = FALSE;
            if(verbose)
            {
               char resspec[MAXRESSPEC];
               blBuildResSpec(start, resspec);
               if(!gotN)
               {
                  fprintf(out, "Residue %s is missing backbone atom N\n",
                          resspec);
               }
               if(!gotCA)
               {
                  fprintf(out, "Residue %s is missing backbone atom CA\n",
                          resspec);
               }
               if(!gotC)
               {
                  fprintf(out, "Residue %s is missing backbone atom C\n",
                          resspec);
               }
               if(!gotO)
               {
                  fprintf(out, "Residue %s is missing backbone atom O\n",
                          resspec);
               }
            }
         }
      }
   }
   return(retval);
}


/************************************************************************/
/*>BOOL CheckBackboneContinuity(FILE *out, PDBSTRUCT *pdbs, BOOL verbose)
   ----------------------------------------------------------------------
*//**
   \param[in]   *out    Output file pointer
   \param[in]   *pdbs   Structured PDB data
   \param[in]   verbose Whether to print information
   \return              PDB OK?

   Checks whether backbone C of one residue is linked to N of next

-  16.08.18 Original   By: ACRM
*/
BOOL CheckBackboneContinuity(FILE *out, PDBSTRUCT *pdbs, BOOL verbose)
{
   PDBCHAIN *pdbchain = NULL;
   BOOL     retval    = TRUE;

   /* For each chain                                                    */
   for(pdbchain=pdbs->chains; pdbchain!=NULL; NEXT(pdbchain))
   {
      PDBRESIDUE *residue = NULL;
   
      /* For each residue                                               */
      for(residue=pdbchain->residues; residue!=NULL; NEXT(residue))
      {
         PDBRESIDUE *nextResidue = residue->next;

         /* If there is a following residue                             */
         if(nextResidue)
         {
            PDB *c = NULL,
                *n = NULL;
            
            /* Find the C of the first residue and the N of the second
               residue
            */
            c = blFindAtomInRes(residue->start,     "C   ");
            n = blFindAtomInRes(nextResidue->start, "N   ");
            if((c==NULL)||(n==NULL))
            {
               retval = FALSE;
            }
            else
            {
               if(DISTSQ(c,n) > MAXBONDDISTSQ)
               {
                  if(verbose)
                  {
                     char resspecC[MAXRESSPEC],
                          resspecN[MAXRESSPEC];

                     blBuildResSpec(c, resspecC);
                     blBuildResSpec(n, resspecN);
                     fprintf(out, 
                             "Residue %s is not joined to residue %s\n",
                             resspecC, resspecN);
                  }
                  retval = FALSE;
               }
            }
         }
      }
   }
   return(retval);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     BOOL *verbose)
   ---------------------------------------------------------------------
*//**
   \param[in]     argc         Argument count
   \param[in]     **argv       Argument array
   \param[out]    *infile      Input file (or blank string)
   \param[out]    *outfile     Output file (or blank string)
   \param[out]    *verbose     Verbose mode
   \return                     Success?

   Parse the command line
   
-  16.08.18 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *verbose)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            *verbose = TRUE;
            break;
         case 'h':
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are <= 2 arguments left                    */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
         }
         
         /* Copy the second to outfile                                  */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
         }
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**
   Displays a usage message

-  16.08.18 Original   By: ACRM
*/
void Usage(void)
{
   printf("\ncheckpdb V1.0 (c) 2018 UCL, Dr. Andrew C.R. Martin\n");

   printf("\nUsage: checkpdb [-v] [in.pdb [out.txt]]\n");
   printf("       -v   Verbose - prints information about errors\n");

   printf("\nThis is the start of a detailed PDB checking program. \
Currently it\n");
   printf("simply checks that all backbone atoms are present and that \
residues\n");
   printf("are all joined as they should be. If run without -v it \
simply\n");
   printf("prints OK or BAD (and returns 0 or 1 respectively). With \
-v it gives\n");
   printf("information about the errors found.\n\n");
}



