/*************************************************************************

   Program:    PDBHAdd
   File:       pdbhadd.c
   
   Version:    V1.1
   Date:       24.11.95
   Function:   Add hydrogens to a PDB file
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
   EMail:      INTERNET: martin@biochem.ucl.ac.uk
               
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  24.08.94 Original
   V1.1  24.11.95 Changes NT atoms to N

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *pgpfile, BOOL *AllH, BOOL *Charmm);
void FixNTerNames(PDB *pdb);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for hydrogen addition

   23.08.94 Original    By: ACRM
   24.11.95 Added call to FixNTerNames()

*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout,
        *pgp = NULL;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        pgpfile[MAXBUFF];
   PDB  *pdb;
   int  natoms, nhyd;
   BOOL AllH   = FALSE,
        Charmm = FALSE;
   
   if(ParseCmdLine(argc, argv, infile, outfile, pgpfile, &AllH, &Charmm))
   {
      if((pgp = OpenPGPFile(pgpfile, AllH)) != NULL)
      {
         if(OpenStdFiles(infile, outfile, &in, &out))
         {
            if((pdb = ReadPDB(in,&natoms)) != NULL)
            {
               FixNTerNames(pdb);
               
               nhyd = HAddPDB(pgp, pdb);
               nhyd += AddNTerHs(&pdb, Charmm);
               
               fprintf(stderr,"%d hydrogens were added.\n",nhyd);

               RenumAtomsPDB(pdb);
               
               WritePDB(out, pdb);
            }
            else
            {
               fprintf(stderr,"No atoms read from PDB file\n");
            }
         }
      }
      else
      {
         fprintf(stderr,"Error: Unable to open proton generation \
parameter file.\n");
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
                     char *pgpfile, BOOL *AllH, BOOL *Charmm)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            char   *pgpfile     PGP filename
            BOOL   *AllH        Add all hydrogens
            BOOL   *Charmm      Do Charmm style N-terminii
   Returns: BOOL                Success?

   Parse the command line
   
   23.08.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *pgpfile, BOOL *AllH, BOOL *Charmm)
{
   argc--;
   argv++;

   infile[0]  = outfile[0] = '\0';
   pgpfile[0]  = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'p':
            argc--;
            argv++;
            strncpy(pgpfile,argv[0],MAXBUFF);
            pgpfile[MAXBUFF-1] = '\0';
            break;
         case 'a':
            *AllH = TRUE;
            break;
         case 'c':
            *Charmm = TRUE;
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
/*>void FixNTerNames(PDB *pdb)
   ---------------------------
   Change atom name NT to N. Fixes problem with Nter H addition

   24.11.95 Original   By: ACRM
*/
void FixNTerNames(PDB *pdb)
{
   PDB *p;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"NT  ",4))
         strcpy(p->atnam,"N   ");
   }
}



/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   23.08.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nPDBHAdd V1.1 (c) 1994, Andrew C.R. Martin, UCL\n\n");
   fprintf(stderr,"Usage: pdbhadd [-p pgpfile] [-a] [-c] [<in.pdb> \
[<out.pdb>]]\n");
   fprintf(stderr,"               -p Specify proton generation \
parameter file\n");
   fprintf(stderr,"               -a Add ALL hydrogens.\n");
   fprintf(stderr,"               -c Do Charmm style N-terminii.\n");
   fprintf(stderr,"\nAdd hydrogens to a PDP file.\n");
   fprintf(stderr,"Note that you should first strip any existing \
hydrogens with hstrip.\n\n");
}




