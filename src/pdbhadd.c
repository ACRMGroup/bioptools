/************************************************************************/
/**

   \file       pdbhadd.c
   
   \version    V1.2
   \date       22.07.14
   \brief      Add hydrogens to a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin 1994-2014
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  24.08.94 Original
-  V1.1  24.11.95 Changes NT atoms to N
-  V1.2  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

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
*//**

   Main program for hydrogen addition

-  23.08.94 Original    By: ACRM
-  24.11.95 Added call to FixNTerNames()

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
      if((pgp = blOpenPGPFile(pgpfile, AllH)) != NULL)
      {
         if(blOpenStdFiles(infile, outfile, &in, &out))
         {
            if((pdb = blReadPDB(in,&natoms)) != NULL)
            {
               FixNTerNames(pdb);
               
               nhyd = blHAddPDB(pgp, pdb);
               nhyd += blAddNTerHs(&pdb, Charmm);
               
               fprintf(stderr,"%d hydrogens were added.\n",nhyd);

               blRenumAtomsPDB(pdb);
               
               blWritePDB(out, pdb);
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
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *pgpfile     PGP filename
   \param[out]     *AllH        Add all hydrogens
   \param[out]     *Charmm      Do Charmm style N-terminii
   \return                     Success?

   Parse the command line
   
-  23.08.94 Original    By: ACRM
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
*//**

   Change atom name NT to N. Fixes problem with Nter H addition

-  24.11.95 Original   By: ACRM
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
*//**

   Prints a usage message

-  23.08.94 Original    By: ACRM
-  22.07.14 V1.2 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\nPDBHAdd V1.2 (c) 1994-2014, Andrew C.R. Martin, UCL\n\n");
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




