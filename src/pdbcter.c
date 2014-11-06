/*************************************************************************

   Program:    PDBCTer
   File:       pdbcter.c
   
   Version:    V1.0
   Date:       24.08.94
   Function:   Set naming for c-terminal oxygens and generate
               coordinates if required.
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0372) 275775
   EMail:      INTERNET: amartin@scitec.adsp.sub.org
                         martin@bsm.bioc.ucl.ac.uk
               UUCP:     ...{uunet|rutgers}!cbmehq!cbmuk!scitec!amartin
               JANET:    martin@uk.ac.ucl.bioc.bsm
               
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
#define MAXBUFF      160
#define STYLE_STD    0
#define STYLE_GROMOS 1
#define STYLE_CHARMM 2

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  int *style);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for hydrogen addition

   23.08.94 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   PDB  *pdb;
   int  natoms,
        style = STYLE_STD;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &style))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = ReadPDB(in,&natoms)) != NULL)
         {
            FixCterPDB(pdb, style);
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
      Usage();
   }
   
   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     int *style)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            int    *style       STYLE_STD, STYLE_GROMOS or STYLE_CHARMM
   Returns: BOOL                Success?

   Parse the command line
   
   23.08.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  int *style)
{
   argc--;
   argv++;

   infile[0]  = outfile[0] = '\0';
   *style = STYLE_STD;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'g':
            *style = STYLE_GROMOS;
            break;
         case 'c':
            *style = STYLE_CHARMM;
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
/*>void Usage(void)
   ----------------
   Prints a usage message

   23.08.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nPDBCTer V1.0 (c) 1994, Andrew C.R. Martin, UCL\n\n");
   fprintf(stderr,"Usage: pdbcter [-g] [-c] [<in.pdb>] [<out.pdb>]\n");
   fprintf(stderr,"               -g Gromos style C-terminii\n");
   fprintf(stderr,"               -c Charmm style C-terminii\n\n");
   fprintf(stderr,"Rename C-terminal oxygens in required style and \
generate second one if required.\n\n");
}

