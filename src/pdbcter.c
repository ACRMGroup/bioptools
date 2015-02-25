/************************************************************************/
/**

   \file       pdbcter.c
   
   \version    V1.2
   \date       25.02.15
   \brief      Set naming for c-terminal oxygens and generate
               coordinates if required.
   
   \copyright  (c) Dr. Andrew C. R. Martin 1994-2015
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

-  V1.1  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.2  25.02.15 Modified for new blRenumAtomsPDB()
                  Supports whole PDB

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
*//**

   Main program for hydrogen addition

-  23.08.94 Original    By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  25.02.15 Supports whole PDB
*/
int main(int argc, char **argv)
{
   FILE     *in  = stdin,
            *out = stdout;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];
   WHOLEPDB *wpdb;
   PDB      *pdb;
   int      style = STYLE_STD;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &style))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb = blReadWholePDB(in)) != NULL)
         {
            pdb = wpdb->pdb;
            
            blFixCterPDB(pdb, style);
            blRenumAtomsPDB(pdb, 1);
               
            blWriteWholePDB(out, wpdb);
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
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *style       STYLE_STD, STYLE_GROMOS or STYLE_CHARMM
   \return                     Success?

   Parse the command line
   
-  23.08.94 Original    By: ACRM
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
*//**

   Prints a usage message

-  23.08.94 Original    By: ACRM
-  22.07.14 V1.1 By: CTP
-  24.02.15 V1.2 and improved help message By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nPDBCTer V1.2 (c) 1994-2015, Andrew C.R. Martin, \
UCL\n\n");
   fprintf(stderr,"Usage: pdbcter [-g] [-c] [in.pdb [out.pdb]]\n");
   fprintf(stderr,"               -g Gromos style C-terminii\n");
   fprintf(stderr,"               -c Charmm style C-terminii\n\n");
   fprintf(stderr,"Rename C-terminal oxygens in required style and \
generate second one\n");
   fprintf(stderr,"if required.\n");
   fprintf(stderr,"\nDefault is to name and generate second oxygen as \
OXT. Gromos names and\n");
   fprintf(stderr,"generates both terminal oxygens as O1 and O2. \
Charmm names first oxygen\n");
   fprintf(stderr,"OT1 and generates a CTER residue containing \
OT2.\n\n");
}

