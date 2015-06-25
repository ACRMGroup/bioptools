/************************************************************************/
/**

   \file       pdbml2pdb.c
   
   \version    V1.0
   \date       25.06.15
   \brief      Convert PDBML format to PDB
   
   \copyright  (c) Dr. Andrew C. R. Martin 2015
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
-  V1.0  25.06.15 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 256

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program

-  25.06.15 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE     *in      = stdin,
            *out     = stdout;
   WHOLEPDB *wpdb;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];

   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb=blReadWholePDB(in))!=NULL)
         {
            FORCEPDB;
            blWriteWholePDB(out, wpdb);
         }
         else
         {
            fprintf(stderr,"No atoms read from PDB file\n");
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
/*>void Usage(void)
   ----------------
*//**
   Prints a usage message

-  26.02.15 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\npdbml2pdb V1.0  (c) 2015 UCL, Andrew C.R. \
Martin\n");
   fprintf(stderr,"Usage: pdbml2pdb [<input.pdb> [<output.pdb>]]\n");
   fprintf(stderr,"I/O is to stdin/stdout if not specified\n\n");
   fprintf(stderr,"Converts a PDBML file to PDB format.\n");
   fprintf(stderr,"The file will only be converted if it is suitable \
for PDB format\n");
   fprintf(stderr,"\nNote that only key header information is retained \
at present\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \return                      Success?

   Parse the command line
   
-  25.06.15 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
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

