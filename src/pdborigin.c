/************************************************************************/
/**

   \file       pdborigin.c
   
   \version    V1.2
   \date       13.02.15
   \brief      Add hydrogens to a PDB file
   
   \copyright  (c) UCL, Dr. Andrew C. R. Martin 1999-2015
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
-  V1.0  28.01.99 Original
-  V1.1  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.2  13.02.15 Added whole PDB support  By: ACRM

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
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for hydrogen addition

-  28.01.99 Original    By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  13.02.15 Added whole PDB support  By: ACRM
*/
int main(int argc, char **argv)
{
   FILE     *in  = stdin,
            *out = stdout;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];
   WHOLEPDB *wpdb;
   PDB      *pdb;
   
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb = blReadWholePDB(in)) != NULL)
         {
            pdb=wpdb->pdb;
            blOriginPDB(pdb);
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \return                     Success?

   Parse the command line
   
-  28.01.99 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
{
   argc--;
   argv++;

   infile[0]  = outfile[0] = '\0';
   
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


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  28.01.99 Original    By: ACRM
-  22.07.14 V1.1 By: CTP
-  13.02.15 V1.2 By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\npdborigin V1.2 (c) 1999-2015, UCL, Andrew C.R. \
Martin\n\n");
   fprintf(stderr,"Usage: pdborigin [in.pdb [out.pdb]]\n");
   fprintf(stderr,"\nMoves a set of PDB coordinates such that the \
centre of geometry\n");
   fprintf(stderr,"is at the origin.\n\n");
}




