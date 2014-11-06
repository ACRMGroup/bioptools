/************************************************************************/
/**

   \file       hetstrip.c
   
   \version    V1.1
   \date       22.07.14
   \brief      Strip het atoms from a PDB file. Acts as filter
   
   \copyright  (c) Dr. Andrew C. R. Martin 1995-2014
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
-  V1.0  17.03.95 Original
-  V1.1  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
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
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for stripping hydrogens

-  04.07.94 Original    By: ACRM
-  15.07.94 Now writes TER cards and returns 0 correctly
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   PDB  *pdb;
   int  natom;
   FILE *in  = stdin, 
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];

   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=blReadPDBAtoms(in,&natom))!=NULL)
         {
            blWritePDB(out,pdb);
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
   ----------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \param[out]     *outfile    Output filename (or blank string)
   \return                     Success

   Parse the command line

-  16.08.94 Original    By: ACRM
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
         case 'h':
            return(FALSE);
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

-  17.03.95 Original    By: ACRM
-  22.07.14 V1.1 By: CTP
*/
void Usage(void)
{            
   fprintf(stderr,"\nHetStrip V1.1 (c) 1994-2014, Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: hetstrip [<in.pdb>] [<out.pdb>]\n\n");
   fprintf(stderr,"Removes het atoms from a PDB file. I/O is through \
stdin/stdout if files\n");
   fprintf(stderr,"are not specified.\n\n");
}
