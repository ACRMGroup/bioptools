/************************************************************************/
/**

   \file       pdbconect.c
   
   \version    V1.0
   \date       26.02.15
   \brief      Rebuild CONECT records for a PDB file
   
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
-  V1.0  26.02.15 Original

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
#define DEF_TOL 0.2

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *tol);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program

-  26.02.15 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE     *in      = stdin,
            *out     = stdout;
   WHOLEPDB *wpdb;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];
   REAL     tol = DEF_TOL;

   if(ParseCmdLine(argc, argv, infile, outfile, &tol))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb=blReadWholePDB(in))!=NULL)
         {
            blBuildConectData(wpdb->pdb, tol);
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

-  26.02.15 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\npdbconect V1.0  (c) 2015 UCL, Andrew C.R. \
Martin\n");
   fprintf(stderr,"Usage: pdbconect [-t x] [<input.pdb> \
[<output.pdb>]]\n");
   fprintf(stderr,"       -t specify tolerance [Default: %.1f]\n", 
           DEF_TOL);

   fprintf(stderr,"\nGenerates CONECT records for a PDB file from the \
covalent radii of the\n");
   fprintf(stderr,"elements involved. Existing CONECT records are \
discarded first\n");
   fprintf(stderr,"I/O is to stdin/stdout if not specified\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     REAL *tol)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *tol         Covalent bond tolerance
   \return                      Success?

   Parse the command line
   
-  26.02.15 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *tol)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *tol      = DEF_TOL;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 't':
            argc--;
            argv++;
            if(argc)
            {
               if(!sscanf(argv[0], "%lf", tol))
                  return(FALSE);
            }
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

