/************************************************************************/
/**

   \file       transpdb.c
   
   \version    V1.1
   \date       22.07.14
   \brief      Simple program to translate PDB files
   
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
-  V1.0  14.09.95 Original
-  V1.1  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

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
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *x, REAL *y, REAL *z);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for PDB rotation

-  17.06.94 Original    By: ACRM
-  21.07.95 Added -m
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE    *in      = stdin,
           *out     = stdout;
   PDB     *pdb;
   int     natoms;
   VEC3F   TVec;
   char    infile[MAXBUFF],
           outfile[MAXBUFF];

   TVec.x = TVec.y = TVec.z = (REAL)0.0;

   if(ParseCmdLine(argc, argv, infile, outfile, 
                   &(TVec.x), &(TVec.y), &(TVec.z)))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natoms))!=NULL)
         {
            blTranslatePDB(pdb, TVec);
            blWritePDB(out, pdb);
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

-  14.09.95 Original    By: ACRM
-  22.07.14 V1.1 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\nTransPDB V1.1  (c) 1995-2014 Andrew C.R. Martin\n");
   fprintf(stderr,"Freely distributable if no profit is made\n\n");
   fprintf(stderr,"Usage: transpdb [-x <x>] [-y <y>] [-z <z>] [-h]\n");
   fprintf(stderr,"              [<input.pdb> [<output.pdb>]]\n");
   fprintf(stderr,"I/O is to stdin/stdout if not specified\n\n");
   fprintf(stderr,"Translates a PDB file\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     REAL *x, REAL *y, REAL *z)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *x           X-translation
   \param[out]     *y           Y-translation
   \param[out]     *z           Z-translation
   \return                     Success?

   Parse the command line
   
-  05.07.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *x, REAL *y, REAL *z)
{
   REAL temp;
   
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'x':
            argc--;
            argv++;
            if(!sscanf(argv[0],"%lf",&temp))
               return(FALSE);
            *x += temp;
            break;
         case 'y':
            argc--;
            argv++;
            if(!sscanf(argv[0],"%lf",&temp))
               return(FALSE);
            *y += temp;
            break;
         case 'z':
            argc--;
            argv++;
            if(!sscanf(argv[0],"%lf",&temp))
               return(FALSE);
            *z += temp;
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

