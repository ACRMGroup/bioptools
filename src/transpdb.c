/*************************************************************************

   Program:    transpdb
   File:       transpdb.c
   
   Version:    V1.0
   Date:       14.09.95
   Function:   Simple program to translate PDB files
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0) 1372 275775
   EMail:      martin@biochem.ucl.ac.uk
               
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
   V1.0  14.09.95 Original

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
   Main program for PDB rotation

   17.06.94 Original    By: ACRM
   21.07.95 Added -m
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
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=ReadPDB(in, &natoms))!=NULL)
         {
            TranslatePDB(pdb, TVec);
            WritePDB(out, pdb);
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
   14.09.95 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nTransPDB V1.0  (c) 1995 Andrew C.R. Martin\n");
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
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            REAL   *x           X-translation
            REAL   *y           Y-translation
            REAL   *z           Z-translation
   Returns: BOOL                Success?

   Parse the command line
   
   05.07.94 Original    By: ACRM
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

