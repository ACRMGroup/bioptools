/*************************************************************************

   Program:    hetstrip
   File:       hetstrip.c
   
   Version:    V1.0
   Date:       17.03.95
   Function:   Strip het atoms from a PDB file. Acts as filter
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0372) 275775
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
   V1.0  17.03.95 Original

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
   Main program for stripping hydrogens

   04.07.94 Original    By: ACRM
   15.07.94 Now writes TER cards and returns 0 correctly
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
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=ReadPDBAtoms(in,&natom))!=NULL)
         {
            WritePDB(out,pdb);
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
   Input:   int    argc        Argument count
            char   **argv      Argument array
   Output:  char   *infile     Input filename (or blank string)
            char   *outfile    Output filename (or blank string)
   Returns: BOOL               Success

   Parse the command line

   16.08.94 Original    By: ACRM
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
   Prints a usage message

   17.03.95 Original    By: ACRM
*/
void Usage(void)
{            
   fprintf(stderr,"\nHetStrip V1.0 (c) 1994, Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: hetstrip [<in.pdb>] [<out.pdb>]\n\n");
   fprintf(stderr,"Removes het atoms from a PDB file. I/O is through \
stdin/stdout if files\n");
   fprintf(stderr,"are not specified.\n\n");
}
