/*************************************************************************

   Program:    nullstrip
   File:       nullstrip.c
   
   Version:    V1.0
   Date:       13.11.96
   Function:   Strips atoms with NULL coordinates
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1996
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)372 275775
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
   V1.0  13.11.96 Original    By: ACRM

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
PDB *StripNulls(PDB *pdb);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for stripping NULL coordinates

   03.11.96 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   PDB  *pdb;
   int  natoms;
   
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = ReadPDB(in,&natoms)) != NULL)
         {
            pdb = StripNulls(pdb);
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   13.11.96 Original    By: ACRM
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

/************************************************************************/
/*>PDB *StripNulls(PDB *pdb)
   -------------------------
   Strip NULL coordinate atoms from a PDB linked list

   13.11.96 Original   By: ACRM
*/
PDB *StripNulls(PDB *pdb)
{
   PDB *p, *prev;
   
   /* Remove from start                                                 */
   for(p=pdb; p!=NULL;)
   {
      if((p->x == (REAL)9999.0) &&
         (p->y == (REAL)9999.0) &&
         (p->z == (REAL)9999.0))
      {
         pdb = p->next;
         free(p);
      }
      else
      {
         break;
      }
      p=pdb;
   }

   /* Remove rest                                                       */
   prev = pdb;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->x == (REAL)9999.0) &&
         (p->y == (REAL)9999.0) &&
         (p->z == (REAL)9999.0))
      {
         prev->next = p->next;
         free(p);
         p=prev;
      }
      prev = p;
   }
   
   return(pdb);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   13.11.96 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nNullStrip V1.0 (c) 1996, Dr. Andrew C.R. Martin, \
UCL\n");

   fprintf(stderr,"\nUsage: nullstrip [in.pdb [out.pdb]]\n");

   fprintf(stderr,"\nRemoves atoms from a PDB file which have NULL \
coordinates (i.e.\n");
   fprintf(stderr,"x = y = z = 9999.0)\n");
}
