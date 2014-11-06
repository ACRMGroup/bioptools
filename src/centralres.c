/*************************************************************************

   Program:    centralres
   File:       centralres.c
   
   Version:    V1.0
   Date:       27.07.12
   Function:   Find the residue nearest the centroid of a protein
   
   Copyright:  (c) Dr. Andrew C. R. Martin 2012
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
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
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"

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
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------

   27.07.12 Original   By: ACRM
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
            VEC3F cg;
            PDB   *p, *closest;
            REAL  bestDistSq = 99999.0,
               dist;
            
            GetCofGPDB(pdb, &cg);
            for(p=pdb; p!=NULL; NEXT(p))
            {
               if(!strncmp(p->atnam, "CA  ",4))
               {
                  dist = DISTSQ(p, (&cg));
                  if(dist < bestDistSq)
                  {
                     bestDistSq = dist;
                     closest = p;
                  }
               }
            }
            fprintf(out, "%c%d%c\n", 
                    closest->chain[0], closest->resnum, closest->insert[0]);
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
   
   27.07.12 Original    By: ACRM
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
/*>void Usage(void)
   ----------------
   Prints a usage message

   08.11.96 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\ncentralres V1.0 (c) 2012 UCL, \
Dr. Andrew C.R. Martin\n");
   fprintf(stderr,"\nUsage: centralres [in.pdb [out.pdb]]\n");

   fprintf(stderr,"Identifies the residue closest to the centroid \
of a protein.\n\n");
}


