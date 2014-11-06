/*************************************************************************

   Program:    pdb2xyz
   File:       pdb2xyz.c
   
   Version:    V1.0
   Date:       24.08.94
   Function:   Convert PDB to Gromos XYZ
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0372) 275775
   EMail:      INTERNET: amartin@scitec.adsp.sub.org
                         martin@bsm.bioc.ucl.ac.uk
               UUCP:     ...{uunet|rutgers}!cbmehq!cbmuk!scitec!amartin
               JANET:    martin@uk.ac.ucl.bioc.bsm
               
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
   V1.0  23.08.94 Original   By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
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
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *title);
void Usage(void);
void WriteXYZ(FILE *out, PDB *pdb, int natoms, char *title);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for format conversion

   23.08.94 Original    By: ACRM
   24.08.94 Changed to call OpenStdFiles()
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        title[MAXBUFF];
   PDB  *pdb;
   int  natoms;
   
   if(ParseCmdLine(argc, argv, infile, outfile, title))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = ReadPDB(in,&natoms)) != NULL)
         {
            WriteXYZ(out, pdb, natoms, title);
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
                     char *title)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            char   *title       Title for file
   Returns: BOOL                Success?

   Parse the command line
   
   23.08.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *title)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   title[0]  = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 't':
            argc--;
            argv++;
            strncpy(title,argv[0],MAXBUFF);
            title[MAXBUFF-1] = '\0';
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

   23.08.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nPDB2XYZ V1.0 (c) 1994, Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: pdb2xyz [-t title] [<in.pdb>] [<out.pdb>]\n\n");
   fprintf(stderr,"Convert PDB format to GROMOS XYZ. N.B. Does NOT \
correct atom order.\n\n");
}

/************************************************************************/
/*>void WriteXYZ(FILE *out, PDB *pdb, int natoms, char *title)
   -----------------------------------------------------------
   Input:    FILE   *out      Output file pointer
             PDB    *pdb      PDB linked list
             int    natoms    Number of atoms
             char   title     Title for XYZ file or blank string

   Write the PDB linked list out in Gromos XYZ format

   23.08.94 Original    By: ACRM
*/
void WriteXYZ(FILE *out, PDB *pdb, int natoms, char *title)
{
   PDB *p;
   int i = 1;

   if(title[0])
      fprintf(out,"%s\n",title);
   else
      fprintf(out,"Gromos XYZ file generated by PDB2XYZ\n");

   fprintf(out,"%5d\n",natoms);
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      fprintf(out,"%5d %-4s %-4s %5d%8.3f%8.3f%8.3f\n",
              p->resnum,
              p->resnam,
              p->atnam,
              i++,
              p->x,
              p->y,
              p->z);
   }
}

