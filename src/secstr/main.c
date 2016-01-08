/*************************************************************************

   Program:    ss
   File:       main.c
   
   Version:    V1.0
   Date:       19.05.99
   Function:   Main program for secondary structure calculation
   
   Copyright:  (c) Inpharmatica, Ltd. 1999
   Author:     Dr. Andrew C.R. Martin
   Address:    60 Charlotte Street,
               London W1P 2AX
   Phone:      +44 (0) 171 631 4644
   EMail:      cvs@inpharmatica.co.uk
               
**************************************************************************
 
   CVS Tags:
   =========
 
   Last modified by:    $Author: andrew $
   Tag:                 $Name:  $
   Revision:            $Revision: 1.3 $
 
**************************************************************************

   This program is copyright. Any copying without permission is illegal.

**************************************************************************

   Description:
   ============


**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   19.05.99 Original   By: ACRM

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "secstr.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
   
/************************************************************************/
/* Globals
*/
BOOL sNewChain = FALSE;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *verbose, BOOL *summary);
void Usage(void);
void WriteSummary(FILE *out, PDB *pdbStart, PDB *pdbStop);

void WriteResults(FILE *out, PDB *pdbStart, PDB *pdbStop);



/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   19.05.99 Original   By: ACRM
   27.05.99 Added error return if blCalcSS out of memory
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   FILE *in = stdin,
        *out = stdout;
   PDB  *pdb;
   int  natoms;
   BOOL verbose = FALSE,
        summary = FALSE;
   
   
   if(!ParseCmdLine(argc, argv, infile, outfile, &verbose, &summary))
   {
      Usage();
      return(0);
   }
   else
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDBAtoms(in, &natoms))!=NULL)
         {
            PDB *start, *stop;
            
            for(start=pdb; start!=NULL; start=stop)
            {
               stop=blFindNextChain(start);

               if(blCalcSecStrucPDB(start, stop, verbose) != 0)
               {
                  return(1);
               }
            
               if(summary)
               {
                  WriteSummary(out, start, stop);
               }
               else
               {
                  WriteResults(out, start, stop);
               }
            }
            
            FREELIST(pdb, PDB);
         }
      }
      else
      {
         return(1);
      }
   }
         
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *verbose, BOOL *summary)
   ---------------------------------------------------------------------
   Input:   int    argc              Argument count
            char   **argv            Argument array
   Output:  char   *infile           Input filename (or blank string)
            char   *outfile          Output filename (or blank string)
            BOOL   *verbose            Verbose?
            BOOL   *summary          Display summary rather than XMAS?
   Returns: BOOL                     Success

   Parse the command line

   19.05.99 Original    By: ACRM
   21.05.99 Added summary
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *verbose, BOOL *summary)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   *verbose = *summary = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            *verbose = TRUE;
            break;
         case 's':
            *summary = TRUE;
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
   19.05.99 Original   By: ACRM
   21.05.99 Added flags
*/
void Usage(void)
{
   fprintf(stderr,"\nss V1.0 (c) 1999, Inpharmatica, Ltd.\n");

   fprintf(stderr,"\nUsage: ss [-v][-s] [in.xmas [out.xmas]]\n");
   fprintf(stderr,"          -v Verbose mode - report dropped \
3rd Hbonds, etc.\n");
   fprintf(stderr,"          -s Summary - write a simple summary file \
rather than XMAS\n");

   fprintf(stderr,"\nCalculates secondary structure assignments \
according to the method of\n");
   fprintf(stderr,"Kabsch and Sander. Reads and writes XMAS format \
files. Input/output is\n");
   fprintf(stderr,"to standard input/output if files are not \
specified.\n\n");
}


/************************************************************************/
/*>void WriteSummary(FILE *out, PDB *pdb)
   --------------------------------------
   Write a summary file with the residue names and secondary structure

   21.05.99 Original   By: ACRM
*/
void WriteSummary(FILE *out, PDB *pdbStart, PDB *pdbStop)
{
   PDB *p;
   
   for(p=pdbStart; p!=pdbStop; p=blFindNextResidue(p))
   {
      fprintf(out, "%s%d%s %s %c\n",
              p->chain,
              p->resnum,
              p->insert,
              p->resnam,
              p->secstr);
   }
}


void WriteResults(FILE *out, PDB *pdbStart, PDB *pdbStop)
{
   WriteSummary(out, pdbStart, pdbStop);
}

