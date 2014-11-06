/*************************************************************************

   Program:    splitchains
   File:       splitchains.c
   
   Version:    V1.1
   Date:       16.01.14
   Function:   Split a PDB file into separate chains
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 1997-2014
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
   V1.0    10.07.97  Original   By: ACRM
   V1.1    16.01.14  Fixed to handle HETATMs properly

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF     160
#define MAXCHAINID    8
#define MAXCHAINS  1024
   

/************************************************************************/
/* Globals
*/
static BOOL gQuiet = FALSE;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, BOOL *current);
BOOL WriteEachPDBChain(char *InFile, PDB *pdb, BOOL current);
void Usage(void);
BOOL BuildFileName(char *OutFile, int maxFileName, char *InFile, 
                   char *chain, BOOL current);
char *mystrncat(char *out, const char *in, size_t len);



/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for extracting selected chains from a PDB file

   07.02.97 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char InFile[MAXBUFF];
   FILE *in  = stdin;
   PDB  *pdb;
   int  natoms;
   BOOL current;
   
   
   if(ParseCmdLine(argc, argv, InFile, &current))
   {
      if(OpenStdFiles(InFile, NULL, &in, NULL))
      {
         if((pdb=ReadPDB(in, &natoms))==NULL)
         {
            if(!gQuiet)
               fprintf(stderr,"No atoms read from input PDB file\n");
            return(1);
         }
         else
         {
            if(!WriteEachPDBChain(InFile,pdb,current))
            {
               if(!gQuiet)
                  fprintf(stderr,"splitchains: Failed to write all \
output files\n");
            }
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, BOOL *current)
   ---------------------------------------------------------------------
   Input:   int    argc        Argument count
            char   **argv      Argument array
   Output:  char   *infile     Input filename (or blank string)
   Returns: BOOL               Success

   Parse the command line

   10.07.97 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, BOOL *current)
{
   argc--;
   argv++;
   
   infile[0] = '\0';
   *current = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         case 'q':
            gQuiet = TRUE;
            break;
         case 'c':
            *current = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 0-1 arguments left                     */
         if(argc > 1)
            return(FALSE);
         if(argc)
         {
            strcpy(infile, argv[0]);
         }

         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL WriteEachPDBChain(char *InFile, PDB *pdb, BOOL current)
   ------------------------------------------------------------
   Input:   char   *InFile     Input filename
            PDB    *pdb        PDB linked list
            BOOL   current     Strip path and write to current directory

   Writes each chain to a separate file

   16.01.14 Rewritten to deal with HETATMs properly   By: ACRM
*/
BOOL WriteEachPDBChain(char *InFile, PDB *pdb, BOOL current)
{
   PDB *p;
   char chain[MAXCHAINID];
   char OutFile[MAXBUFF];
   FILE *fp;
   char chains[MAXCHAINS][MAXCHAINID];
   int  nchains = 0,
        i;

   /* Build a list of the observed chains                               */
   chain[0] = '\0'; /* Dummy initialization                             */
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(strncmp(p->chain, chain, MAXCHAINID))
      {
         /* Chain has changed, has it been seen before?                 */
         BOOL found = FALSE;
         for(i=0; i<nchains; i++)
         {
            if(!strncmp(p->chain, chains[i], MAXCHAINID))
            {
               found=TRUE;
               break;
            }
         }

         if(!found)
         {
            strncpy(chains[nchains], p->chain, MAXCHAINID);
            nchains++;
         }
         
         strncpy(chain, pdb->chain, MAXCHAINID);
      }
   }

   /* Go through each chain and write the residues that match that chain*/
   for(i=0; i<nchains; i++)
   {
      /* Open a file to store the previous chain                        */
      if(BuildFileName(OutFile, MAXBUFF, InFile, chains[i], current))
      {
         if((fp=fopen(OutFile, "w"))==NULL)
         {
            if(!gQuiet)
               fprintf(stderr,"splitchains: Could not write output \
file: %s\n", OutFile);
            return(FALSE);
         }
      }
      else
      {
         if(!gQuiet)
            fprintf(stderr,"splitchains: No memory to build output \
filename\n");
         return(FALSE);
      }

      /* Opened output file OK, so write the records to it              */
      for(p=pdb; p!=NULL; NEXT(p))
      {
         /* If the chain has changed                                    */
         if(!strncmp(p->chain, chains[i], MAXCHAINID))
         {
            WritePDBRecord(fp,p);
         }
      }
      fprintf(fp, "TER   \n");
      fclose(fp);
   }
   
            
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   10.07.97 Original
   16.01.14 V1.1
*/
void Usage(void)
{
   fprintf(stderr,"SplitChains V1.1 (c) 1997-2014 \
Dr. Andrew C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: splitchains [-c][-q] [in.pdb]\n");
   fprintf(stderr,"       -c  Output to current directory\n");
   fprintf(stderr,"       -q  Quiet - no error messages\n");

   fprintf(stderr,"\nSplitChains takes a PDB file (of the specified name \
or from stdin if a\n");
   fprintf(stderr,"filename is not given) and creates separate output \
files for each\n");
   fprintf(stderr,"chain.\n");

   fprintf(stderr,"\nIf a filename is specified, the output files will \
be will be given the\n");
   fprintf(stderr,"basename of the file followed by the chain name and \
a .pdb extension.\n");
   fprintf(stderr,"For example, if the file is given as pdb3hfl.ent, the \
output files would\n");
   fprintf(stderr,"be pdb3hflL.pdb, pdb3hflH.pdb, pdb3hflY.pdb (3hfl has \
L,H and Y chains).\n");

   fprintf(stderr,"\nIf the -c flag is given, any path specified for the \
file will be removed\n");
   fprintf(stderr,"before the output filename is created such that files \
are written to the\n");
   fprintf(stderr,"current directory.\n");

   fprintf(stderr,"\nIf no filename is given (input comes from stdin), \
output will simply\n");
   fprintf(stderr,"be the chain name with the .pdb extension and will be \
placed in the\n");
   fprintf(stderr,"current directory.\n");

   fprintf(stderr,"\nA blank chain name is converted to the digit 0\n\n");
}


/************************************************************************/
/*>BOOL BuildFileName(char *OutFile, int maxOutFileName, char *InFile, 
                      char *chain, BOOL current)
   -------------------------------------------------------------------
   Input:    int   maxOutFileName    Max buffer size for output filename
             char  *InFile           Input filename
             char  *chain            Chain name to put into filename
             BOOL  current           Strip path info so we write to 
                                     current directory
   Output:   char  *OutFile          Output filename
   Returns:  BOOL                    success?

   Constructs a filename from a current filename plus a chain name

   10.07.97 Original
   16.01.14 Uses mystrncat() and chain is handled as a string
*/
BOOL BuildFileName(char *OutFile, int maxOutFileName, char *InFile, 
                   char *chain, BOOL current)
{
   char *chp,
        *chpt,
        *WorkFile,
        *stem,
        *path;

   if(chain[0] == ' ')
      strcpy(chain, "0");

   if(InFile[0])
   {
      /* Make a working copy of the input filename                      */
      if((WorkFile=(char *)malloc((1+strlen(InFile))*sizeof(char)))==NULL)
         return(FALSE);
      strcpy(WorkFile,InFile);
      
      chp = WorkFile;
      
      /* Find the start of the file stem (after all directory specs)    */
      while((chpt=strchr(chp, '/'))!=NULL)
         chp = chpt+1;
      stem = chp;
      
      /* Terminate at the end of the path                               */
      if((chp!=WorkFile) && (*(chp-1) == '/'))
      {
         path = WorkFile;
         *(chp-1) = '\0';
      }
      else
      {
         path = NULL;
      }
      
      /* Remove any extension                                           */
      chp = stem + strlen(stem);
      while(chp > stem)
      {
         if(*chp == '.')
         {
            *chp = '\0';
            break;
         }
         chp--;
      }
      
      /* Build the final filename                                       */
      OutFile[0] = '\0';
      if(!current && path!=NULL)         /* The path                    */
      {
         mystrncat(OutFile, path, maxOutFileName);
         mystrncat(OutFile, "/",  maxOutFileName);
      }
      mystrncat(OutFile, stem, maxOutFileName);  /* The filestem        */
      mystrncat(OutFile, chain, maxOutFileName); /* The chain name      */
      mystrncat(OutFile,".pdb", maxOutFileName); /* The extension       */
      
      /* Free allocated memory                                          */
      free(WorkFile);
   }
   else
   {
      strcpy(OutFile, chain);
      mystrncat(OutFile,".pdb", maxOutFileName);
   }

   return(TRUE);
}


/************************************************************************/
/*>char *mystrncat(char *out, const char *in, size_t len)
   ------------------------------------------------------
   Input:   const char   *in     Input string to be appended
            size_t       len     Length of output string
   I/O:     char         *out    String that we are appending to

   A simpler version of strncat. strncat takes the max number of chars
   to be appended whereas this takes the max number of chars that 'out'
   can hold.

   16.01.14  Original   By: ACRM
*/
char *mystrncat(char *out, const char *in, size_t len)
{
   int lenOut, lenIn, cpLen;
   
   lenOut = strlen(out);
   lenIn  = strlen(in);
   cpLen  = len - lenOut - lenIn;
   
   strncat(out, in, cpLen);
   return(out);
}
