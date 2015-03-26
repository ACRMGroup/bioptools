/************************************************************************/
/**

   \file       pdbsplitchains.c
   
   \version    V2.0
   \date       26.03.15
   \brief      Split a PDB file into separate chains
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 1997-2015
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
-  V1.0    10.07.97  Original   By: ACRM
-  V1.1    16.01.14  Fixed to handle HETATMs properly
-  V1.2    22.07.14  Renamed deprecated functions with bl prefix.
                     Added doxygen annotation. By: CTP
-  V1.3    06.11.14  Renamed from splitchains By: ACRM
-  V1.4    12.03.15  Checks blank chain as string
-  V2.0    26.03.15  Major rewrite to handled whole PDB

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
BOOL WriteEachPDBChain(char *InFile, WHOLEPDB *wpdb, BOOL current);
void Usage(void);
BOOL BuildFileName(char *OutFile, int maxFileName, char *InFile, 
                   char *chain, BOOL current);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for extracting selected chains from a PDB file

-  07.02.97 Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   char InFile[MAXBUFF];
   FILE *in  = stdin;
   WHOLEPDB *wpdb;
   BOOL current;

   
   
   if(ParseCmdLine(argc, argv, InFile, &current))
   {
      if(blOpenStdFiles(InFile, NULL, &in, NULL))
      {
         if((wpdb=blReadWholePDB(in))==NULL)
         {
            if(!gQuiet)
               fprintf(stderr,"No atoms read from input PDB file\n");
            return(1);
         }
         else
         {
            if(!WriteEachPDBChain(InFile,wpdb,current))
            {
               if(!gQuiet)
                  fprintf(stderr,"pdbsplitchains: Failed to write all \
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
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \return                     Success

   Parse the command line

-  10.07.97 Original    By: ACRM
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
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  10.07.97 Original
-  16.01.14 V1.1
-  22.07.14 V1.2 By: CTP
-  12.03.15 V1.4
-  26.03.15 V2.0
*/
void Usage(void)
{
   fprintf(stderr,"pdbsplitchains V2.0 (c) 1997-2015 \
Dr. Andrew C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: pdbsplitchains [-c][-q] [in.pdb]\n");
   fprintf(stderr,"       -c  Output to current directory\n");
   fprintf(stderr,"       -q  Quiet - no error messages\n");

   fprintf(stderr,"\npdbsplitchains takes a PDB file (of the specified name \
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
*//**

   \param[in]      maxOutFileName    Max buffer size for output filename
   \param[in]      *InFile           Input filename
   \param[in]      *chain            Chain name to put into filename
   \param[in]      current           Strip path info so we write to 
   \param[in]      directory
   \param[out]     *OutFile          Output filename
   \return                     success?

   Constructs a filename from a current filename plus a chain name

-  10.07.97 Original
-  16.01.14 Uses blStrncat() and chain is handled as a string
-  12.03.15 Checks blank chain as string
*/
BOOL BuildFileName(char *OutFile, int maxOutFileName, char *InFile, 
                   char *chain, BOOL current)
{
   char *chp,
        *chpt,
        *WorkFile,
        *stem,
        *path;

   if(CHAINMATCH(chain, " "))
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
         blStrncat(OutFile, path, maxOutFileName);
         blStrncat(OutFile, "/",  maxOutFileName);
      }
      blStrncat(OutFile, stem,  maxOutFileName); /* The filestem        */
      blStrncat(OutFile, chain, maxOutFileName); /* The chain name      */
      blStrncat(OutFile,".pdb", maxOutFileName); /* The extension       */
      
      /* Free allocated memory                                          */
      free(WorkFile);
   }
   else
   {
      strcpy(OutFile, chain);
      blStrncat(OutFile,".pdb", maxOutFileName);
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL WriteEachPDBChain(char *InFile, WHOLEPDB *wpdb, BOOL currentDir)
   ---------------------------------------------------------------------
*//**

   \param[in]      *InFile     Input filename
   \param[in]      *wpdb       Whole PDB linked list
   \param[in]      currentDir  Strip path and write to current directory

   Writes each chain to a separate file

   Builds a list of the chain IDs using a STRINGLIST
   Uses blGetPDBChainAsCopy() to copy each chain in turn, attaching it to
   the whole PDB structure and writing it to a file.
   Then frees the copy and moves onto the next chain.

-  16.01.14 Rewritten to deal with HETATMs properly   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  25.03.15 Complete rewrite to support whole PDB  By: ACRM
*/
BOOL WriteEachPDBChain(char *InFile, WHOLEPDB *wpdb, BOOL currentDir)
{
   PDB        *pdb,
              *chain;
   FILE       *fp;
   char       **chainLabels = NULL;
   int        nChains, chainNum;

   pdb  = wpdb->pdb;

   /* Build a list of chain labels that are used                        */
   if((chainLabels = blGetPDBChainLabels(pdb, &nChains))==NULL)
      return(FALSE);

   /* Step through the list of chain labels                             */
   for(chainNum=0; chainNum<nChains; chainNum++)
   {
      if(!gQuiet)
         fprintf(stderr,"Writing chain '%s'\n", chainLabels[chainNum]);
      
      /* Make a copy of this chain                                      */
      if((chain=blGetPDBChainAsCopy(pdb, chainLabels[chainNum]))!=NULL)
      {
         char OutFile[MAXBUFF];

         /* Link this into the whole pdb structure, build a filename
            and write the chain
         */
         wpdb->pdb = chain;
         if(BuildFileName(OutFile, MAXBUFF, InFile, 
                          chainLabels[chainNum], currentDir))
         {
            if((fp=fopen(OutFile, "w"))!=NULL)
            {
               blWriteWholePDB(fp, wpdb);
               fclose(fp);
            }
            else
            {
               if(!gQuiet)
                  fprintf(stderr,"pdbsplitchains: Could not write output \
file: %s\n", OutFile);
               return(FALSE);
            }
         }
         else
         {
            if(!gQuiet)
               fprintf(stderr,"pdbsplitchains: No memory to build output \
filename\n");
            return(FALSE);
         }

         /* Free the storage for the copy of this chain and restore the
            whole PDB structure
         */
         FREELIST(chain, PDB);
         wpdb->pdb = pdb;
      }
   }
   return(TRUE);
}


