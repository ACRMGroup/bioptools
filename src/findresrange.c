/************************************************************************/
/**

   \file       findresrange.c
   
   \version    V1.2
   \date       22.07.14
   \brief      Find a residue range given a key residue and a number
               of residues on either side
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2010-2014
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
-  V1.0   19.05.10  Original   By: ACRM
-  V1.1   02.06.10  Fixed junk printing when out of range
-  V1.2   22.07.14  Renamed deprecated functions with bl prefix.
                    Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>

#include "bioplib/pdb.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 320
#define RESBUFF 40

#define NO_STARTRES    1
#define NO_ENDRES      2
#define NO_KEYRESPARSE 3
#define NO_MEMORY      4
#define NO_KEYRES      5

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *keyres, int *width);
int GetResidueRange(PDB *pdb, char *keyres, int width, char *startres, 
                    char *endres);
void Usage(void);
int main(int argc, char **argv);


/************************************************************************/
int main(int argc, char **argv)
{
   char keyres[RESBUFF],
        startres[RESBUFF],
        endres[RESBUFF];
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   int  width, error, natoms;
   FILE *in  = stdin,
        *out = stdout;
   PDB  *pdb;
   
   if(ParseCmdLine(argc, argv, infile, outfile, keyres, &width))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDB(in, &natoms)) == NULL)
         {
            fprintf(stderr, "Unable to read PDB file\n");
            return(1);
         }
         
         if((error=GetResidueRange(pdb, keyres, width, startres, endres))
            == 0)
         {
            fprintf(out, "%s %s\n", startres, endres);
            return(0);
         }
         else
         {
            switch(error)
            {
            case NO_MEMORY:
               fprintf(out, "Error: No memory for PDB structure\n");
               break;
            case NO_KEYRESPARSE:
               fprintf(out, "Error: Illegal key residue specification: %s\n",
                       keyres);
               break;
            case NO_KEYRES:
               fprintf(out, "Error: Key residue %s not found\n", keyres);
               break;
            case NO_STARTRES:
               fprintf(out, "Error: No residue %d before key residue %s\n", 
                       width, keyres);
               break;
            case NO_ENDRES:
               fprintf(out, "Error: No residue %d after key residue %s\n", 
                       width, keyres);
               break;
            }
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
/*>int GetResidueRange(PDB *pdb, char *keyres, int width, char *startres, 
                       char *endres)
   -----------------------------------------------------------------------
*//**

-  02.06.10 Returns NO_KEYRES properly
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int GetResidueRange(PDB *pdb, char *keyres, int width, char *startres, 
                    char *endres)
{
   PDBSTRUCT  *pdbs;
   PDBCHAIN   *pdbc;
   PDBRESIDUE *pdbr, *endr, *startr;
   char       chain[8], insert[8];
   int        resnum, i;
   BOOL       found = FALSE;
   
   
   if((pdbs = blAllocPDBStructure(pdb))==NULL)
      return(NO_MEMORY);

   if(blParseResSpec(keyres, chain, &resnum, insert))
   {
      for(pdbc = pdbs->chains; pdbc!=NULL; NEXT(pdbc))
      {
         if(pdbc->chain[0] == chain[0])
         {
            for(pdbr = pdbc->residues; pdbr!=NULL; NEXT(pdbr))
            {
               if((pdbr->resnum == resnum) && 
                  (pdbr->insert[0] == insert[0]))
               {
                  startr = endr = pdbr;
                  found = TRUE;
                  for(i=0; i<width; i++)
                  {
                     if(endr!=NULL)   endr   = endr->next;
                     if(startr!=NULL) startr = startr->prev;
                  }
                  if(startr == NULL) return(NO_STARTRES);
                  if(endr   == NULL) return(NO_ENDRES);

                  strncpy(startres, startr->resid, RESBUFF);
                  strncpy(endres,   endr->resid,   RESBUFF);

                  return(0);
               }
            }
            break;
         }
      }
      return(NO_KEYRES);
   }
   else
   {
      return(NO_KEYRESPARSE);
   }

   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     char *keyres, int *width)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *patchfile   Output file (or blank string)
   \return                     Success?

   Parse the command line
   
-  19.05.10 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *keyres, int *width)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';

   if(argc < 3)
      return(FALSE);
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2--4 arguments left                    */
         if(argc < 2 || argc > 4)
            return(FALSE);
         
         /* Copy the first to keyres                                    */
         strcpy(keyres, argv[0]);
         argc--;
         argv++;

         /* Copy the next to width                                      */
         sscanf(argv[0], "%d", width);
         argc--;
         argv++;
         
         if(argc)
         {
            /* Copy the next to infile                                  */
            strcpy(infile, argv[0]);
            argc--;
            argv++;
            
            /* If there's another, copy it to outfile                   */
            if(argc)
               strcpy(outfile, argv[0]);
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

-  22.07.14 V1.2 By: CTP
*/
void Usage(void)
{
   printf("\nfindresrange V1.2 (c) 2010-2014 UCL, Andrew C.R. Martin\n");
   printf("\nUsage: findresrange keyres width [input.pdb \
[output.txt]]\n");
   printf("\nTakes a PDB file as input and given:\n");
   printf("1. a key residue (keyres) specified in the format \
[chain]resnum[insert]\n");
   printf("(where chain and insert are optional and chain may be \
followed by a '.'\n");
   printf("if it is numeric)\n");
   printf("2. a number of residues (width)\n");
   printf("will return the residue identifiers for the residues width \
before and\n");
   printf("width after the key residue.\n\n");
}

