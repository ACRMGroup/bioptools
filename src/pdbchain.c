/************************************************************************/
/**

   \file       pdbchain.c
   
   \version    V2.3
   \date       09.08.19
   \brief      Insert chain labels into a PDB file
   
   \copyright  (c) UCL, Prof. Andrew C. R. Martin 1994-2019
   \author     Prof. Andrew C. R. Martin
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

   N.B. This version only supports protein entries (not DNA).
        This version only handles 1-letter chain labels

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  12.07.94 Original
-  V1.1  25.07.94 Fixed bug for C compilers which don't initialise
                  arrays to null characters
-  V1.2  24.08.94 Changed to call OpenStdFiles()
-  V1.3  05.01.95 Added check on HETATM when printing missing atoms
-  V1.4  27.01.95 Chain labels restricted to A--
-  V1.5  16.10.95 Was failing to renumber chains immediately following
                  HETATM records.
-  V1.6  30.05.02 Changed PDB field from 'junk' to 'record_type'
-  V1.7  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.8  06.11.14 Renamed from chainpdb By: ACRM
-  V1.9  13.02.15 Added whole PDB suport
-  V1.10 05.03.15 Replaced blFindEndPDB() with blFindNextResidue()
-  V2.0  10.03.15 Now generates more than 26 chain labels and, if 
                  specfied with -c expects a comma-separated list of
                  chain labels instead of a simple string. i.e chains
                  L and H are now specified as L,H instead of LH
-  V2.1  13.03.15 Now supports old chain specification method if
                  called as chainpdb
-  V2.2  28.01.18 Increased MAXCHAINLABEL from 8 to 16
-  V2.3  09.08.19 Added -v flag

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/array.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define MAXCHAIN 160
#define MAXCHAINLABEL 16

#define CNDISTSQ  3.5
#define CADISTSQ 16.0

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char ***chains, BOOL *BumpChainOnHet, BOOL *verbose);
void Usage(void);
void DoChain(PDB *pdb, char **chains, BOOL BumpChainOnHet, BOOL verbose);
char *GetChainLabel(int ChainNum);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for PDB chain making

-  12.07.94 Original    By: ACRM
-  24.08.94 Changed to call OpenStdFiles()
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE     *in  = stdin,
            *out = stdout;
   char     infile[MAXBUFF],
            outfile[MAXBUFF],
            **chains;
   PDB      *pdb;
   WHOLEPDB *wpdb;
   BOOL     BumpChainOnHet = FALSE;  /* Flag to indicate ATOMs after 
                                        HETATMs should be start of a new
                                        residue                         */
   BOOL     verbose = FALSE;         /* Flag to print CA-CA distances   */


   if(ParseCmdLine(argc, argv, infile, outfile, &chains, &BumpChainOnHet,
                   &verbose))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb=blReadWholePDB(in))==NULL)
         {
            fprintf(stderr,"No atoms read from input file\n");
         }
         else
         {
            pdb = wpdb->pdb;
            DoChain(pdb, chains, BumpChainOnHet, verbose);
            blWriteWholePDB(out, wpdb);
         }
      }
      else
      {
         Usage();
         return(1);
      }
   }
   else
   {
      Usage();
      return(1);
   }
   
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     char ***chains, BOOL *BumpChainOnHet, BOOL *verbose)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc             Argument count
   \param[in]      **argv           Argument array
   \param[out]     *infile          Input filename (or blank string)
   \param[out]     *outfile         Output filename (or blank string)
   \param[out]     ***chains        Array of chain labels (if specified)
                                    NULL if no chains allocated
   \param[out]     *BumpChainOnHet  Bump chain label on ATOMs after HETs
   \param[out]     *verbose         Print CA-CA distance of chain breaks
   
   \return                          Success

   Parse the command line

-  12.07.94 Original    By: ACRM
-  16.10.95 Sets BumpChainOnHet
-  10.03.15 Changed chains to strings. Now returns that array
-  13.03.15 Gone back to returning BOOL
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char ***chains, BOOL *BumpChainOnHet, BOOL *verbose)
{
   BOOL oldStyle;

   oldStyle = blCheckProgName(argv[0], "chainpdb");
   *chains = NULL;
  
   argc--;
   argv++;
   
   *BumpChainOnHet = FALSE;
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'c':
            argc--;
            argv++;

            if(oldStyle)
            {
               if((*chains = blSplitStringOnChars(argv[0]))
                  ==NULL)
               {
                  fprintf(stderr,"No memory for storing chain labels: %s\n",
                          argv[0]);
                  exit(1);
               }
            }
            else
            {
               if((*chains = blSplitStringOnCommas(argv[0], MAXCHAINLABEL))
                  ==NULL)
               {
                  fprintf(stderr,"No memory for storing chain labels: %s\n",
                          argv[0]);
                  exit(1);
               }
            }
            break;
         case 'b':
            *BumpChainOnHet = TRUE;
            break;
         case 'v':
            *verbose        = TRUE;
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
*//**

   Print a usage message

-  12.07.94 Original    By: ACRM
-  25.07.94 V1.1
-  24.08.94 V1.2
-  05.01.95 V1.3
-  27.01.95 V1.4
-  16.10.95 V1.5
-  30.05.02 V1.6
-  22.07.14 V1.7 By: CTP
-  06.11.14 V1.8 By: ACRM
-  05.03.15 V1.10
-  10.03.15 V2.0
-  13.03.15 V2.1
-  28.01.18 V2.2
-  09.08.19 V2.3
*/
void Usage(void)
{
   fprintf(stderr,"\npdbchain V2.3 (c) 1994-2019 Prof. Andrew C.R. \
Martin, UCL\n");

   fprintf(stderr,"\nUsage: pdbchain [-c chain[,chain[...]]][-b][-v] \
[in.pdb [out.pdb]]\n");
   fprintf(stderr,"       -c Specify chain names to use\n");
   fprintf(stderr,"       -b If ATOM records follow HETATM records they \
start a new chain\n");
   fprintf(stderr,"       -v Print CA-CA distance of each chain break\n");

   fprintf(stderr,"\nSplits a PDB file into chains using distance \
criteria\n\n");
   fprintf(stderr,"If files are not specified, stdin and stdout are \
used.\n");
   fprintf(stderr,"If a chain is to be skipped with -c, use a - \
instead of the label or\nnumber.\n\n");
   fprintf(stderr,"Note that chain labels used in the headers will \
not be updated as this is\n");
   fprintf(stderr,"designed to be used with models and partial PDB \
files.\n");
   fprintf(stderr,"If called using the old name, 'chainpdb', the chain \
labels are\n");
   fprintf(stderr,"supplied without commas (e.g. chains L and H as LH \
instead of L,H)\n\n");
}


/************************************************************************/
/*>void DoChain(PDB *pdb, char **chains, BOOL BumpChainOnHet,
                BOOL verbose)
   ----------------------------------------------------------
*//**

   \param[in,out]  *pdb            PDB linked list
   \param[in]      *chains         Chain labels (or blank string)
   \param[in]      BumpChainOnHet  Bump the chain label when a HETATM
                                   is found
   \param[in]      verbose         Print CA-CA distance of chain breaks

   Do the actual chain naming.

-  12.07.94 Original    By: ACRM
-  25.07.94 Only increments ch if *ch != \0
-  04.01.95 Added check on HETATM records
-  27.01.95 ChainNum count now mod 26 so labels will cycle A-Z
-  16.10.95 Handles BumpChainOnHet
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  05.03.15 Replaced blFindEndPDB() with blFindNextResidue()
-  10.03.15 Chains is now an array
-  09.08.19 Added verbose handling to print CA-CA distances
*/
void DoChain(PDB *pdb, char **chains, BOOL BumpChainOnHet, BOOL verbose)
{
   PDB  *p,
        *start,
        *end,
        *LastStart = NULL,
        *N         = NULL,
        *C         = NULL,
        *CPrev     = NULL,
        *CAPrev    = NULL,
        *CA        = NULL;
   int  ChainNum   = 0,
        ChainIndex = 0;
   char chain[MAXCHAINLABEL];
   BOOL NewChain;
   

   if((chains!=NULL) && chains[ChainIndex][0])
      strcpy(chain, chains[ChainIndex++]);
   else
      strcpy(chain, "A");
   
   for(start=pdb; start!=NULL; start=end)
   {
      NewChain = FALSE;
      end = blFindNextResidue(start);

      CA = N = C = NULL;
      
      for(p=start; p!=end; NEXT(p))
      {
         if(!strncmp(p->atnam,"CA  ",4)) CA  = p;
         if(!strncmp(p->atnam,"N   ",4)) N   = p;
         if(!strncmp(p->atnam,"C   ",4)) C   = p;
      }

      if(CPrev != NULL && N != NULL)
      {
         /* A C was defined in the last residue and an N in this one
            Calc C-N distance
         */
         if(DISTSQ(CPrev, N) > CNDISTSQ)
            NewChain = TRUE;
      }
      else if(CAPrev != NULL && CA != NULL)
      {
         /* No C-N connection, but a CAs found
            Calc CA-CA distance
         */
         if(DISTSQ(CAPrev, CA) > CADISTSQ)
            NewChain = TRUE;
      }
      else if(LastStart != NULL)
      {
         char buffer[80],
              atoms[80];
         
         /* Build string specifying faulty residues                     */
         if((CPrev == NULL || CAPrev == NULL) &&
            (N     == NULL || CA     == NULL))
            sprintf(buffer,"residues %s.%d%c and %s.%d%c",
                    LastStart->chain,LastStart->resnum,
                    LastStart->insert[0],
                    start->chain,start->resnum,start->insert[0]);
         else if(CPrev == NULL || CAPrev == NULL)
            sprintf(buffer,"residue %s.%d%c",
                    LastStart->chain,LastStart->resnum,
                    LastStart->insert[0]);
         else
            sprintf(buffer,"residue %s.%d%c",
                    start->chain,start->resnum,start->insert[0]);

         /* Build string specifying faulty atoms                        */
         atoms[0] = '\0';
         if(CAPrev == NULL || CA == NULL) strcat(atoms,"CA ");
         if(N      == NULL)               strcat(atoms,"N ");
         if(CPrev  == NULL)               strcat(atoms,"C ");

         /* Print warning message                                       */
         if((LastStart != NULL && 
             strncmp(LastStart->record_type, "HETATM", 6)) &&
            (start     != NULL && 
             strncmp(start->record_type,     "HETATM", 6)))
            fprintf(stderr, "Warning: Atoms missing in %s: %s\n",
                    buffer, atoms);

         if(BumpChainOnHet &&
            (LastStart != NULL && 
             !strncmp(LastStart->record_type, "HETATM", 6)) &&
            (start     != NULL && 
             !strncmp(start->record_type,     "ATOM  ",6)))
            NewChain = TRUE;
      }

      /* If we've changed chain, set the new chain name                 */
      if(NewChain)
      {
         ChainNum++;
         
         if((chains!=NULL) && chains[ChainIndex][0])
         {
            strcpy(chain,chains[ChainIndex++]);
         }
         else
         {
            strcpy(chain, GetChainLabel(ChainNum));
         }

         if(verbose)
         {
            if(CAPrev != NULL && CA != NULL)
            {
               REAL CaCaDist = DIST(CAPrev, CA);
               fprintf(stderr, "CA-CA distance to start of chain %s: \
%.3f\n", chain, CaCaDist);
            }
         }
      }

      /* Copy the name into this residue                                */
      for(p=start; p!=end; NEXT(p))
         strcpy(p->chain, chain);
      
      /* Set pointers for next residue                                  */
      CAPrev    = CA;
      CPrev     = C;
      LastStart = start;
   }
}


/************************************************************************/
/*>char *GetChainLabel(int ChainNum)
   ---------------------------------
*//**
   \param[in]  ChainNum    Chain number
   \return                 Chain label 

   Converts a chain number (>=0) into a chain label. Chain labels run
   from A-Z, a-z, 1-9, 0, and then 63 onwards as multi-character strings

-  10.03.15 Original   By: ACRM
*/
char *GetChainLabel(int ChainNum)
{
   static char chain[MAXCHAINLABEL];
   
   if(ChainNum < 26)
   {
      chain[0] = (char)(65 + ChainNum);
      chain[1] = '\0';
   }
   else if(ChainNum < 52)
   {
      chain[0] = (char)(97 + (ChainNum-26));
      chain[1] = '\0';
   }
   else if(ChainNum < 61)
   {
      sprintf(chain,"%d", ChainNum-51);
   }
   else if(ChainNum == 61)
   {
      strcpy(chain,"0");
   }
   else
   {
      sprintf(chain,"%d", ChainNum);
   }
   
   return(chain);
}

