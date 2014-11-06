/************************************************************************/
/**

   \file       ChainPDB.c
   
   \version    V1.7
   \date       22.07.14
   \brief      Insert chain labels into a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin 1994-2014
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

   N.B. This version only supports protein entries (not DNA).


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

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define MAXCHAIN 64

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
                  char *chains, BOOL *BumpChainOnHet);
void Usage(void);
void DoChain(PDB *pdb, char *chains, BOOL BumpChainOnHet);


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
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        chains[MAXCHAIN];
   PDB  *pdb;
   int  natoms;
   BOOL BumpChainOnHet = FALSE;  /* Flag to indicate ATOMs after HETATMs
                                    should be start of a new residue    */

        
   if(ParseCmdLine(argc, argv, infile, outfile, chains, &BumpChainOnHet))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natoms))==NULL)
         {
            fprintf(stderr,"No atoms read from input file\n");
         }
         else
         {
            DoChain(pdb,chains,BumpChainOnHet);
            blWritePDB(out, pdb);
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
                     char *chains, BOOL *BumpChainOnHet)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc             Argument count
   \param[in]      **argv           Argument array
   \param[out]     *infile          Input filename (or blank string)
   \param[out]     *outfile         Output filename (or blank string)
   \param[out]     *chains          Chain labels
   \param[out]     *BumpChainOnHet  Bump chain label on ATOMs after HETs
   \return                          Success

   Parse the command line

-  12.07.94 Original    By: ACRM
-  16.10.95 Sets BumpChainOnHet
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *chains, BOOL *BumpChainOnHet)
{
   argc--;
   argv++;
   
   *BumpChainOnHet = FALSE;
   infile[0] = outfile[0] = '\0';
   chains[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'c':
            argc--;
            argv++;
            strncpy(chains,argv[0],MAXCHAIN);
            chains[MAXCHAIN-1] = '\0';
            UPPER(chains);
            break;
         case 'b':
            *BumpChainOnHet = TRUE;
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
*/
void Usage(void)
{
   fprintf(stderr,"\nChainPDB V1.7 (c) 1994-2014 Dr. Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Freely distributable if no profit is made\n");
   fprintf(stderr,"Splits a PDB file into chains using distance \
criteria\n\n");
   
   fprintf(stderr,"Usage: chainpdb [-c <chains>] [<infile> \
[<outfile>]]\n");
   fprintf(stderr,"       -c Specify chain names to use\n");
   fprintf(stderr,"       -b If ATOM records follow HETATM records they \
start a new chain\n\n");
   fprintf(stderr,"If files are not specified, stdin and stdout are \
used.\n");
   fprintf(stderr,"If a chain is to be skipped with -c, use a - \
instead of the label or\nnumber.\n\n");
}


/************************************************************************/
/*>void DoChain(PDB *pdb, char *chains, BOOL BumpChainOnHet)
   ---------------------------------------------------------
*//**

   \param[in,out]  *pdb            PDB linked list
   \param[in]      *chains         Chain labels (or blank string)

   Do the actual chain naming.

-  12.07.94 Original    By: ACRM
-  25.07.94 Only increments ch if *ch != \0
-  04.01.95 Added check on HETATM records
-  27.01.95 ChainNum count now mod 26 so labels will cycle A-Z
-  16.10.95 Handles BumpChainOnHet
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
void DoChain(PDB *pdb, char *chains, BOOL BumpChainOnHet)
{
   PDB  *p,
        *start,
        *end,
        *LastStart = NULL,
        *N        = NULL,
        *C        = NULL,
        *CPrev    = NULL,
        *CAPrev   = NULL,
        *CA       = NULL;
   int  ChainNum  = 0;
   char chain,
        *ch;
   BOOL NewChain;
   
   ch = chains;

   chain = (*ch) ? *ch : 'A';
   
   for(start=pdb; start!=NULL; start=end)
   {
      NewChain = FALSE;
      end = blFindEndPDB(start);

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
            sprintf(buffer,"residues %c%d%c and %c%d%c",
                    LastStart->chain[0],LastStart->resnum,
                    LastStart->insert[0],
                    start->chain[0],start->resnum,start->insert[0]);
         else if(CPrev == NULL || CAPrev == NULL)
            sprintf(buffer,"residue %c%d%c",
                    LastStart->chain[0],LastStart->resnum,
                    LastStart->insert[0]);
         else
            sprintf(buffer,"residue %c%d%c",
                    start->chain[0],start->resnum,start->insert[0]);

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
         if(*ch) ch++;
         
         if(*ch)
         {
            chain = *ch;
         }
         else
         {
            chain = (char)(65 + (ChainNum%26));
         }
      }

      /* Copy the name into this residue                                */
      for(p=start; p!=end; NEXT(p))
      {
         p->chain[0] = chain;
         p->chain[1] = '\0';
      }
      
      /* Set pointers for next residue                                  */
      CAPrev    = CA;
      CPrev     = C;
      LastStart = start;
   }
}
