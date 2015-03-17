/************************************************************************/
/**

   \file       pdbgetchain.c
   
   \version    V2.1
   \date       13.03.15
   \brief      Extract chains from a PDB file
   
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
-  V1.0  07.02.97 Original   By: ACRM
-  V1.1  24.07.97 Modified such that if the chain is specified as `0'
                  gets that chain if it exists; otherwise gets everything
-  V1.2  06.01.99 Added -n parameter which gets chains by sequential
                  number
-  V1.3  06.04.09 Added -l parameter which retains lower case chain names
-  V1.4  21.05.09 Changed to use ReadWholePDB() with -k
-  V1.5  29.06.09 Added -a parameter for ATOMs only (discards HETATMs)
-  V1.6  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.7  06.11.14 Renamed from getchain  By: ACRM
-  V1.8  13.02.15 Removed -k option - headers are always kept
-  V2.0  04.03.15 Major rewrite to modify the PDB linked list in place
                  and write it as a whole rather than writing records
                  as we go. Also added support for multi-character
                  chain labels. If called as pdbgetchain, the chains must
                  be comma separated. If called as getchain, only the old
                  single character chain labels (not comma-separated) are
                  supported.
-  V2.1  13.03.15 Modified to use bioplib routines for list parsing

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "bioplib/array.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF  160
#define MAXCHAIN 160
#define MAXCHAINLABEL 8

/* Converts an integer to a character representation of that integer.
   The input integer is expected to be between 1 and 10. Behaviour
   is undefined outside this range.
-  06.01.99 Original  By: ACRM
*/
#define ITOCHAR(x) ((x)==10?'0':'0'+(char)(x))

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
char **ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                    BOOL *numeric, BOOL *atomsOnly);
PDB *FindEndOfChain(PDB *chain);
void SelectPDBChains(WHOLEPDB *wpdb, char **chains, BOOL numeric);
BOOL ValidChain(PDB *pdb, char **chains, BOOL numeric);



/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for extracting selected chains from a PDB file

-  07.02.97 Original   By: ACRM
-  06.04.09 Added lowercase option
-  22.05.09 Added keepHeader
-  29.06.09 Added atomsOnly
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  13.02.15 Now always keeps header  By: ACRM
- Removed lowercase
*/
int main(int argc, char **argv)
{
   char InFile[MAXBUFF],
        OutFile[MAXBUFF],
        **chains = NULL;
   FILE *in  = stdin,
        *out = stdout;
   BOOL numeric = FALSE,
        atomsOnly = FALSE;
   WHOLEPDB *wpdb = NULL;
   
   if((chains = ParseCmdLine(argc, argv, InFile, OutFile, &numeric,
                             &atomsOnly))!=NULL)
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         /* 29.06.09 Added atomsOnly option                             */
         if(atomsOnly)
         {
            wpdb=blReadWholePDBAtoms(in);
         }
         else
         {
            wpdb=blReadWholePDB(in);
         }
         
         if((wpdb == NULL)||
            (wpdb->pdb == NULL))
         {
            fprintf(stderr,"No atoms read from input PDB file\n");
            return(1);
         }
         else
         {
            SelectPDBChains(wpdb, chains, numeric);
            blWriteWholePDB(out, wpdb);
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
void SelectPDBChains(WHOLEPDB *wpdb, char **chains, BOOL numeric)
{
   PDB  *chainStart    = NULL,
        *endOfChain    = NULL,
        *nextChain     = NULL,
        *keptChainsEnd = NULL;
   
   for(chainStart=wpdb->pdb; chainStart!=NULL; chainStart=nextChain)
   {
      endOfChain       = FindEndOfChain(chainStart);
      nextChain        = endOfChain->next;
      endOfChain->next = NULL;

      if(ValidChain(chainStart, chains, numeric))
      {
         if(keptChainsEnd!=NULL)
            keptChainsEnd->next = chainStart;
         keptChainsEnd = endOfChain;
      }
      else
      {
         if(chainStart == wpdb->pdb)
            wpdb->pdb = nextChain;

         FREELIST(chainStart, PDB);
      }
   }
}

/************************************************************************/
BOOL ValidChain(PDB *pdb, char **chains, BOOL numeric)
{
   int i;
   static int sChainCount = 0;

   if(numeric)
   {
      sChainCount++;
      for(i=0; ((chains[i] != NULL) && (chains[i][0] != '\0')); i++)
      {
         int numericChain;
         
         if(sscanf(chains[i], "%d", &numericChain))
         {
            if(sChainCount == numericChain)
               return(TRUE);
         }
      }
   }
   else
   {
      for(i=0; ((chains[i] != NULL) && (chains[i][0] != '\0')); i++)
      {
         if(CHAINMATCH(pdb->chain, chains[i]))
            return(TRUE);
      }
   }
   
   return(FALSE);
}

/************************************************************************/
PDB *FindEndOfChain(PDB *chain)
{
   PDB *p;
   if(chain==NULL)
      return(NULL);
   
   for(p=chain; 
       ((p->next !=NULL) && CHAINMATCH(p->chain, p->next->chain));
       NEXT(p))
   {
      continue;
   }
   return(p);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  07.02.97 Original   By: ACRM
-  24.07.97 V1.1
-  06.01.98 V1.2
-  06.04.09 V1.3
-  22.05.09 V1.4
-  29.06.09 V1.5
-  22.07.14 V1.6 By: CTP
-  06.11.14 V1.7 By: ACRM
-  13.02.15 V1.8
-  04.03.15 V2.0
-  13.03.15 V2.1
*/
void Usage(void)
{
   fprintf(stderr,"\npdbgetchain V2.1 (c) 1997-2015 Dr. Andrew C.R. \
Martin, UCL\n");

   fprintf(stderr,"\nUsage: pdbgetchain [-n] [-a] \
chain[,chain[...]] [in.pdb [out.pdb]]\n");
   fprintf(stderr,"       -n Specify chains numerically: 1 is the first \
chain, 2 the\n");
   fprintf(stderr,"          second, etc.\n");
   fprintf(stderr,"       -a ATOMs only (discard HETATMs)\n");

   fprintf(stderr,"\npdbgetchain reads a PDB file and write out only \
those chains specified\n");
   fprintf(stderr,"on the command line. If input and output filenames \
are not given\n");
   fprintf(stderr,"I/O is through standard input/output.\n");
   fprintf(stderr,"\nThe -k (keep headers) and -l (take lowercase chain \
names) options in\n");
   fprintf(stderr,"previous versions are now deprecated.\n");
   fprintf(stderr,"\nHeaders may contain references to chains that are \
no longer present.\n\n");
   fprintf(stderr,"If the program is called as getchain rather than \
pdbgetchain, the old\n");
   fprintf(stderr,"behaviour of only accepting one-character chain names \
and taking them\n");
   fprintf(stderr,"as a non-comma separated set is used. e.g. chains L \
and H, would be\n");
   fprintf(stderr,"specified as LH rather than L,H\n\n");
} 


/************************************************************************/
/*>char **ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                       BOOL *numeric, BOOL *atomsOnly)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \param[out]     *outfile    Output filename (or blank string)
   \param[out]     *numeric    Chains are specified numerically
   \param[out]     *atomsOnly  Discard HETATMs
   \return                     2D array of chain labels

   Parse the command line

-  07.02.97 Original    By: ACRM
-  06.01.99 Added -n and numeric parameter
-  06.04.09 Added -l lowercase option
-  22.05.09 Added -k option
-  29.06.09 Added -a option
-  04.03.15 Removed -l option
-  13.03.15 Switched to new bioplib list parsing routines
*/
char **ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                    BOOL *numeric, BOOL *atomsOnly)
{
   char **chains = NULL;
   BOOL oldStyle = FALSE;
   
   oldStyle = blCheckProgName(argv[0], "getchain");

   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'n':
            *numeric = TRUE;
            break;
         case 'l':
            fprintf(stderr,"The -l option is now deprecated\n");
            break;
         case 'k':
            fprintf(stderr,"The -k option is now deprecated\n");
            break;
         case 'a':
            *atomsOnly = TRUE;
            break;
         case 'h':
            return(NULL);
            break;
         default:
            return(NULL);
            break;
         }
      }
      else
      {
         /* Check that there are 1-3 arguments left                     */
         if(argc > 3 || argc < 1)
            return(NULL);
         
         /* Parse the chains out of the first argument                  */
         if(oldStyle)
         {
            if((chains = blSplitStringOnChars(argv[0]))
               ==NULL)
            {
               fprintf(stderr,"No memory for storing chain labels: %s\n",
                       argv[0]);
               exit(1);
            }
         }
         else
         {
            if((chains = blSplitStringOnCommas(argv[0], MAXCHAINLABEL))
               ==NULL)
            {
               fprintf(stderr,"No memory for storing chain labels: %s\n",
                       argv[0]);
               exit(1);
            }
         }
         
         /* If there's another, copy it to infile                       */
         argc--;
         argv++;
         if(argc)
         {
            strcpy(infile, argv[0]);
         
            /* If there's another, copy it to outfile                   */
            argc--;
            argv++;
            if(argc)
               strcpy(outfile, argv[0]);
         }

         return(chains);
      }
      argc--;
      argv++;
   }
   
   return(chains);
   
}

