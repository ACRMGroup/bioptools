/************************************************************************/
/**

   \file       pdbgetchain.c
   
   \version    V1.8
   \date       13.02.15
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

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF  160
#define MAXCHAIN 160

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
void WritePDBChains(FILE *out, WHOLEPDB *wpdb, char *chains, 
                    BOOL numeric);
int  WritePDBChainsByNumber(FILE *out, PDB *pdb, char *chains);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *chains, BOOL *numeric, BOOL *lowercase,
                  BOOL *atomsOnly);


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
*/
int main(int argc, char **argv)
{
   char InFile[MAXBUFF],
        OutFile[MAXBUFF],
        chains[MAXCHAIN];
   FILE *in  = stdin,
        *out = stdout;
   BOOL numeric = FALSE,
        lowercase = FALSE,
        atomsOnly = FALSE;
   WHOLEPDB *wpdb = NULL;
   
   if(ParseCmdLine(argc, argv, InFile, OutFile, chains, &numeric,
                   &lowercase, &atomsOnly))
   {
      /* 06.04.09 Added lowercase option                                */
      if(!lowercase)
      {
         UPPER(chains);
      }
      
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
            WritePDBChains(out, wpdb, chains, numeric);
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
/*>void WritePDBChains(FILE *out, WHOLEPDB *wpdb, char *chains, 
                       BOOL numeric)
   ------------------------------------------------------------
*//**

   \param[in]      *out      Output file pointer
   \param[in]      *pdb      WHOLEPDB structure containing PDB linked 
   \param[in]      *chains   Chain names to be written to output file
   \param[in]      numeric   Chains are specified numerically

   Like blWritePDB(), but takes a list of chain names to be written. Other
   chains are not written to the output.

   If chain `0' is specified and no chain 0 exists, then the empty chain
   name is written

-  07.02.97 Original   By: ACRM
-  24.07.97 Modified to handle chain 0 as being the empty chain if
            chain 0 doesn't exist
-  06.01.99 Added numeric parameter
-  22.05.09 Now takes WHOLEPDB rather than PDB with keepHeader
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  13.02.15 Now always keeps header  By: ACRM
-  23.02.15 Modified to use blWriteTerCard() and new version of 
            blWriteWholePDBTrailer()
*/
void WritePDBChains(FILE *out, WHOLEPDB *wpdb, char *chains, BOOL numeric)
{
   PDB *p;
   char PrevChain;
   BOOL FoundZero = FALSE,
        Written   = FALSE;
   PDB  *pdb = wpdb->pdb;
   int  numTer = 0;
   
   blWriteWholePDBHeader(out, wpdb);

   if(numeric)
   {
      numTer = WritePDBChainsByNumber(out, pdb, chains);
   }
   else
   {
      PDB *prev = NULL;
      PrevChain = '\0';
      
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(strchr(chains,p->chain[0]))
         {
            if((PrevChain != p->chain[0]) && (PrevChain != '\0'))
            {
               /* Chain change, insert TER card                         */
               blWriteTerCard(out, prev);
               numTer++;
            }
            blWritePDBRecord(out,p);
            PrevChain = p->chain[0];
            Written = TRUE;
         }
         if(p->chain[0] == '0')
            FoundZero = TRUE;

         prev=p;
      }
      if(Written)
      {
         Written = FALSE;
         blWriteTerCard(out, prev);
         numTer++;
      }
      
      /* If we were looking for a chain zero and there wasn't a chain 
         with that name, then we write the chain with a blank chain name
      */
      if(strchr(chains,'0') && !FoundZero)
      {
         for(p=pdb; p!=NULL; NEXT(p))
         {
            if(p->chain[0] == ' ')
            {
               blWritePDBRecord(out,p);
               Written = TRUE;
            }
         }
         if(Written)
         {
            Written = FALSE;
            blWriteTerCard(out, prev);
            numTer++;
         }
      }
   }
   
   blWriteWholePDBTrailer(out, wpdb, numTer);
}


/************************************************************************/
/*>BOOL WritePDBChainsByNumber(FILE *out, PDB *pdb, char *chains)
   --------------------------------------------------------------
*//**

   \param[in]      *out     Output file pointer
   \param[in]      *pdb     PDB linked list
   \param[in]      *chains  Chain numbers to be written to output file

   Like blWritePDB(), but takes a list of chain numbers to be written. 1
   is the first chain, 2 the second up to 0 being the 10th. Other
   chains are not written to the output.

-  06.01.99 Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  23.02.15 Modified to use blWriteTerCard() and new version of 
            blWriteWholePDBTrailer()
            Now returns number of TER cards written
*/
int WritePDBChainsByNumber(FILE *out, PDB *pdb, char *chains)
{
   PDB  *p,
        *prev=NULL;
   char PrevChain,
        ThisChainChar;
   int  ThisChainNum = 1,
        numTer       = 0;
   BOOL Written      = FALSE,
        Warned       = FALSE;
   

   PrevChain = pdb->chain[0];
   ThisChainChar = ITOCHAR(ThisChainNum);
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(PrevChain != p->chain[0])
      {
         /* Chain change, insert TER card                               */
         if(Written)
         {
            blWriteTerCard(out, prev);
            numTer++;
            Written = FALSE;
         }
         PrevChain = p->chain[0];
         ThisChainNum++;
         ThisChainChar = ITOCHAR(ThisChainNum);
         if(!Warned && (ThisChainChar > '9'))
         {
            fprintf(stderr,"Warning: More than 10 chains in PDB file\n");
            Warned = TRUE;
         }
      }
      if(strchr(chains,ThisChainChar))
      {
         blWritePDBRecord(out,p);
         Written = TRUE;
      }
      prev = p;
   }

   if(Written)
   {
      blWriteTerCard(out, prev);
      numTer++;
   }

   return(numTer);
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
*/
void Usage(void)
{
   fprintf(stderr,"\npdbgetchain V1.8 (c) 1997-2015 Dr. Andrew C.R. \
Martin, UCL\n");

   fprintf(stderr,"\nUsage: pdbgetchain [-n] [-l] [-a] chains \
[in.pdb [out.pdb]]\n");
   fprintf(stderr,"       -n Specify chains numerically: 1 is the first \
chain, 2 the second,\n");
   fprintf(stderr,"          etc. Maximum chain number is therefore 0 \
(representing 10) such\n");
   fprintf(stderr,"          that 120 would represent chains 1, 2 and \
10\n");
   fprintf(stderr,"       -l Retain lower case chain names\n");
   fprintf(stderr,"       -a ATOMs only (discard HETATMs)\n");

   fprintf(stderr,"\npdbgetchain reads a PDB file and write out only \
those chains specified\n");
   fprintf(stderr,"on the command line. If input and output filenames \
are not given\n");
   fprintf(stderr,"I/O is through standard input/output\n");
   fprintf(stderr,"\nIf chain 0 is specified and there is no chain of \
that name then any chain\n");
   fprintf(stderr,"with a blank chain name will be written.\n\n");
   fprintf(stderr,"The headers are kept (the -k option in previous \
versions is now deprecated)\n");
   fprintf(stderr,"but may contain references to chains that are no \
longer present.\n\n");
} 


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     char *chains, BOOL *numeric, BOOL *lowercase,
                     BOOL *atomsOnly)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \param[out]     *outfile    Output filename (or blank string)
   \param[out]     *chains     Chain names to be written
   \param[out]     *numeric    Chains are specified numerically
   \param[out]     *lowercase  Chain names may be in lower case
   \param[out]     *atomsOnly  Discard HETATMs
   \return                     Success

   Parse the command line

-  07.02.97 Original    By: ACRM
-  06.01.99 Added -n and numeric parameter
-  06.04.09 Added -l lowercase option
-  22.05.09 Added -k option
-  29.06.09 Added -a option
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *chains, BOOL *numeric, BOOL *lowercase,
                  BOOL *atomsOnly)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = chains[0] = '\0';
   
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
            *lowercase = TRUE;
            break;
         case 'k':
            fprintf(stderr,"The -k option is now deprecated\n");
            break;
         case 'a':
            *atomsOnly = TRUE;
            break;
         case 'h':
            return(FALSE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 1-3 arguments left                     */
         if(argc > 3)
            return(FALSE);
         
         /* Copy the first to chains                                    */
         strcpy(chains, argv[0]);
         
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

         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   if(!chains[0])
      return(FALSE);

   return(TRUE);
}

