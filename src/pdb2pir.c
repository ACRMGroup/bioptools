/************************************************************************/
/**

   \file       pdb2pir.c
   
   \version    V2.15
   \date       13.03.19
   \brief      Convert PDB to PIR sequence file
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL 1994-2019
   \author     Dr. Andrew C. R. Martin
   \par
               Biomolecular Structure and Modelling,
               University College London,

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

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
-  V1.0  10.05.94 Original   By: ACRM
-  V1.0a 11.02.97 Reformatted Usage message only
-  V2.0  22.08.97 Now also considers SEQRES records
-  V2.1  26.08.97 Reports the label (if specified) when a chain
                  from SEQRES is not found in ATOM records. Added -p
                  and fixed a few small bugs
-  V2.2  10.08.98 Rewrite of WritePIR() to fix bug with -p -c which
                  resulted in nothing printed for header of a chain
                  following a non-protein chain.
-  V2.3  02.10.00 Added -x option
-  V2.4  18.10.00 Added -f option
-  V2.5  08.02.01 Added -n option
-  V2.6  07.03.07 Now looks up MODRES records to find out original amino
                  acids. Also no longer does a rewind after reading
                  SEQRES so will work with stdin
-  V2.7  12.03.07 Revised gap penalty from 10 to 2
-  V2.8  22.05.09 Added -i option. Also chains in SEQRES but not in ATOM
                  now appear in lower case if -u not specified
-  V2.9  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V2.10 26.08.14 Use renamed macros blPDB2SeqXNoX() and blPDB2SeqX(). 
                  By: CTP
-  V2.11 07.11.14 Initialized a variable  By: ACRM
-  V2.12 25.11.14 Initialized a variable  By: ACRM
-  V2.13 10.03.15 Improved multi-character chain support
-  V2.14 11.06.15 Moved generally useful code into Bioplib
-  V2.15 13.03.19 Now valgrind clean

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "bioplib/array.h"


/************************************************************************/
/* Defines
*/
#define MAXTITLE  160
#define MAXBUFF   160
#define GAPPEN     2   /* 12.03.08 Gap penalty was 10!!!                */
#define MAXCHAINS 160

#define safetoupper(x) ((islower(x))?toupper(x):(x))
#define safetolower(x) ((isupper(x))?tolower(x):(x))

/************************************************************************/
/* Globals
*/
BOOL gQuiet     = FALSE,
     gUpper     = FALSE,
     gFASTA     = FALSE;

char gLabel[blMAXPIRLABEL];

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
void Usage(void);
char *FixSequence(char *seqres,    char *sequence, 
                  char **seqchains, char **atomchains,
                  char **outchains, BOOL IgnoreSEQRES,
                  int nAtomChains);
char *CombineSequence(char *align1, char *align2, int align_len);
int GetPDBChains(PDB *pdb, char *chains);
void PrintNumbering(FILE *out, PDB *pdb, MODRES *modres);
char *strdup(const char *s);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for PDB to PIR conversion

-  10.05.94 Original    By: ACRM
-  22.08.97 Various new flags to handle SEQRES
-  26.08.97 Option to ignore nucleic acids
-  02.10.00 Added -x
-  18.10.00 Added -f
-  08.02.01 Added -n
-  07.03.07 Calls ReadWholePDBAtoms()
-  22.05.09 Added -i option
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  26.08.14 Use renamed macros blPDB2SeqXNoX() and blPDB2SeqX(). By: CTP
-  25.11.14 Initialized seqres  By: ACRM
-  11.06.15 Changed to blGetModresWholePDB() and blGetSeqresWholePDB() 
            By: ACRM
-  13.03.19 Adds 2 rather than 1 to the length of the sequence in
            realloc()
*/
int main(int argc, char **argv)
{
   int  filecount = 0,
      len1,
      nAtomChains;
   FILE *in       = stdin,
        *out      = stdout;
   char *sequence,
        *fixedsequence,
        *seqres = NULL,
        **seqchains,
        **atomchains,
        **outchains,
        title[MAXTITLE];
   PDB  *pdb;
   WHOLEPDB *wpdb;
   BOOL ByChain      = FALSE,
        UseSEQRES    = FALSE,
        SkipX        = FALSE,
        IgnoreSEQRES = FALSE,
        doNucleic    = TRUE,
        DoNumbering  = FALSE;
   MODRES *modres = NULL;
   
   if((outchains = (char **)blArray2D(sizeof(char), MAXCHAINS, blMAXCHAINLABEL))==NULL)
   {
      fprintf(stderr,"Error: No memory for outchains array\n");
      return(1);
   }
   if((seqchains = (char **)blArray2D(sizeof(char), MAXCHAINS, blMAXCHAINLABEL))==NULL)
   {
      fprintf(stderr,"Error: No memory for seqchains array\n");
      return(1);
   }
   
   gLabel[0] = gLabel[blMAXPIRLABEL-1]  = '\0';
   title[0]  = title[MAXTITLE-1] = '\0';
   
   argv++;
   argc--;
   
   /* Parse command line and open files                                 */
   while(argc)
   {
      if(argv[0][0] == '-')    /* Handle switches                       */
      {
         switch(argv[0][1])
         {
         case 'h':
         case 'H':
         case '?':
            Usage();
            return(0);
            break;
         case 'l':
         case 'L':
            argc--;
            argv++;
            strncpy(gLabel,argv[0],blMAXPIRLABEL-1);
            break;
         case 't':
         case 'T':
            argc--;
            argv++;
            strncpy(title,argv[0],MAXTITLE-1);
            break;
         case 'q':
         case 'Q':
            gQuiet = TRUE;
            break;
         case 'i':
         case 'I':
            IgnoreSEQRES = TRUE;
            break;
         case 'c':
         case 'C':
            ByChain = TRUE;
            break;
         case 'x':
         case 'X':
            SkipX = TRUE;
            break;
         case 'u':
         case 'U':
            gUpper = TRUE;
            break;
         case 's':
         case 'S':
            UseSEQRES = TRUE;
            break;
         case 'p':
         case 'P':
            doNucleic = FALSE;
            break;
         case 'f':
         case 'F':
            gFASTA = TRUE;
            ByChain = TRUE;
            break;
         case 'n':
         case 'N':
            DoNumbering = TRUE;
            break;
         default:
            Usage();
            return(1);
         }
      }
      else                     /* Filenames                             */
      {
         switch(++filecount)
         {
         case 1:
            if((in = fopen(argv[0],"r"))==NULL)
            {
               fprintf(stderr,"Error: Unable to open input file: %s\n",
                       argv[0]);
               return(1);
            }
            break;
         case 2:
            if((out = fopen(argv[0],"w"))==NULL)
            {
               fprintf(stderr,"Error: Unable to open output file: %s\n",
                       argv[0]);
               return(1);
            }
            break;
         default:
            Usage();
            return(1);
         }
      }
      argc--;
      argv++;
   }

   /* Read PDB file                                                     */
   if(((wpdb = blReadWholePDBAtoms(in)) == NULL)||(wpdb->pdb==NULL))
   {
      fprintf(stderr,"Error: Unable to read atoms from input file%s%s\n",
              ((gLabel[0])?" Label: ":""),
              ((gLabel[0])?gLabel:""));
      return(1);
   }
   pdb = wpdb->pdb;

   if(UseSEQRES)
   {
      /* Read MODRES and SEQRES records                                 */
      modres = blGetModresWholePDB(wpdb);
      seqres = blGetSeqresAsStringWholePDB(wpdb, seqchains, modres, 
                                           doNucleic);
   }

   /* Extract sequence from PDB linked list                             */
   if((atomchains = blGetPDBChainLabels(pdb, &nAtomChains)) == NULL)
   {
      fprintf(stderr,"Error: No memory for atom chain labels\n");
      return(1);
   }

   /* Convert PDB linked list to a sequence                             */
   if(SkipX)
   {
      if(doNucleic)
      {
         if((sequence = blPDB2SeqXNoX(pdb))==NULL)
         {
            fprintf(stderr,"Error: No memory for sequence data\n");
            return(1);
         }
      }
      else
      {
         if((sequence = blPDBProt2SeqXNoX(pdb))==NULL)
         {
            fprintf(stderr,"Error: No memory for sequence data\n");
            return(1);
         }
      }
   }
   else
   {
      if(doNucleic)
      {
         if((sequence = blPDB2SeqX(pdb))==NULL)
         {
            fprintf(stderr,"Error: No memory for sequence data\n");
            return(1);
         }
      }
      else
      {
         if((sequence = blPDBProt2SeqX(pdb))==NULL)
         {
            fprintf(stderr,"Error: No memory for sequence data\n");
            return(1);
         }
      }
   }
      
   /* Append a * since blPDB2Seq() doesn't do this; note that this 
      will have to change if we fix blPDB2Seq() in future
   */
   len1 = strlen(sequence);
   if((sequence=(char *)realloc((void *)sequence,
                                (len1+2)*sizeof(char)))==NULL)
   {
      fprintf(stderr,"Error: No memory to expand sequence data\n");
      return(1);
   }
   strcat(sequence,"*");

   if(UseSEQRES)
   {
      /* Fiddle with sequences to combine information from SEQRES and ATOM
         records
      */
      if((fixedsequence = FixSequence(seqres,sequence,seqchains,
                                      atomchains,outchains,IgnoreSEQRES,
                                      nAtomChains))
         ==NULL)
         return(1);

      /* Write out the PIR file                                         */
      blWriteOneStringPIR(out,gLabel,title,fixedsequence,outchains,
                          ByChain,gFASTA);
   }
   else
   {
      blWriteOneStringPIR(out,gLabel,title,sequence,atomchains,
                          ByChain,gFASTA);
   }

   if(DoNumbering)
   {
      PrintNumbering(out, pdb, modres);
   }

   return(0);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message.

-  10.05.94 Original    By: ACRM
-  11.02.97 Reformatted message.
-  22.08.97 V2.0
-  26.08.97 V2.1
-  10.08.98 V2.2
-  02.10.00 V2.3
-  18.10.00 V2.4
-  08.03.01 V2.5
-  07.03.07 V2.6
-  12.03.07 V2.7
-  22.05.09 V2.8
-  22.07.14 V2.9 By: CTP
-  07.11.14 V2.11 By: ACRM
-  25.11.14 V2.12
-  10.03.15 V2.13
-  11.06.15 V2.14
-  13.03.19 V2.15
*/
void Usage(void)
{
   fprintf(stderr,"\npdb2pir V2.15 (c) 1994-2019 Dr. Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"\nUsage: pdb2pir [-h][-l label][-t title][-s][-c][-x]\
[-u][-p][-q]\n");
   fprintf(stderr,"               [-f][-n][-i] [infile [outfile]]\n");
   fprintf(stderr,"       -h      This help message\n");
   fprintf(stderr,"       -q      Quiet - no warning messages\n");
   fprintf(stderr,"       -x      Do not include X characters for unknown \
amino acids\n");
   fprintf(stderr,"               Simply skip them instead\n");
   fprintf(stderr,"       -c      Separate header for each chain\n");
   fprintf(stderr,"       -s      Use data from SEQRES records\n");
   fprintf(stderr,"       -i      Ignore SEQRES records where there are \
no ATOMs\n");
   fprintf(stderr,"       -u      All sequence output in upper case\n");
   fprintf(stderr,"       -p      Handle protein only. When combined\n");
   fprintf(stderr,"               with -c DNA/RNA chains skipped\n");
   fprintf(stderr,"       -f      Output FASTA format (implies -c)\n");
   fprintf(stderr,"       -n      Output the numbering\n");
   fprintf(stderr,"       -h      This help message\n");
   fprintf(stderr,"       -l      Specify the PIR label identifier\n");
   fprintf(stderr,"       -t      Specify the PIR title line (use double \
inverted commas\n");
   fprintf(stderr,"               if more than one word)\n");
   fprintf(stderr,"       infile  Input PDB file (stdin if not \
specified)\n");
   fprintf(stderr,"       outfile Output PDB file (stdout if not \
specified)\n\n");
   fprintf(stderr,"Extracts a PIR sequence file from a PDB file.\n\n");
   fprintf(stderr,"Normally just extracts the sequence from the ATOM \
records. By specifying\n");
   fprintf(stderr,"the -s flag, the SEQRES records will also be \
considered. The two sequences\n");
   fprintf(stderr,"will be aligned and the ATOM records will be taken as \
correct, but any\n");
   fprintf(stderr,"additional residues from the SEQRES records will be \
added in lower case\n");
   fprintf(stderr,"(or upper case if the -u flag is given). If -i is \
also added then SEQRES\n");
   fprintf(stderr,"sequences with no matching chain in the ATOM records \
are skipped.\n");

   fprintf(stderr,"\nThe -n option causes the sequence to be output \
again in records of the\n");
   fprintf(stderr,"form:\n");
   fprintf(stderr,"># pos resspec aa\n");
   fprintf(stderr,"\nwhere:\n");
   fprintf(stderr,"pos is the position in the sequence (starting \
from 1)\n");
   blPrintResSpecHelp(stderr);
   fprintf(stderr,"aa is the 1-letter amino acid code.\n");
   fprintf(stderr,"\nNote that when used with -s, only the amino acids \
specified in the\n");
   fprintf(stderr,"ATOM coordinate records will be listed in this \
way.\n");
   fprintf(stderr,"\n");
}


/************************************************************************/
/*>int GetPDBChains(PDB *pdb, char *chains)
   ----------------------------------------
*//**

   Extracts a list of chains from a PDB linked list

-  22.08.97 Original   By: ACRM
-  10.03.15 Checks number of chains
*/
int GetPDBChains(PDB *pdb, char *chains)
{
   PDB  *p;
   char lastchain[8];
   int  nchain = 0;
   
   lastchain[0] = '\0';
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!CHAINMATCH(p->chain, lastchain))
      {
         if(nchain >= MAXCHAINS)
            return(-1);
         
         chains[nchain++] = p->chain[0];
         strcpy(lastchain, p->chain);
      }
   }

   chains[nchain] = '\0';

   return(nchain);
}


/************************************************************************/
/*>char *FixSequence(char *seqres, char *sequence, char **seqchains, 
                     char **atomchains, char **outchains, 
                     BOOL IgnoreSEQRES, int nAtomChains)
   -----------------------------------------------------------------
*//**

   Create a final output sequence by combining the information from the
   ATOM and SEQRES records.

-  21.08.97 Original   By: ACRM
-  22.08.97 Added chain information
-  26.08.97 A couple of bug fixes in initialising allocated memory
            Warning/error messages give the label if specified
-  22.05.09 Added IgnoreSEQRES to ignore chains that are in SEQRES but
            not in ATOM records
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
char *FixSequence(char *seqres, char *sequence, char **seqchains, 
                  char **atomchains, char **outchains, BOOL IgnoreSEQRES,
                  int nAtomChains)
{
   int  i, j, len, len1, len2,
        nchain[2],
        align_len,
        NOutChain = 0,
        ArraySize = 0;
   char *ptr,
        *buffer,
        *outseq = NULL,
        *combseq = NULL,
        *align1,
        *align2,
        **seqs[2];
   BOOL DoneSEQRES[MAXCHAINS],
        DoneATOM[MAXCHAINS],
        DoInit;

   /* Set flags to say we haven't handled the sequences yet             */
   for(i=0; i<MAXCHAINS; i++)
   {
      DoneSEQRES[i] = FALSE;
      DoneATOM[i]   = FALSE;
   }
   
   /* If the sequences and chains are identical just copy one of them
      and return
   */
   if(!strcmp(seqres,sequence) && !strcmp(seqchains[0],atomchains[0])) /* FIXME! */
   {
      for(i=0; i<nAtomChains; i++)
      {
         strcpy(outchains[i], seqchains[i]);
      }
      return(strdup(sequence));
   }

   /* Create a temporary buffer to store a sequence                     */
   len1 = strlen(seqres);
   len2 = strlen(sequence);
   if((buffer = (char *)malloc((1+MAX(len1, len2))
                               * sizeof(char)))==NULL)
   {
      return(NULL);
   }

   /* See how many chains there are and create arrays of char pointers
      to store the separate chains
   */
   nchain[0] = blCountchar(seqres,   '*');        /* SEQRES sequence    */
   if(len1 && seqres[len1-1] != '*')
      nchain[0]++;
   nchain[1] = blCountchar(sequence, '*');        /* ATOM sequence      */
   if(len2 && sequence[len2-1] != '*')
      nchain[1]++;
   for(i=0; i<2; i++)
   {
      if((seqs[i] = (char **)malloc(nchain[i] * sizeof(char *)))==NULL)
         return(NULL);
   }

   /* Transfer the individual chains into the split arrays              */
   for(i=0; i<2; i++)
   {
      strcpy(buffer,((i==0)?seqres:sequence));
      ptr = buffer;
      
      for(j=0; j<nchain[i]; j++)
      {
         TERMAT(ptr,'*');
         len = strlen(ptr);
         if((seqs[i][j] = (char *)malloc((1+len)*sizeof(char)))==NULL)
            return(NULL);
         strcpy(seqs[i][j], ptr);
         seqs[i][j][len] = '\0';

         ptr += strlen(ptr) + 1;
      }
   }

   /* Now align the sequences of the matching chains                    */
   for(i=0; i<nchain[0]; i++)
   {
      for(j=0; j<nchain[1]; j++)
      {
         if(CHAINMATCH(seqchains[i], atomchains[j]))
         {
            DoneSEQRES[i] = TRUE;
            DoneATOM[j]   = TRUE;
            strcpy(outchains[NOutChain++], seqchains[i]);
            
            if(!strcmp(seqs[0][i], seqs[1][j]))
            {
               /* If they are identical, copy to the output array
                  (+2 in the array size for * and \0)
               */
               DoInit = (ArraySize) ? FALSE : TRUE;
               ArraySize += strlen(seqs[0][i]) + 2;
               if((outseq = (char *)realloc(outseq, 
                                            ArraySize*sizeof(char)))
                  ==NULL)
               {
                  return(NULL);
               }
               if(DoInit)
                  outseq[0] = '\0';
               
               strcat(outseq,seqs[0][i]);
               strcat(outseq,"*");
            }
            else
            {
               /* The sequences are non-identical so we align them      */
               len1 = strlen(seqs[0][i]);
               len2 = strlen(seqs[1][j]);
               if((align1=(char *)malloc((len1+len2)*sizeof(char)))==NULL)
                  return(NULL);
               if((align2=(char *)malloc((len1+len2)*sizeof(char)))==NULL)
                  return(NULL);
            
               if(!blAlign(seqs[0][i], len1,
                           seqs[1][j], len2,
                           FALSE, TRUE, GAPPEN,
                           align1, align2, &align_len))
                  return(NULL);

               if((combseq=CombineSequence(align1,align2,align_len))
                  ==NULL)
                  return(NULL);

               free(align1); 
               free(align2);

               /* Allocate memory for the output sequence
                  (+2 in the array size for * and \0)
               */
               if(ArraySize==0)
               {
                  ArraySize = align_len + 2;
                  if((outseq = (char *)malloc(ArraySize*sizeof(char)))
                     ==NULL)
                     return(NULL);
                  outseq[0] = '\0';
               }
               else
               {
                  ArraySize += align_len + 2;
                  if((outseq = (char *)realloc(outseq,
                                               ArraySize*sizeof(char)))
                     ==NULL)
                     return(NULL);
                  outseq[ArraySize-1] = '\0';
               }

               strcat(outseq,combseq);
               strcat(outseq,"*");
               free(combseq);
            }
            break;
         }
      }
   }

   /* Add any chains from the ATOM records not yet handled              */
   for(i=0; i<nchain[1]; i++)
   {
      if(!DoneATOM[i])
      {
         /* Allocate memory for the output sequence
            (+2 in the array size for * and \0)
         */
         if(ArraySize==0)
         {
            ArraySize = strlen(seqs[1][i]) + 2;
            if((outseq = (char *)malloc(ArraySize*sizeof(char)))
               ==NULL)
               return(NULL);
            outseq[0] = '\0';
         }
         else
         {
            ArraySize += strlen(seqs[1][i]) + 2;
            if((outseq = (char *)realloc(outseq,
                                         ArraySize*sizeof(char)))
               ==NULL)
               return(NULL);
            outseq[ArraySize-1] = '\0';
         }
         
         strcat(outseq,seqs[1][i]);
         strcat(outseq,"*");
         strcpy(outchains[NOutChain++], atomchains[i]);
      }
   }

   /* Add any chains from the SEQRES records not yet handled.
      We issue a warning about these ones
   */
   if(!IgnoreSEQRES) /* 22.05.09 Added this check                       */
   {
      for(i=0; i<nchain[0]; i++)
      {
         if(!DoneSEQRES[i])
         {
            /* Allocate memory for the output sequence
               (+2 in the array size for * and \0)
            */
            if(ArraySize==0)
            {
               ArraySize = strlen(seqs[0][i]) + 2;
               if((outseq = (char *)malloc(ArraySize*sizeof(char)))
                  ==NULL)
                  return(NULL);
               outseq[0] = '\0';
            }
            else
            {
               ArraySize += strlen(seqs[0][i]) + 2;
               if((outseq = (char *)realloc(outseq,
                                            ArraySize*sizeof(char)))
                  ==NULL)
                  return(NULL);
               outseq[ArraySize-1] = '\0';
            }

            /* 22.05.09 Added this                                      */
            if(!gUpper)
            {
               LOWER(seqs[0][i]);
            }
            
            strcat(outseq,seqs[0][i]);
            strcat(outseq,"*");
            strcpy(outchains[NOutChain++], seqchains[i]);
            
            if(!gQuiet)
            {
               fprintf(stderr,"Warning: Chain %s from SEQRES records \
not found in ATOM records%s%s\n",
                       seqchains[i],
                       ((gLabel[0])?" Label: ":""),
                       ((gLabel[0])?gLabel:""));
            }
         }
      }
   }

   /*** NEEDS TO FREE seqs[][] too                                    ***/
   free(buffer);         /*  11.06.15                                   */
   return(outseq);
}


/************************************************************************/
/*>char *CombineSequence(char *align1, char *align2, int align_len)
   ----------------------------------------------------------------
*//**

   Combine the information from the two sequences

-  22.08.97 Original   By: ACRM
*/
char *CombineSequence(char *align1, char *align2, int align_len)
{
   static char *outseq = NULL;
   int         i;

   if((outseq=(char *)malloc((align_len+1) + sizeof(char)))==NULL)
      return(NULL);

#ifdef DEBUG
   align1[align_len] = '\0';
   fprintf(stderr,"%s\n", align1);
   align2[align_len] = '\0';
   fprintf(stderr,"%s\n", align2);
#endif
   
   
   for(i=0; i<align_len; i++)
   {
      if((align1[i] == align2[i]) || (align1[i] == '-'))
      {
         outseq[i] = safetoupper(align2[i]);
      }
      else
      {
         if(gUpper)
            outseq[i] = safetoupper(align1[i]);
         else
            outseq[i] = safetolower(align1[i]);
      }
   }
   outseq[align_len] = '\0';
   
   return(outseq);
}




/************************************************************************/
/*>void PrintNumbering(FILE *out, PDB *pdb, MODRES *modres)
   --------------------------------------------------------
*//**

   \param[in]      *out       File to which to send output
   \param[in]      *pdb       PDB linked list
   \param[in]      *modres    MODRES linked list

   Simply works through the PDB linked list printing lines of the form
   ># pos c.nnni aa
   where pos is the position in the sequence (all chains numbered through
   sequentially), c.nnni is the chain/resnum/insert code and aa is the
   1-letter amino acid code.

-  08.03.01 Original   By: ACRM
-  07.03.07 Added modres stuff
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  10.03.15 Updated for multicharacter and numeric chain labels By: ACRM
-  11.06.15 Changed to blFindOriginalResType()
*/
void PrintNumbering(FILE *out, PDB *pdb, MODRES *modres)
{
   PDB  *p;
   int  pos = 0;
   int  lastresnum = -999999;
   char lastchain[8],
        lastinsert = ' ',
        resid[16],
        one;
   
   lastchain[0] = '\0';
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum    != lastresnum) ||
         (p->insert[0] != lastinsert) ||
         !CHAINMATCH(p->chain, lastchain))
      {
         pos++;
         lastresnum = p->resnum;
         strcpy(lastchain, p->chain);
         lastinsert = p->insert[0];
         
         sprintf(resid,"%s.%d%c", p->chain, p->resnum, p->insert[0]);
         one = blThrone(p->resnam);

         /* 07.03.07 Added code to check for modified amino acids    */
         if(one == 'X')
         {
            char tmpthree[8];
            blFindOriginalResType(p->resnam, tmpthree, modres);
            one = blThrone(tmpthree);
         }
               
         fprintf(out, "># %d %s %c\n", pos, resid, one);
      }
   }
}

