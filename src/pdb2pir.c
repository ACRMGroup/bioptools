/************************************************************************/
/**

   \file       pdb2pir.c
   
   \version    V2.10
   \date       26.08.14
   \brief      Convert PDB to PIR sequence file
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL 1994-2014
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
#include "bioplib/fsscanf.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines
*/
#define MAXLAB     64
#define MAXTITLE  160
#define ALLOCSIZE  80
#define MAXBUFF   160
#define GAPPEN     2   /* 12.03.08 Gap penalty was 10!!!                */
#define MAXCHAINS  80

#define safetoupper(x) ((islower(x))?toupper(x):(x))
#define safetolower(x) ((isupper(x))?tolower(x):(x))

typedef struct mres
{
   struct mres *next;
   char modres[8];
   char origres[8];
}  MODRES;

/************************************************************************/
/* Globals
*/
BOOL gQuiet     = FALSE,
     gUpper     = FALSE,
     gDoNucleic = TRUE,
     gFASTA     = FALSE;

char gLabel[MAXLAB];

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
void WritePIR(FILE *out, char *label, char *title, char *sequence,
              char *chains, BOOL ByChain);
void Usage(void);
char *ReadSEQRES(WHOLEPDB *wpdb, char *chains, MODRES *modres);
MODRES *ReadMODRES(WHOLEPDB *wpdb);
char *FixSequence(char *seqres,    char *sequence, 
                  char *seqchains, char *atomchains,
                  char *outchains, BOOL IgnoreSEQRES);
char *CombineSequence(char *align1, char *align2, int align_len);
int GetPDBChains(PDB *pdb, char *chains);
void PrintNumbering(FILE *out, PDB *pdb, MODRES *modres);
char *strdup(const char *s);
void LookupModres(char *orig, char *new, MODRES *modres);


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
*/
int main(int argc, char **argv)
{
   int  filecount = 0,
        len1;
   FILE *in       = stdin,
        *out      = stdout;
   char *sequence,
        *fixedsequence,
        *seqres,
        seqchains[MAXCHAINS],
        atomchains[MAXCHAINS],
        outchains[MAXCHAINS],
        title[MAXTITLE];
   PDB  *pdb;
   WHOLEPDB *wpdb;
   BOOL ByChain      = FALSE,
        UseSEQRES    = FALSE,
        SkipX        = FALSE,
        IgnoreSEQRES = FALSE,
        DoNumbering  = FALSE;
   MODRES *modres = NULL;
   

   gLabel[0] = gLabel[MAXLAB-1]  = '\0';
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
            strncpy(gLabel,argv[0],MAXLAB-1);
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
            gDoNucleic = FALSE;
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
   if((wpdb = blReadWholePDBAtoms(in)) == NULL)
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
      modres = ReadMODRES(wpdb);
      seqres = ReadSEQRES(wpdb, seqchains, modres);
   }

   /* Extract sequence from PDB linked list                             */
   GetPDBChains(pdb, atomchains);

   /* Convert PDB linked list to a sequence                             */
   if(SkipX)
   {
      if(gDoNucleic)
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
      if(gDoNucleic)
      {
         if((sequence = PDB2SeqX(pdb))==NULL)
         {
            fprintf(stderr,"Error: No memory for sequence data\n");
            return(1);
         }
      }
      else
      {
         if((sequence = PDBProt2SeqX(pdb))==NULL)
         {
            fprintf(stderr,"Error: No memory for sequence data\n");
            return(1);
         }
      }
   }
      
   /* Append a * since PDB2Seq() doesn't do this; note that this will
      have to change if we fix PDB2Seq() in future
   */
   len1 = strlen(sequence);
   if((sequence=(char *)realloc((void *)sequence,
                                (len1+1)*sizeof(char)))==NULL)
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
                                      atomchains,outchains,IgnoreSEQRES))
         ==NULL)
         return(1);

      /* Write out the PIR file                                         */
      WritePIR(out,gLabel,title,fixedsequence,outchains,ByChain);
   }
   else
   {
      WritePIR(out,gLabel,title,sequence,atomchains,ByChain);
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
*/
void Usage(void)
{
   fprintf(stderr,"\npdb2pir V2.9 (c) 1994-2014 Dr. Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"\nUsage: pdb2pir [-h] [-l label] [-t title] [-s] [-c] \
[-u] [-p] [-q] [-x] [-f] [-n] [-i] [infile [outfile]]\n");
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
also added then SEQRES sequences\n");
   fprintf(stderr,"with no matching chain in the ATOM records are \
skipped.\n");

   fprintf(stderr,"\nThe -n option causes the sequence to be output \
again in records of the\n");
   fprintf(stderr,"form:\n");
   fprintf(stderr,"># pos resnum aa\n");
   fprintf(stderr,"where pos is the position in the sequence (starting \
from 1), resnum\n");
   fprintf(stderr,"is the residue number in the PDB file in the form \
[c]nnn[i] (where\n");
   fprintf(stderr,"c is an optional chain label, nnn is the residue \
number and i is \n");
   fprintf(stderr,"an optional insert code) and aa is the 1-letter \
amino acid code.\n");
   fprintf(stderr,"Note that when used with -s, only the amino acids \
specified in the\n");
   fprintf(stderr,"ATOM coordinate records will be listed in this \
way.\n");
   fprintf(stderr,"\n");
}


/************************************************************************/
/*>char *ReadSEQRES(WHOLEPDB *wpdb, char *chains, MODRES *modres)
   --------------------------------------------------------------
*//**

   Reads sequence from SEQRES records into a character string in 1-letter
   code. Chains are terminated by * characters.

-  21.08.97 Original   by: ACRM
-  22.08.97 Added chains parameter
-  26.08.97 No longer reads DNA/RNA
-  07.03.07 Added code to check for modified amino acids
            Now reads from wpdb rather than from the file
*/
char *ReadSEQRES(WHOLEPDB *wpdb, char *chains, MODRES *modres)
{
   static char *sequence = NULL;
   char        buffer[MAXBUFF],
               chain,
               lastchain,
               seq3[13][4];
   int         i,
               nchain    = 0,
               nres      = 0,
               ArraySize = ALLOCSIZE;
   BOOL        AddStar   = FALSE;
   STRINGLIST  *s;
   

   if((sequence=(char *)malloc(ArraySize * sizeof(char)))==NULL)
   {
      return(NULL);
   }
   sequence[0] = '\0';
   
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      strncpy(buffer, s->string, MAXBUFF);
      TERMINATE(buffer);
      if(!strncmp(buffer,"SEQRES",6))
      {
         fsscanf(buffer, 
                 "%11x%c%7x%3s%1x%3s%1x%3s%1x%3s%1x%3s%1x%3s%1x%3s%1x\
%3s%1x%3s%1x%3s%1x%3s%1x%3s%1x%3s",
                 &chain,
                 seq3[0],  seq3[1],  seq3[2],  seq3[3],  seq3[4], 
                 seq3[5],  seq3[6],  seq3[7],  seq3[8],  seq3[9],
                 seq3[10], seq3[11], seq3[12]);

         if((nres == 0) && !AddStar)
         {
            /* This is the first line so we set the lastchain           */
            lastchain = chain;
            if(chains!=NULL)
               chains[nchain++] = chain;
         }
         else if(nres+15 >= ArraySize)
         {
            /* Allocate more space if needed                            */
            ArraySize += ALLOCSIZE;
            if((sequence=(char *)realloc((void *)sequence, 
                                         ArraySize*sizeof(char)))
               == NULL)
            {
               return(NULL);
            }
         }

         if(chain != lastchain)
         {
            sequence[nres++] = '*';
            lastchain = chain;
            if(chains!=NULL)
               chains[nchain++] = chain;
         }

         for(i=0; i<13; i++)
         {
            AddStar = TRUE;
            if(!strncmp(seq3[i],"   ",3))
               break;
            sequence[nres] = blThronex(seq3[i]);

            /* 07.03.07 Added code to check for modified amino acids    */
            if(sequence[nres] == 'X')
            {
               char tmpthree[8];
               LookupModres(seq3[i], tmpthree, modres);
               sequence[nres] = blThronex(tmpthree);
            }
               
            if(!gBioplibSeqNucleicAcid || gDoNucleic)
               nres++;
         }
      }
   }

   if(AddStar)
   {
      sequence[nres++] = '*';
   }
   sequence[nres++] = '\0';
   if(chains!=NULL)
      chains[nchain] = '\0';

   return(sequence);
}

/************************************************************************/
/*>int GetPDBChains(PDB *pdb, char *chains)
   ----------------------------------------
*//**

   Extracts a list of chains from a PDB linked list

-  22.08.97 Original   By: ACRM
*/
int GetPDBChains(PDB *pdb, char *chains)
{
   PDB  *p;
   char lastchain = '\0';
   int  nchain = 0;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->chain[0] != lastchain)
      {
         chains[nchain++] = p->chain[0];
         lastchain = p->chain[0];
      }
   }

   chains[nchain] = '\0';

   return(nchain);
}


/************************************************************************/
/*>char *FixSequence(char *seqres, char *sequence, char *seqchains, 
                     char *atomchains, char *outchains, BOOL IgnoreSEQRES)
   -----------------------------------------------------------------------
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
char *FixSequence(char *seqres, char *sequence, char *seqchains, 
                  char *atomchains, char *outchains, BOOL IgnoreSEQRES)
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
   if(!strcmp(seqres,sequence) && !strcmp(seqchains,atomchains))
   {
      strcpy(outchains, seqchains);
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
         if(seqchains[i] == atomchains[j])
         {
            DoneSEQRES[i] = TRUE;
            DoneATOM[j]   = TRUE;
            outchains[NOutChain++] = seqchains[i];
            
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
         outchains[NOutChain++] = atomchains[i];
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
            outchains[NOutChain++] = seqchains[i];
            
            if(!gQuiet)
            {
               fprintf(stderr,"Warning: Chain %c from SEQRES records \
not found in ATOM records%s%s\n",
                       seqchains[i],
                       ((gLabel[0])?" Label: ":""),
                       ((gLabel[0])?gLabel:""));
            }
         }
      }
   }

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
/*>void WritePIR(FILE *out, char *label, char *title, char *sequence,
                 char *chains, BOOL ByChain)
   ------------------------------------------------------------------
*//**

   \param[in]      *out        File pointer
   \param[in]      *sequence   1-letter amino acid code

   Writes a PIR sequence file from a 1-letter code array.
   Adds a terminating * if required.

-  10.05.94 Original    By: ACRM
-  22.08.97 Can now handle chains separately
-  26.08.97 If chains are handled seprately, don't bother writing out
            an empty chain
-  10.08.98 Basically a total rewrite to fix a bug which caused the 
            header not to be printed with -c -p for a chain after one
            which was non-protein. Much simplified the code by printing
            the header at the beginning of a chain rather than end
            of previous chain.
-  18.10.00 Added code to write FASTA as well
*/
void WritePIR(FILE *out, char *label, char *title, char *sequence,
              char *chains, BOOL ByChain)
{
   int  i, 
        count,
        chaincount = 0;
   BOOL GotStar = FALSE,
        DoneHeader = FALSE,
        Printed;
   char chn[8],
        outlabel[MAXLAB];

   strcpy(outlabel,label);

   /* If we are not going by chain create a label                       */
   if(!ByChain)
   {
      strcpy(outlabel, "PDBPIR");
      if(label[0])
         strcpy(outlabel, label);
   }

   /* If we are not doing it by-chain, then simple display the label
      and header now otherwise set flag to say we haven't done header
   */
   if(ByChain)
   {
      DoneHeader = FALSE;
   }
   else
   {
      if(gFASTA)
      {
         fprintf(out,">%s\n",outlabel);
      }
      else
      {
         fprintf(out,">P1;%s\n",outlabel);
         fprintf(out,"Sequence extracted from PDB file - %s\n",
                 (title[0]?title:"By pdb2pir"));
      }
      DoneHeader = TRUE;
   }

   /* Start by setting flag to say nothing has been printed             */
   Printed = FALSE;
   
   /* Loop through the sequence with i. Use count to count number of
      residues printed on a line
   */
   for(i=0,count=1; sequence[i]; i++)
   {
      if(ByChain)
      {
         /* If we don't have a header, then print one                   */
         if((sequence[i] != '*') && (!DoneHeader))
         {
            if(label[0])
            {
               strcpy(outlabel,label);
               chn[0] = chains[chaincount];
               chn[1] = '\0';
               strcat(outlabel,chn);
            }
            else
            {
               sprintf(outlabel,"Chain%c",chains[chaincount]);
            }
            if(gFASTA)
            {
               fprintf(out,">%s\n",outlabel);
            }
            else
            {
               fprintf(out,">P1;%s\n",outlabel);
               fprintf(out,"Sequence extracted from PDB file - %s\n",
                       (title[0]?title:"By pdb2pir"));
            }
            
            DoneHeader = TRUE;
            count = 0;
         }
      }
      
      if(count==30)
      {
         fputc('\n',out);
         count = 0;
      }

      if(sequence[i] == '*')
      {
         chaincount++;
         if(Printed)
         {
            if(gFASTA)
            {
               fprintf(out,"\n");
            }
            else
            {
               fprintf(out,"*\n");
            }
            
            GotStar = TRUE;
         }
         
         DoneHeader = FALSE;
         Printed    = FALSE;
      }
      else
      {
         fputc(sequence[i],out);
         GotStar = FALSE;
         Printed = TRUE;
      }

      if(Printed)
         count++;
   }

   if(!GotStar)
   {
      if(gFASTA)
      {
         fprintf(out,"\n");
      }
      else
      {
         fprintf(out,"*\n");
      }
   }
}


/************************************************************************/
/*>void PrintNumbering(FILE *out, PDB *pdb, MODRES *modres)
   --------------------------------------------------------
*//**

   \param[in]      *out       File to which to send output
   \param[in]      *pdb       PDB linked list
   \param[in]      *modres    MODRES linked list

   Simply works through the PDB linked list printing lines of the form
   ># pos cnnni aa
   where pos is the position in the sequence (all chains numbered through
   sequentially), cnnni is the chain/resnum/insert code and aa is the
   1-letter amino acid code.

-  08.03.01 Original   By: ACRM
-  07.03.07 Added modres stuff
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
void PrintNumbering(FILE *out, PDB *pdb, MODRES *modres)
{
   PDB  *p;
   int  pos = 0;
   int  lastresnum = -999999;
   char lastchain  = ' ',
        lastinsert = ' ',
        resid[16],
        one;
   
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->resnum    != lastresnum) ||
         (p->insert[0] != lastinsert) ||
         (p->chain[0]  != lastchain))
      {
         pos++;
         lastresnum = p->resnum;
         lastchain  = p->chain[0];
         lastinsert = p->insert[0];
         
         sprintf(resid,"%c%d%c", p->chain[0], p->resnum, p->insert[0]);
         one = blThrone(p->resnam);

         /* 07.03.07 Added code to check for modified amino acids    */
         if(one == 'X')
         {
            char tmpthree[8];
            LookupModres(p->resnam, tmpthree, modres);
            one = blThrone(tmpthree);
         }
               
         fprintf(out, "># %d %s %c\n", pos, resid, one);
      }
   }
}

void LookupModres(char *orig, char *new, MODRES *modres)
{
   MODRES *m;
   for(m=modres; m!=NULL; NEXT(m))
   {
      if(!strncmp(orig, m->modres, 3))
      {
         strncpy(new, m->origres, 3);
         PADCHARMINTERM(new, ' ', 4);
         return;
      }
   }
}

MODRES *ReadMODRES(WHOLEPDB *wpdb)
{
   STRINGLIST *s;
   char *ch;
   MODRES *modres = NULL,
          *m = NULL;
   
   
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, "MODRES", 6))
      {
         if(m==NULL)
         {
            INIT(modres, MODRES);
            m = modres;
         }
         else
         {
            ALLOCNEXT(m, MODRES);
         }
         if(m==NULL)
         {
            fprintf(stderr,"pdb2pir: Error! No memory for modres\n");
            exit(1);
         }
         
         ch = s->string+12;
         strncpy(m->modres, ch, 3);
         PADCHARMINTERM(m->modres, ' ', 4);
         
         ch = s->string+24;
         strncpy(m->origres, ch, 3);
         PADCHARMINTERM(m->origres, ' ', 4);
         if(m->origres[0] == ' ')
         {
            strncpy(m->origres, "XXX ", 4);
         }
      }
   }
   return(modres);
}


