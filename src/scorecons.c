/*************************************************************************

   Program:    scorecons
   File:       scorecons.c
   
   Version:    V1.9
   Date:       13.03.25
   Function:   Scores conservation from a PIR sequence alignment
               Not to be confused with the program of the same name
               by Will Valdar (this one predates his!)
   
   Copyright:  (c) Prof. Andrew C. R. Martin 1996-2025
   Author:     Prof. Andrew C. R. Martin
               Tom Northey (implemented Valdar01 scoring)
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
   Calculates a conservation score from a PIR alignment file.
   Two major methods are implemented:

   1. A matrix-based score which is calculated as:

           ---       /
       S = \ D      / ij
           /  ij   /
           ---    /
            ij
    By default, this uses the pet91 matrix, the Jones et al (1991) update 
    of the Dayhoff matrix, normalised such that all scores on the diagonal
    are equal.

    2. A statistical-entropy based score:

              N
             ---          /
       S = - \ p log p   / log(MIN(NSeq, N))
             /  i     i /
             ---       /
             i=1

    where N is the number of groups (i.e. 20 amino acids) and NSeq is the
    number of sequences being compared. The division by log(MIN(NSeq, N))
    simply ensures that the entropy is in the range 0-1
    To turn this value into a conservation score, the result is subtracted
    from 1 such that 1 is conserved and 0 is absolutely non-conserved

    Two version of this score are provided; one where the 20 amino acids
    are all counted as different; the other where they are classes into
    8 groups.

    In addition, a combined score may be calculated where:

              --                        --
              |                          |
    S = S   x | (1 - 8/20) x S   +  8/20 |
         20   |               8          |
              --                        --

    This equation has a number of favourable properties:
    a. If S20 == S8 == 1, then S = 1
    b. If S8 == 0, then S = S20 x 8/20
    c. As the information content of the reduced set decreases (i.e.
       we use something like S3 instead of S8), so the influence of
       the reduced set on the value of S decreases.

    Also considered:

            --------------------------
           / S'  x  (S' + (20/8 - 1))
          /   20      8
    S' = /   ------------------------
        /           (20/8)
      \/
            (where S' = (1.0-S) )



         S     +  (8/20) x S
          20                8
    S =  --------------------
              1 + (8/20)

   3. The Valdar01 scoring method

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   11.09.96 Original
   V1.1   17.09.96 Tidied up and improved command line options
                   Added entropy scoring
   V1.2   18.09.96 Corrected calcuation of combined entropy score
   V1.3   15.07.08 Added -x flag
   V1.4   11.08.15 Modified for new bioplib
   V1.5   24.08.15 Implemented the Valdar01 scoring By: TCN
   V1.6   11.10.19 Fixed reading of a specified matrix - it was ignoring
                   -m before!
   V1.7   10.08.22 Added -s option to score a single position's 
                   distribution of amino acids
   V1.7.1 13.09.22 Modified int to LONG for very big data setps.
   V1.8   04.10.22 Added -l (log) option
   V1.8.1 22.11.22 Fixed bug in command parsing - was failing when no
                   flags given.
   V1.9   13.03.25 Added -i flag to skip deletions

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/seq.h"
#include "bioplib/array.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
typedef struct _seqdata
{
   struct _seqdata *next;
   char *seqs[8];
}  SEQDATA;

typedef struct 
{
   char res;
   int  NGroup,
        group[2];
}  AMINOACID;

typedef long long int VLONG;
#define VLONGMAX 9223372036854775000
#define LONGMAX  2147483600
#define MAXCHECK VLONGMAX
typedef long double LREAL;   

#define DATADIR "DATADIR"
#define MUTMAT  "pet91.mat"
#define MAXBUFF 160
#define MAXDATA 1000
#define LOGSCALE 100
#define TINY    0.000001
#define METH_MDM       0
#define METH_ENTROPY   1
#define METH_ENTROPY20 2
#define METH_ENTROPY8  3
#define METH_VALDAR    4

#define MINSINLEN      8

#define LAMBDA_UNITIALIZED 0
#define SEQWEIGHTS_UNITIALIZED NULL

/************************************************************************/
/* Globals
*/

static REAL sLambda = LAMBDA_UNITIALIZED;
static REAL *sSeqWeights=NULL;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *matrix, int *Method, BOOL *extended,
                  char *single, BOOL *doLog, REAL *maxFraction,
                  BOOL *reduceData, BOOL *ignoreGaps);
void Usage(void);
BOOL ReadAndScoreSeqs(FILE *fp, FILE *out, int MaxInMatrix, int Method,
                      BOOL Extended, BOOL ignoreGaps);
SEQDATA *ReadAllSeqs(FILE *fp);
char **ListToTable(SEQDATA *SeqList, int *nseq, int *seqlen);
void PadSeqs(char **SeqTable, int nseq, int seqlen);
void DisplayScores(FILE *fp, char **SeqTable, int nseq, int seqlen, 
                   int MaxInMatrix, int Method, BOOL Extended,
                   BOOL ignoreGaps);
REAL CalcScore(char **SeqTable, int nseq, int seqlen, int pos,
               int MaxInMatrix, int Method, BOOL ignoreGaps);
REAL MDMBasedScore(char **SeqTable, int nseq, int pos, int MaxInMatrix,
                   BOOL ignoreGaps);
REAL EntropyScore(char **SeqTable, int nseq, int pos, 
                  AMINOACID *aminoacids, int NGroups);

REAL valdarScore(char **SeqTable, int pos, int numSeqs, int seqlen, 
                 int MaxInMatrix);
void initLambda(int numSeqs);
BOOL initSequenceWeights(char **SeqTable, int numSeqs, int seqlen, 
                         int MaxInMatrix);
REAL **getSeqDistTable(char **SeqTable, int numSeqs, int seqlen, 
                       int MaxInMatrix);
REAL getInterSeqDistance(char **SeqTable, int seqAindex, int seqBindex, 
                         int seqlen, int MaxInMatrix);
int getNonGapPosCount(char **SeqTable, int seqAindex, int seqBindex, 
                      int seqlen);
REAL valdarMatrixScore(char res1, char res2, int MaxInMatrix);
BOOL ReadAndScoreSingle(char *single, FILE *out, int MaxInMatrix,
                        int Method, BOOL Extended, BOOL doLog,
                        REAL maxFraction, BOOL reduceData);
char **ParseSingle(char *single, int *nseq, BOOL doLog, REAL maxFraction,
                   BOOL reduceData);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for scoring variability on the basis of a mutation matrix
   across a sequence alignment stored in a PDB file.

   11.09.96 Original   By: ACRM
   17.09.96 Rewritten. Zeros the MDM
   18.09.96 Added check on environment variable if ReadMDM() failed.
   15.07.08 Added -x/Extended handling
   11.10.19 Fixed code to actually read a different matrix if specified!
   10.08.22 Added code to handle single column mode
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char InFile[MAXBUFF],
        OutFile[MAXBUFF],
        matrix[MAXBUFF],
        single[MAXBUFF];
   int  MaxInMatrix,
        Method      = METH_MDM;
   BOOL Extended    = FALSE,
        doLog       = FALSE,
        ignoreGaps  = FALSE,
        reduceData  = FALSE;
   REAL maxFraction = (REAL)0.0;
   

   strncpy(matrix, MUTMAT, MAXBUFF-1);
   
   if(ParseCmdLine(argc, argv, InFile, OutFile, matrix, &Method,
                   &Extended, single, &doLog, &maxFraction,
                   &reduceData, &ignoreGaps))
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         if(!blReadMDM(matrix))
         {
            fprintf(stderr,"Unable to read mutation matrix: %s\n",matrix);
            if(getenv(DATADIR) == NULL)
            {
               fprintf(stderr,"Environment variable (%s) not set.\n",
                       DATADIR);
            }
            return(1);
         }
         MaxInMatrix = blZeroMDM();
         if(single[0] != '\0')
         {
            return(ReadAndScoreSingle(single, out, MaxInMatrix, Method,
                                      Extended, doLog, maxFraction,
                                      reduceData)?0:1);
         }
         else
         {
            return(ReadAndScoreSeqs(in, out, MaxInMatrix, Method,
                                    Extended, ignoreGaps)?0:1);
         }
      }
      else
      {
         return(1);
      }
   }
   else
   {
      Usage();
   }
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     char *matrix, int *Method, BOOL *Extended,
                     char *single, BOOL *doLog, REAL *maxFraction,
                     BOOL *reduceData, BOOL *ignoreGaps)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            char   *matrix      Mutation matrix name
            int    *Method      Scoring method
            BOOL   *Extended    Extended precision printing
            char   *single      Single position distribution (-s)
            BOOL   *doLog       -l Used with -s to take logs of all the
                                counts
            REAL   *maxFraction First count is scaled down to be no 
                                more than the specified fraction of the
                                data
            BOOL   *reduceData  Reduce the dataset size to MAXDATA
            BOOL   *ignoreGaps  Ignore gaps in the alignment scoring
   Returns: BOOL                Success?

   Parse the command line
   
   17.09.96 Original    By: ACRM
   15.07.08 Added -x
   24.08.15 Added -d    By: TCN
   10.08.22 Added -s    By: ACRM
   04.10.22 Added -l, -f, -r
   22.11.22 Fixed bug when no flags given
   13.03.25 Added -i
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *matrix, int *Method, BOOL *Extended,
                  char *single, BOOL *doLog, REAL *maxFraction,
                  BOOL *reduceData, BOOL *ignoreGaps)
{
   int nFlags = 0;
   
   argc--;
   argv++;

   strncpy(matrix, MUTMAT, MAXBUFF);
   infile[0]    = outfile[0] = '\0';
   single[0]    = '\0';
   *Extended    = FALSE;
   *reduceData  = FALSE;
   *ignoreGaps  = FALSE;
   *maxFraction = (REAL)0.0;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'm':
            if(--argc < 0)
               return(FALSE);
            argv++;
            strncpy(matrix,argv[0],MAXBUFF);
            break;
         case 'e':
            *Method = METH_ENTROPY;
            break;
         case 'a':
            *Method = METH_ENTROPY20;
            break;
         case 'g':
            *Method = METH_ENTROPY8;
            break;
         case 'd':
            *Method = METH_VALDAR;
            break;
         case 'x':
            *Extended = TRUE;
            break;
         case 'i':
            *ignoreGaps = TRUE;
            break;
         case 'l':
            *doLog = TRUE;
            nFlags++;
            break;
         case 'r':
            *reduceData = TRUE;
            nFlags++;
            break;
         case 's':
            if(--argc < 0)
               return(FALSE);
            argv++;
            strncpy(single,argv[0],MAXBUFF);
            break;
         case 'f':
            nFlags++;
            *maxFraction = (REAL)0.5;
            if(argv[0][2] == '=')
            {
               if(!sscanf(((argv[0])+3), "%lf", maxFraction))
                  return(FALSE);
               if(*maxFraction >= (REAL)1.0)
                  return(FALSE);
            }
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         if(single[0] != '\0')
         {
            /* Check there is only 1 argument left                      */
            if(argc > 1)
               return(FALSE);

            strcpy(outfile, argv[0]);
         }
         else
         {
            if(*doLog)
            {
               fprintf(stderr,"Error: -l must be used with -s\n");
               return(FALSE);
            }
            if(*reduceData)
            {
               fprintf(stderr,"Error: -r must be used with -s\n");
               return(FALSE);
            }
            if(*maxFraction > TINY)
            {
               fprintf(stderr,"Error: -f must be used with -s\n");
               return(FALSE);
            }
            
            /* Check that there are only 1 or 2 arguments left          */
            if(argc > 2)
               return(FALSE);
            
            /* Copy the first to infile                                 */
            strcpy(infile, argv[0]);
            
            /* If there's another, copy it to outfile                   */
            argc--;
            argv++;
            if(argc)
               strcpy(outfile, argv[0]);
         }

         if(nFlags > 1)
            return(FALSE);
         
         return(TRUE);
      }
      argc--;
      argv++;
   }

   if(nFlags > 1)
      return(FALSE);
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL ReadAndScoreSeqs(FILE *fp, FILE *out, int MaxInMatrix, int Method,
                         BOOL Extended, BOOL ignoreGaps)
   -----------------------------------------------------------------------
   Routine which reads in files, calculates and displays the variability
   scores.

   11.09.96 Original   By: ACRM
   17.09.96 Added out parameter and MaxInMatrix
   15.07.08 Added Extended parameter
   13.03.25 Added ignoreGaps parameter
*/
BOOL ReadAndScoreSeqs(FILE *fp, FILE *out, int MaxInMatrix, int Method,
                      BOOL Extended, BOOL ignoreGaps)
{
   SEQDATA *SeqList;
   char    **SeqTable;
   int     seqlen, nseq;
   
   /* Read sequences into a linked list                                 */
   if((SeqList=ReadAllSeqs(fp))==NULL)
      return(FALSE);

   /* Reorganise into a 2D table                                        */
   if((SeqTable = ListToTable(SeqList, &nseq, &seqlen))==NULL)
   {
      FREELIST(SeqList, SEQDATA);
      return(FALSE);
   }

   /* Free up linked list                                               */
   FREELIST(SeqList, SEQDATA);

   /* Calculate and print scores                                        */
   DisplayScores(out, SeqTable, nseq, seqlen, MaxInMatrix, Method,
                 Extended, ignoreGaps);

   /* Free memory from sequence table                                   */
   blFreeArray2D((char **)SeqTable, nseq, seqlen);

   return(TRUE);
}


/************************************************************************/
/*>SEQDATA *ReadAllSeqs(FILE *fp)
   ------------------------------
   Reads in all the sequences in a PIR file into a linked list of 
   sequences.

   11.09.96 Original   By: ACRM
*/
SEQDATA *ReadAllSeqs(FILE *fp)
{
   SEQDATA *start = NULL, 
           *prev  = NULL,
           *p     = NULL;
   SEQINFO seqinfo;
   BOOL    punct,
           error;

   INIT(start, SEQDATA);
   if((p=start)==NULL)
      return(NULL);
   
   while(blReadPIR(fp, TRUE, p->seqs, 2, &seqinfo, &punct, &error))
   {
      if(error)
      {
         FREELIST(start, SEQDATA);
         return(NULL);
      }
      
      prev = p;
      ALLOCNEXT(p,SEQDATA);
      if(p==NULL)
      {
         FREELIST(start, SEQDATA);
         return(NULL);
      }
   }

   /* Free up last one which was unused                                 */
   free(p);
   prev->next = NULL;

   return(start);
}

   
/************************************************************************/
/*>char **ListToTable(SEQDATA *SeqList, int *nseq, int *seqlen)
   ------------------------------------------------------------
   Converts the linked list of sequences to a 2D table

   11.09.96 Original   By: ACRM
*/
char **ListToTable(SEQDATA *SeqList, int *nseq, int *seqlen)
{
   SEQDATA *p;
   char **table;
   int  i, j, len;

   *seqlen=strlen(SeqList->seqs[0]);

   /* First count the number of sequences and find the longest seq      */
   for(p=SeqList, *nseq=0; p!=NULL; NEXT(p))
   {
      (*nseq)++;
      if(strlen(p->seqs[0]) > *seqlen)
         *seqlen = strlen(SeqList->seqs[0]);
   }

   /* Allocate memory                                                   */
   if((table = (char **)blArray2D(sizeof(char), *nseq, (*seqlen + 1)))
      ==NULL)
      return(NULL);
   
   for(p=SeqList, i=0; p!=NULL; NEXT(p), i++)
   {
      len = strlen(p->seqs[0]);
      for(j=0; j<len; j++)
         table[i][j] = p->seqs[0][j];
      table[i][j] = '\0';
   }

   PadSeqs(table, *nseq, *seqlen);
   
   return(table);
}


/************************************************************************/
/*>void PadSeqs(char **SeqTable, int nseq, int seqlen)
   ---------------------------------------------------
   Pads all sequences in a table to the same length

   11.09.96 Original   By: ACRM
*/
void PadSeqs(char **SeqTable, int nseq, int seqlen)
{
   int i;
   for(i=0; i<nseq; i++)
   {
      if(strlen(SeqTable[i]) < seqlen)
         blPadchar(SeqTable[i], seqlen, '-');
   }
}


/************************************************************************/
/*>void DisplayScores(FILE *fp, char **SeqTable, int nseq, int seqlen,
                      int MaxInMatrix, int Method, BOOL Extended,
                      BOOL ignoreGaps)
   -------------------------------------------------------------------
   Display the variability scores for each position in the alignment

   11.09.96 Original   By: ACRM
   17.09.96 Added MaxInMatrix and prints amino acid list
   15.07.08 Added Extended parameter and printing
   13.03.25 Added ignoreGaps
*/
void DisplayScores(FILE *fp, char **SeqTable, int nseq, int seqlen, 
                   int MaxInMatrix, int Method, BOOL Extended,
                   BOOL ignoreGaps)
{
   int i, j;

   for(i=0; i<seqlen; i++)
   {
      if(Extended)
      {
         fprintf(fp,"%4d %9.6f ",
                 i+1,
                 CalcScore(SeqTable, nseq, seqlen, i, 
                           MaxInMatrix, Method, ignoreGaps));
      }
      else
      {
         fprintf(fp,"%4d %6.3f ",
                 i+1,
                 CalcScore(SeqTable, nseq, seqlen, i, 
                           MaxInMatrix, Method, ignoreGaps));
      }
      
      for(j=0; j<nseq; j++)
      {
         fprintf(fp,"%c",SeqTable[j][i]);
      }
      fprintf(fp,"\n");
   }
}


/************************************************************************/
/*>REAL CalcScore(char **SeqTable, int nseq, int seql, int pos, 
                  int MaxInMatrix, int Method, BOOL ignoreGaps)
   ------------------------------------------------------------
   Calculate the score for a given position in the alignment

   11.09.96 Original   By: ACRM
   17.09.96 Changed score to LONG rather than ULONG since return value
            from CalcMDMScore() can be -ve!
            Added MaxInMatrix
   18.09.96 Changed calculation of combined score
   11.08.15 Initialize e
   24.08.15 Add seql parameter and valdar01 method.  By: TCN
   13.03.25 Added ignoreGaps.   By: ACRM
*/
REAL CalcScore(char **SeqTable, int nseq, int seql, int pos, 
               int MaxInMatrix, int Method, BOOL ignoreGaps)
{
   REAL e = 0.0,
        e9, e21;
   
   /* Defines group membership for the amino acid types                 */
   static AMINOACID AA21Groups[] =
   {  { 'A', 1, { 0,  0}},
      { 'C', 1, { 1,  0}},
      { 'D', 1, { 2,  0}},
      { 'E', 1, { 3,  0}},
      { 'F', 1, { 4,  0}},
      { 'G', 1, { 5,  0}},
      { 'H', 1, { 6,  0}},
      { 'I', 1, { 7,  0}},
      { 'K', 1, { 8,  0}},
      { 'L', 1, { 9,  0}},
      { 'M', 1, {10,  0}},
      { 'N', 1, {11,  0}},
      { 'P', 1, {12,  0}},
      { 'Q', 1, {13,  0}},
      { 'R', 1, {14,  0}},
      { 'S', 1, {15,  0}},
      { 'T', 1, {16,  0}},
      { 'V', 1, {17,  0}},
      { 'W', 1, {18,  0}},
      { 'Y', 1, {19,  0}},
      { 'B', 2, { 2, 11}},
      { 'Z', 2, { 3, 13}},
      { '-', 1, {20,  0}},
      { 'X', 1, {20,  0}},
      { ' ', 0, { 0,  0}}
   };

   static AMINOACID AA9Groups[] =
   {  { 'A', 1, { 5,  0}},
      { 'C', 1, { 7,  0}},
      { 'D', 1, { 3,  0}},
      { 'E', 1, { 3,  0}},
      { 'F', 1, { 1,  0}},
      { 'G', 1, { 5,  0}},
      { 'H', 1, { 1,  0}},
      { 'I', 1, { 0,  0}},
      { 'K', 1, { 2,  0}},
      { 'L', 1, { 0,  0}},
      { 'M', 1, { 7,  0}},
      { 'N', 1, { 4,  0}},
      { 'P', 1, { 6,  0}},
      { 'Q', 1, { 4,  0}},
      { 'R', 1, { 2,  0}},
      { 'S', 1, { 4,  0}},
      { 'T', 1, { 4,  0}},
      { 'V', 1, { 0,  0}},
      { 'W', 1, { 1,  0}},
      { 'Y', 1, { 1,  0}},
      { 'B', 2, { 3,  4}},
      { 'Z', 2, { 3,  4}},
      { '-', 1, { 8,  0}},
      { 'X', 1, { 8,  0}},
      { ' ', 0, { 0,  0}}
   };

   switch(Method)
   {
   case METH_MDM:
      return(MDMBasedScore(SeqTable, nseq, pos, MaxInMatrix, ignoreGaps));
   case METH_ENTROPY20:
      return((REAL)1.0 -
             EntropyScore(SeqTable, nseq, pos, AA21Groups, 21));
   case METH_ENTROPY8:
      return((REAL)1.0 -
             EntropyScore(SeqTable, nseq, pos, AA9Groups, 9)); 
   case METH_ENTROPY:
/*
      e21 = (REAL)1.0 - EntropyScore(SeqTable, nseq, pos, AA21Groups, 21);
      e9  = (REAL)1.0 - EntropyScore(SeqTable, nseq, pos, AA9Groups,  9);
      e   = (REAL)sqrt((double)((e21 * (e9 -(REAL)1.0 + 
                                        (REAL)20.0/(REAL)8.0)) / 
                                ((REAL)20.0/(REAL)8.0)));
      return(e);
*/
      e21 = EntropyScore(SeqTable, nseq, pos, AA21Groups, 21);
      e9  = EntropyScore(SeqTable, nseq, pos, AA9Groups,  9);
      e   = e21 * ((1.0 - (8.0/20.0))*e9 + (8.0/20.0));
      e   = 1.0 - e;
      return((REAL)e);
   case METH_VALDAR:
      return(valdarScore(SeqTable, pos, nseq, seql, MaxInMatrix));
   default:
      return((REAL)0.0);
   }

   return((REAL)0.0);
}


/************************************************************************/
/*>REAL MDMBasedScore(char **SeqTable, int nseq, int pos, int MaxInMatrix,
                      BOOL ignoreGaps)
   -----------------------------------------------------------------------
   Calculate the score for a given position in the alignment using the
   MDM Method

   11.09.96 Original   By: ACRM
   17.09.96 Changed score to LONG rather than ULONG since return value
            from CalcMDMScore() can be -ve!
            Added MaxInMatrix
   13.09.22 Changed to LONG for count as well as score - can change 
            these to VLONG if needed
   13.03.25 Added ignoreGaps
*/
REAL MDMBasedScore(char **SeqTable, int nseq, int pos, int MaxInMatrix,
                   BOOL ignoreGaps)
{
   LONG  i, j;
   VLONG count,
         score;
   char  res1,
         res2;
   
   count = 0;
   score = 0L;
   for(i=0; i<nseq-1; i++)
   {
      res1 = SeqTable[i][pos];
      if(res1 == ' ')
         res1 = '-';

      if(!ignoreGaps || (res1 != '-'))
      {
         for(j=i+1; j<nseq; j++)
         {
            res2 = SeqTable[j][pos];
            if(res2 == ' ')
               res2 = '-';

            if(!ignoreGaps || (res2 != '-'))
            {
               if(++count > MAXCHECK)
               { 
                  fprintf(stderr, "Score count too large to store!\n");
                  exit(1);
               }
               
               score += blCalcMDMScore(res1, res2);
               if(score > MAXCHECK)
               { 
                  fprintf(stderr, "Score (%Ld) too large to store!\n",
                          score);
                  exit(1);
               }
            }
         }
      }
   }

#ifdef DEBUG
   printf("Score: %Ld Count: %Ld\n", score, count);
#endif
   
   return(((REAL)score/(REAL)count)/(REAL)MaxInMatrix);
}

/************************************************************************/
/*>REAL EntropyScore(char **SeqTable, int nseq, int pos, 
                     AMINOACID *aminoacids, int NGroups)
   -----------------------------------------------------
   Calculates an entropy score based on the equation:
   S = - \sum_i p_i \log p_i
   where p_i is n_i/N

   A table of classifications for the amino acids is supplied to the 
   routine such that the amino acid types may be grouped. This table
   also allows residues to belong to 2 groups to account for B (ASX)
   and Z (GLX).

   The result is scaled such that the entropy runs from 0 (all conserved)
   to 1.0 (maximum variability)

   17.09.96 Original   By: ACRM
*/
REAL EntropyScore(char **SeqTable, int nseq, int pos, 
                  AMINOACID *aminoacids, int NGroups)
{
   REAL entropy = (REAL)0.0,
        *count;
   int  type, i, j;

   /* Allocate memory to store the counts and zero them                 */
   if((count = (REAL *)malloc(NGroups * sizeof(REAL)))==NULL)
      return((REAL)9999.0);
   for(i=0; i<NGroups; i++)
      count[i] = (REAL)0.0;
   
   /* For each recognised amino acid type                               */
   for(type=0; aminoacids[type].NGroup != 0; type++)
   {
      /* Count through the sequences to see how many amino acids are of 
         this type
      */
      for(i=0; i<nseq; i++)
      {
         /* We've got an amino acid of this type                        */
         if(SeqTable[i][pos] == aminoacids[type].res)
         {
            /* We allow amino acids to belong to more than one group to
               handle B (ASX) and Z (GLX).

               For each group to which this residue belongs             
            */
            for(j=0; j<aminoacids[type].NGroup; j++)
            {
               /* Increment the count by 1 over the number of groups to
                  which this residue belongs
               */
               count[aminoacids[type].group[j]] += 
                  (REAL)1.0/(REAL)aminoacids[type].NGroup;
            }
         }
      }
   }

   /* Now run through all the counts and convert them to fractions      */
   for(i=0; i<NGroups; i++)
      count[i] /= nseq;

   /* Add up the entropy score                                          */
   for(i=0; i<NGroups; i++)
   {
      if(count[i] > (REAL)0.0)
      {
         entropy -= count[i] * 
            (REAL)log((double)count[i]);
      }
   }

   /* At this stage entropy runs from 0 (absolutely conserved) up.
      We now divide by log of the number of sequences or number of
      groups (whichever is smaller) such that the maximum value is 1.0
   */
   entropy /= log(MIN(nseq,NGroups));
   
   return(entropy);
}

/************************************************************************/
/*>REAL valdarScore(char **SeqTable, int pos, int numSeqs, int seqlen,
                    int MaxInMatrix)
   -----------------------------------------------------------------------
   Calculate the conservation score of an alignment position, using the
   valdar01 method.
   
   Returns 9999.0 on error

   20.08.15 Original   By: TCN
*/
REAL valdarScore(char **SeqTable, int pos, int numSeqs, int seqlen,
                 int MaxInMatrix) 
{
   int  i, j;
   LREAL matrixScore;
   LREAL weightedSum = 0;

   if(sSeqWeights == SEQWEIGHTS_UNITIALIZED)
   {
      if(initSequenceWeights(SeqTable, numSeqs, seqlen, 
                             MaxInMatrix)==FALSE)
      {
         return((REAL)9999.0);
      }
   }
   
   if(sLambda == LAMBDA_UNITIALIZED)
   {
      initLambda(numSeqs);
   }
   
   for(i=0; i<numSeqs; i++)
   {
      for(j=i+1; j<numSeqs; j++)
      {
         matrixScore = valdarMatrixScore(SeqTable[i][pos], 
                                         SeqTable[j][pos], 
                                         MaxInMatrix);
         weightedSum += sSeqWeights[i] * sSeqWeights[j] * matrixScore;
      }
   }

   return(sLambda * weightedSum);
}

/************************************************************************/
/*>void initLambda(int numSeqs)
   ----------------------------
   Initializes the global static REAL lambda, a scalar used for
   calculating conscores using the valdar01 method.
   
   20.08.2015 Original   By: TCN
*/
void initLambda(int numSeqs) 
{
   REAL weightSum = 0;
   int  i, j;
    
   for(i=0; i<numSeqs; i++)
   {
      for(j=i+1; j<numSeqs; j++)
      {
         weightSum += sSeqWeights[i] * sSeqWeights[j];
      }
   }

   sLambda =  ((REAL)1 / weightSum); 
}

/************************************************************************/
/*>REAL initSequenceWeights(char **SeqTable, int numSeqs, int seqlen,
                            int MaxInMatrix)
   -----------------------------------------------------------------------
   Initializes the global static array sSeqWeights, used for the
   valdar01 method. Each sequence is given a weight according to its
   evolutionary distance from the other in the alignment.
   
   20.08.15 Original   By: TCN
*/
BOOL initSequenceWeights(char **SeqTable, int numSeqs, int seqlen,
                         int MaxInMatrix) 
{
   int  i, j;
   REAL seqDistSum;
   REAL **seqDistTable = NULL;
   
   if((seqDistTable = getSeqDistTable(SeqTable, numSeqs, seqlen, 
                                      MaxInMatrix))==NULL)
      return(FALSE);
   
   if((sSeqWeights = malloc(sizeof(REAL) * numSeqs))==NULL)
      return(FALSE);
   
   for(i=0; i<numSeqs; i++)
   {
      seqDistSum = (REAL)0.0;
      for(j=0; j<numSeqs; j++)
      {
         if(i == j)
            continue;
         seqDistSum += seqDistTable[i][j];
      }
      sSeqWeights[i] = seqDistSum / ((REAL)numSeqs - (REAL)1.0);
   }
   
   return(TRUE);
}

/************************************************************************/
/*>REAL **getSeqDistTable(char **SeqTable, int numSeqs, int seqlen, 
                          int MaxInMatrix)
   ----------------------------------------------------------------
   Get 2d table of inter-sequence evolutionary distances.
   
   20.08.15 Original   By: TCN
   10.08.22 Changed to use blArray2D() and fixed bug in second dimension
            (which was seqlen instead of numSeqs)  By: ACRM
*/
REAL **getSeqDistTable(char **SeqTable, int numSeqs, int seqlen, 
                       int MaxInMatrix)
{
   int i, j;
   REAL **seqDistTable;

   if((seqDistTable = (REAL **)blArray2D(sizeof(REAL),
                                         numSeqs, numSeqs))==NULL)
      return(NULL);
    
   for(i=0; i<numSeqs; i++)
   {
      for(j=0; j<numSeqs; j++)
      {
         if(i == j)
            continue;
         seqDistTable[i][j] = getInterSeqDistance(SeqTable, i, j, 
                                                  seqlen, MaxInMatrix);  
      }
   }

   return seqDistTable;
}


/************************************************************************/
/*>REAL getInterSeqDistance(char **SeqTable, int seqAindex, int seqBindex,
                            int seqlen, int MaxInMatrix)
   -----------------------------------------------------------------------
   Calculates the evolutionary distance between two sequences.
   
   20.08.2015 Original   By: TCN
*/
REAL getInterSeqDistance(char **SeqTable, int seqAindex, int seqBindex,
                         int seqlen, int MaxInMatrix)
{
   int  nonGapPosCount,
        matrixScoreSum = 0,
        pos;
   char res1, res2;

   nonGapPosCount = getNonGapPosCount(SeqTable, seqAindex, 
                                      seqBindex, seqlen);
    
   for(pos=0; pos<seqlen; pos++)
   {
      res1 = SeqTable[seqAindex][pos];
      res2 = SeqTable[seqBindex][pos];
        
      if (res1 == ' ')
         res1 = '-';
      if (res2 == ' ')
         res2 = '-';
       
      if(! (res1 == '-' && res2 == '-'))
      {
         matrixScoreSum += valdarMatrixScore(res1, res2, MaxInMatrix);
      }
   }

   return ((REAL)1.0 - ((REAL)matrixScoreSum / (REAL)nonGapPosCount));
}

/************************************************************************/
/*>int getNonGapPosCount(char **SeqTable, int seqAindex, int seqBindex, 
                         int seqlen) 
   --------------------------------------------------------------------
   Given two sequences, counts the number of positions where at least 
   one sequence is non-gap.
   
   20.08.2015 Original   By: TCN
*/
int getNonGapPosCount(char **SeqTable, int seqAindex, int seqBindex, 
                      int seqlen)
{
   int  pos,
        nonGapPosCount = 0;
   char res1, res2;
    
   for(pos=0; pos<seqlen; pos++)
   {
      res1 = SeqTable[seqAindex][pos];
      res2 = SeqTable[seqBindex][pos];
        
      if (res1 == ' ')
         res1 = '-';
      if (res2 == ' ')
         res2 = '-';
       
      if(! (res1 == '-' && res2 == '-'))
      {
         nonGapPosCount += 1;
      }
   }

   return nonGapPosCount;
}


/************************************************************************/
/*>REAL valdarMatrixScore(char res1, char res2, int MaxInMatrix)
   -------------------------------------------------------------
   Calculate the score for a given position in the alignment for valdar01
   method.
   
   20.08.2015 Original   By: TCN
*/
REAL valdarMatrixScore(char res1, char res2, int MaxInMatrix)
{
    LONG score = 0L;

    if (res1 == '-' || res2 == '-')
    {
       return ((REAL)0.0);
    }
    else
    {
        score = blCalcMDMScore(res1, res2);
        return ((REAL)score / (REAL)(MaxInMatrix));
    }

    return ((REAL)0.0);
}


/************************************************************************/
/*>BOOL ReadAndScoreSingle(char *single, FILE *out, int MaxInMatrix, 
                           int Method, BOOL Extended, BOOL doLog,
                           REAL maxFraction, BOOL reduceData)
   -----------------------------------------------------------------------
   Routine which takes the residue distribution for a single position,
   and calculates and displays the variability scores.

   10.08.22 Original   By: ACRM
   04.10.22 Added doLog, maxFraction and reduceData
*/
BOOL ReadAndScoreSingle(char *single, FILE *out, int MaxInMatrix,
                        int Method, BOOL Extended, BOOL doLog,
                        REAL maxFraction, BOOL reduceData)
{
   char    **SeqTable;
   int     seqlen, nseq;

   seqlen = 1;

   if((SeqTable = ParseSingle(single, &nseq, doLog, maxFraction,
                              reduceData))==0)
      return(FALSE);
   
   /* Calculate and print scores                                        */
   DisplayScores(out, SeqTable, nseq, seqlen, MaxInMatrix, Method,
                 Extended, FALSE);

   /* Free memory from sequence table                                   */
   blFreeArray2D((char **)SeqTable, nseq, seqlen);

   return(TRUE);
}

/************************************************************************/
/*>char **ParseSingle(char *single, int *nSeq, BOOL doLog, 
                      REAL maxFraction, BOOL reduceData)
   -------------------------------------------------------
   Parses the input single column data (in the form "A:m,C:n,D:o,...")
   and creates a sequence table which is 1 column wide and m+n+o+...
   rows long (i.e. in the normal format for full sequences).

   10.08.22  Original   By: ACRM
   04.10.22  Added doLog code; 
             Added maxFraction code
             Added reduceData code
*/
char **ParseSingle(char *single, int *nSeq, BOOL doLog, REAL maxFraction,
                   BOOL reduceData)
{
   char **fields   = NULL,
        **seqTable = NULL;
   int  nFields,
        numAA,
        row,
        firstAA = 0,
        otherAAs;
   BOOL error = FALSE;
   REAL reductionFactor = 0.0;
   

   /* Split the 'single' string on commas                               */
   if((fields = blSplitStringOnCommas(single, MINSINLEN))==NULL)
      return(NULL);
   
   /* Count the number of amino acids                                   */
   numAA=0;
   for(nFields=0; fields[nFields][0] != '\0'; nFields++)
   {
      char *ptr;
      int  nAA;
      
      /* Look for the colon                                             */
      if((ptr=strchr(fields[nFields],':'))==NULL)
      {
         error = TRUE;
         break;
      }
      /* Find the number after the colon                                */
      if(!sscanf(ptr+1, "%d", &nAA))
      {
         error = TRUE;
         break;
      }
      if(nAA <= 0)
      {
         fprintf(stderr, "Error: counts must be >0 (%s)\n",
                 fields[nFields]);
         exit(1);
      }
      
      if(doLog)
      {
         REAL logAA = 1 + LOGSCALE*log(nAA);
         nAA = (int)logAA;
      }
      if(nFields == 0)
      {
         firstAA = nAA;
      }
         
      /* Add to the amino acid count                                    */
      numAA += nAA;
   }
   
   if(reduceData && (numAA > MAXDATA))
   {
      reductionFactor = (REAL)MAXDATA / (REAL)numAA;
   }
   else
   {
      reduceData = FALSE;
   }
   

   if(error)
      return(NULL);

   if(maxFraction > TINY)
   {
      numAA   -= firstAA;
      otherAAs = numAA;
      numAA    = (int)(0.5+((REAL)numAA / (1-maxFraction)));
      firstAA  = (int)(0.5+((maxFraction * otherAAs) / (1-maxFraction)));
   }

   /* Allocate the output 2D array                                      */
   if((seqTable = (char **)blArray2D(sizeof(char), numAA, 1))==NULL)
      return(NULL);
   
   /* Populate the output 2D array. We don't need the error checking
      now since we have already checked the data
   */
   row   = 0;
   numAA = 0;
   for(nFields=0; fields[nFields][0] != '\0'; nFields++)
   {
      char *ptr;
      int  nAA, i;
      /* Look for the colon                                             */
      ptr=strchr(fields[nFields],':');
      *ptr = '\0';
      
      /* Find the number after the colon                                */
      sscanf(ptr+1, "%d", &nAA);

      if(doLog)
      {
         REAL logAA = 1 + LOGSCALE*log(nAA);
         nAA = (int)logAA;
      }
      if((nFields == 0) && (maxFraction > TINY))
      {
         nAA = firstAA;
      }
      if(reduceData)
      {
         nAA = (int)(0.5 + (reductionFactor * nAA));
      }
      numAA += nAA;
      
      /* Populate                                                       */
      for(i=0; i<nAA; i++)
      {
         seqTable[row][0] = fields[nFields][0];
         row++;
      }
   }

   /* Free the field storage                                            */
   blFreeArray2D((char **)fields, nFields, MINSINLEN);

   /* Output and return                                                 */
   *nSeq = numAA;
   return(seqTable);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   17.09.96 Original   By: ACRM
   15.07.08 V1.3 - added -x
   11.08.15 V1.4
   24.08.15 V1.5 (added -d Valdar method)
   11.10.19 V1.6 Fixed reading of matrix with -m
   10.08.22 V1.7 Added -s
   13.09.22 V1.7.1
   04.10.22 V1.8
   22.11.22 V1.8.1
   13.03.25 V1.9
*/
void Usage(void)
{
   fprintf(stderr,"\nScoreCons V1.9 (c) 1996-2025 Prof. Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"          valdar01 scoring implemented by Tom \
Northey\n");

   fprintf(stderr,"\nUsage: scorecons [-m matrixfile] [-a|-g|-e|-d] \
[-x] [-i] [alignment.pir [output.dat]]\n");
   fprintf(stderr," -or-  scorecons -s A:n,C:n,D:n,... [-m matrixfile] \
[-a|-g|-e|-d] [-i] [-r|-l|-f[=n]]\n");
   fprintf(stderr,"                 [-x] [output.dat]\n");
   fprintf(stderr,"       -m Specify the mutation matrix (Default: %s)\n",
           MUTMAT);
   fprintf(stderr,"       -a Score by entropy method per residue\n");
   fprintf(stderr,"       -g Score by entropy method, 8 groups of \
residues\n");
   fprintf(stderr,"       -e Score by combined entropy method\n");
   fprintf(stderr,"       -d Score by the valdar01 method\n");
   fprintf(stderr,"       -x Extended precision output\n");
   fprintf(stderr,"       -i Ignore gaps\n");
   
   fprintf(stderr,"       -s Score a single column of an alignment \
specifying residue counts\n");
   fprintf(stderr,"          on the command line\n");
   fprintf(stderr,"       -r Reduce dataset sizes for speed \
(used with -s) - maximum\n");
   fprintf(stderr,"          datapoints %d\n", MAXDATA);
   fprintf(stderr,"       -l Scale counts by taking logs \
(used with -s)\n");
   fprintf(stderr,"       -f Set the count for the first AA to the \
specified fraction of the total\n");
   fprintf(stderr,"          <1.0, (used with -s) [Default: 0.5]\n");

   fprintf(stderr,"\nCalculates a conservation score between 0 and 1 \
for a PIR format\n");
   fprintf(stderr,"sequence alignment file. Output consists of the \
alignment position,\n");
   fprintf(stderr,"the score and the residues seen at that \
position.\n");

   fprintf(stderr,"\nBy default, the conservation score is calculated \
from an updated version\n");
   fprintf(stderr,"of the Dayhoff mutation matrix. Alternatively, a \
statistical entropy\n");
   fprintf(stderr,"scoring method or the valdar01 method may be \
employed.\n");

   fprintf(stderr,"\nThe grouped entropy method places amino acids into \
groups: \n");
   fprintf(stderr,"ILV, FHWY, KR, DE, NQST, AG, P, CM. Deletions and X \
residues form another\n");
   fprintf(stderr,"group while B(ASX) and Z(GLX) are placed in both the \
DE and NQST groups.\n");

   fprintf(stderr,"\nThe combined entropy method uses a combination of \
the grouped and\n");
   fprintf(stderr,"ungrouped entropy scores. The values are multiplied \
together, but\n");
   fprintf(stderr,"weighted such that the the grouped score contributes \
less owing to the\n");
   fprintf(stderr,"loss of information content.\n");

   fprintf(stderr,"\nThe valdar01 method is a scoring method developed \
by Will Valdar.\n");
   fprintf(stderr,"It uses a weighting scheme based on the evolutionary \
distance\n");
   fprintf(stderr,"between aligned sequences.\n");

   fprintf(stderr,"\nNote that to use -r, -l or -f, you must be using \
-s and that only one\n");
   fprintf(stderr,"of -r, -l and -f may be used\n");
   
   fprintf(stderr,"\nThis program should not be confused by the program \
of the same name by\n");
   fprintf(stderr,"Will Valdar. This program was written first and \
Will's program was\n");
   fprintf(stderr,"inspired by this one.\n\n");
}
