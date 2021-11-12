/************************************************************************/
/**
   \file       pdbrepair.c
   
   \version    
   \date       29.10.21
   \brief      
   
   \copyright  (c) Prof. Andrew C. R. Martin 2021
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0    Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/general.h"
#include "bioplib/array.h"
#include "bioplibnew.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF   160
#define MAXCHAINS 240
#define GAPPEN    2
#define safetoupper(x) ((islower(x))?toupper(x):(x))
#define safetolower(x) ((isupper(x))?tolower(x):(x))
#define CONECT_TOL 0.2

/************************************************************************/
/* Prototypes
*/
static char *sCombineSequence(char *align1, char *align2, int align_len,
                              BOOL upper);
char *strdup(const char *s);

#ifdef UNUSED_JUNK
static int sGetPDBChains(PDB *pdb, char *chains);
/************************************************************************/
/*>static int sGetPDBChains(PDB *pdb, char *chains)
   ----------------------------------------
*//**

   Extracts a list of chains from a PDB linked list

-  22.08.97 Original   By: ACRM
-  10.03.15 Checks number of chains
*/
static int sGetPDBChains(PDB *pdb, char *chains)
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
#endif


/************************************************************************/
/*>char *blFixSequence(char *seqresSequence, char *atomSequence, 
                       char **seqresChains, char **atomChains,
                       char **outchains, BOOL IgnoreSEQRES, int nAtomChains)
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
-  27.07.21 Returns the atomSequence if SEQRES not read
*/
char *blFixSequence(char *seqresSequence, char *atomSequence,
                    char **seqresChains, char **atomChains,
                    char **outchains, BOOL IgnoreSEQRES, int nAtomChains,
                    BOOL upper, BOOL quiet, char *label)
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


   if((seqresSequence == NULL) || (atomSequence == NULL))
   {
      for(i=0; i<nAtomChains; i++)
      {
         strcpy(outchains[i], atomChains[i]);
      }
      return(strdup(atomSequence));
   }
   
   
   /* Set flags to say we haven't handled the sequences yet             */
   for(i=0; i<MAXCHAINS; i++)
   {
      DoneSEQRES[i] = FALSE;
      DoneATOM[i]   = FALSE;
   }
   
   /* If the sequences and chains are identical just copy one of them
      and return
   */
   if(!strcmp(seqresSequence,atomSequence) &&
      !strcmp(seqresChains[0],atomChains[0])) /* FIXME! */
   {
      for(i=0; i<nAtomChains; i++)
      {
         strcpy(outchains[i], seqresChains[i]);
      }
      return(strdup(atomSequence));
   }

   /* Create a temporary buffer to store a sequence                     */
   len1 = strlen(seqresSequence);
   len2 = strlen(atomSequence);
   if((buffer = (char *)malloc((1+MAX(len1, len2))
                               * sizeof(char)))==NULL)
   {
      return(NULL);
   }

   /* See how many chains there are and create arrays of char pointers
      to store the separate chains
   */
   nchain[0] = blCountchar(seqresSequence,   '*');
   if(len1 && seqresSequence[len1-1] != '*')
      nchain[0]++;
   nchain[1] = blCountchar(atomSequence, '*');
   if(len2 && atomSequence[len2-1] != '*')
      nchain[1]++;

   for(i=0; i<2; i++)
   {
      if((seqs[i] = (char **)malloc(nchain[i] * sizeof(char *)))==NULL)
         return(NULL);
   }

   /* Transfer the individual chains into the split arrays              */
   for(i=0; i<2; i++)
   {
      strcpy(buffer,((i==0)?seqresSequence:atomSequence));
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
         if(CHAINMATCH(seqresChains[i], atomChains[j]))
         {
            DoneSEQRES[i] = TRUE;
            DoneATOM[j]   = TRUE;
            strcpy(outchains[NOutChain++], seqresChains[i]);
            
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

               if((combseq=sCombineSequence(align1,align2,align_len, upper))
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
         strcpy(outchains[NOutChain++], atomChains[i]);
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
            if(!upper)
            {
               LOWER(seqs[0][i]);
            }
            
            strcat(outseq,seqs[0][i]);
            strcat(outseq,"*");
            strcpy(outchains[NOutChain++], seqresChains[i]);
            
            if(!quiet)
            {
               if(label != NULL)
               {
                  fprintf(stderr,"Warning: Chain %s from SEQRES records \
not found in ATOM records%s%s\n",
                          seqresChains[i],
                          ((label[0])?" Label: ":""),
                          ((label[0])?label:""));
               }
               else
               {
                  fprintf(stderr,"Warning: Chain %s from SEQRES records \
not found in ATOM records\n",
                          seqresChains[i]);
               }
            }
         }
      }
   }

#ifdef DEBUG
   fprintf(stderr, "\n=============================================\n\n");
#endif

   /*** NEEDS TO FREE seqs[][] too                                    ***/
   free(buffer);         /*  11.06.15                                   */
   return(outseq);
}

/************************************************************************/
/*>static char *sCombineSequence(char *align1, char *align2, int align_len)
   ----------------------------------------------------------------
*//**

   Combine the information from the two sequences

-  22.08.97 Original   By: ACRM
*/
static char *sCombineSequence(char *align1, char *align2, int align_len,
                              BOOL upper)
{
   static char *outseq = NULL;
   int         i;

   if((outseq=(char *)malloc((align_len+1) + sizeof(char)))==NULL)
      return(NULL);

#ifdef DEBUG
   fprintf(stderr,"\n----------------------------------------------\n");
   fprintf(stderr,"Alignment:\n");
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
#ifdef DEBUG
         fprintf(stderr, "AL1: %c AL2: %c\n", align1[i], align2[i]);
#endif
         if(align2[i] == '-')
         {
            if(upper)
               outseq[i] = safetoupper(align1[i]);
            else
               outseq[i] = safetolower(align1[i]);
         }
         else
         {
            outseq[i] = safetoupper(align2[i]);
         }
      }
   }
   outseq[align_len] = '\0';
#ifdef DEBUG
   fprintf(stderr, "\n(Sequence now assembled)\n");
#endif
   return(outseq);
}


/************************************************************************/
void blRenumResiduesPDB(PDB *pdb, int offset)
{
   PDB *p,
       *prev      = pdb;
   int resnum     = offset,
       prevResnum = pdb->resnum;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!CHAINMATCH(p->chain, prev->chain))
      {
         resnum = offset;
      }
      else if((p->resnum != prevResnum) ||
              !INSERTMATCH(p->insert, prev->insert)
             )
      {
         resnum++;
      }

      prevResnum = p->resnum;
      prev       = p;
      p->resnum  = resnum;
   }
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->insert[0] = ' ';
      p->insert[1] = '\0';
   }
   
}


/************************************************************************/
STRINGLIST *blCreateSEQRES(PDB *pdb)
{
   HASHTABLE  *seqByChain = NULL;
   STRINGLIST *seqres     = NULL;
   int        nChains     = 0;
   char       **chains    = NULL,
              buffer[MAXBUFF],
              aa[8];
   
   if((seqByChain = blPDB2SeqXByChain(pdb))!=NULL)
   {
      if((chains = blGetPDBChainLabels(pdb, &nChains))==NULL)
      {
         return(NULL);
      }
      else
      {
         int i;
         
         for(i=0; i<nChains; i++)
         {
            int  chainLen,
                 lineNum   = 1,
                 resNum    = 0;
            char *sequence = NULL;
            BOOL Stored    = TRUE;

            sequence = blGetHashValueString(seqByChain, chains[i]);
            if(sequence != NULL)
            {
               chainLen = strlen(sequence);
               for(resNum=0; resNum<chainLen; resNum++)
               {
                  if(!(resNum%13))
                  {
                     if(!Stored)
                     {
                        strcat(buffer, "\n");
                        seqres = blStoreString(seqres, buffer);
                        Stored = TRUE;
                     }
                     
                     sprintf(buffer, "SEQRES%4d %c%5d  ",
                             lineNum++,
                             chains[i][0],
                             chainLen);
                  }
                  sprintf(aa, "%-4s", blOnethr(sequence[resNum]));
                  strcat(buffer, aa);
                  Stored = FALSE;
               }
               if(!Stored)
               {
                  strcat(buffer, "\n");
                  seqres = blStoreString(seqres, buffer);
                  Stored = TRUE;
               }
            }
         }
      }
   }
   
   if(seqByChain != NULL)
      blFreeHash(seqByChain);
   
   return(seqres);
}

/************************************************************************/
void blReplacePDBHeader(WHOLEPDB *wpdb, char *recordType,
                        STRINGLIST *replacement)
{
   STRINGLIST *s,
              *previousRecord = NULL,
              *firstRecord    = NULL,
              *prev           = NULL,
              *next           = NULL,
              *nextRecord     = NULL;
   BOOL       gotHeader       = FALSE;

   /* Find the records before and after the type we are looking for     */
   for(s=wpdb->header; s!=NULL; NEXT(s))
   {
      if(!strncmp(s->string, recordType, 6))
      {
         if(!gotHeader)
         {
            /* This is the first record of the header type we are
               looking for
            */
            firstRecord    = s;
            previousRecord = prev;
            gotHeader      = TRUE;
         }
      }
      else if(gotHeader)
      {
         /* This is the first record after the geader type we are
            looking for
         */
         nextRecord = s;
         break;
      }
      prev = s;
   }

   /* Return if we didn't find the relevant header                      */
   if(!gotHeader)
      return;

   /* Free the ones we don't need                                       */
   for(s=firstRecord; s!=nextRecord; s=next)
   {
      next = s->next;
      FREE(s->string);
      FREE(s);
   }

   /* Patch in the replacement                                          */
   if(replacement != NULL)
   {
      if(previousRecord == NULL)
      {
         wpdb->header = replacement;
      }
      else
      {
         previousRecord->next = replacement;
      }
      
      s=replacement;
      LAST(s);
      s->next = nextRecord;
   }
   else  /* No replacement: just link the record before to the one after*/
   {
      if(previousRecord == NULL)
      {
         wpdb->header = nextRecord;
      }
      else
      {
         previousRecord->next = nextRecord;
      }
   }
}

