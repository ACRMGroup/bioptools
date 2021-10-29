/************************************************************************/
/**

   \file       pdbrepair.c
   
   \version    V0.1
   \date       29.10.21
   \brief      Count residues and atoms in a PDB file
   
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
/* #include <ctype.h> */
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/general.h"
#include "bioplib/array.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF   160
#define MAXCHAINS 240
#define GAPPEN    2
#define safetoupper(x) ((islower(x))?toupper(x):(x))
#define safetolower(x) ((isupper(x))?tolower(x):(x))


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);
PDB *RepairPDB(PDB *pdb, char *fixedseq);
int GetPDBChains(PDB *pdb, char *chains);
char *FixSequence(char *seqres, char *sequence, char **seqchains, 
                  char **atomchains, char **outchains, BOOL IgnoreSEQRES,
                  int nAtomChains);
char *CombineSequence(char *align1, char *align2, int align_len);
void PrintNumbering(FILE *out, PDB *pdb, MODRES *modres);
char *strdup(const char *s);



/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for counting residues & atoms

-  29.10.21 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
        
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      FILE *in  = stdin,
           *out = stdout;
      
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         WHOLEPDB *wpdb;

         if(((wpdb=blReadWholePDB(in))==NULL) || (wpdb->pdb==NULL))
         {
            fprintf(stderr,"No atoms read from input file\n");
         }
         else
         {
            PDB  *fixedPDB,
                 *pdb;
            MODRES *modres = NULL;
            char *seqres = NULL,
               *sequence = NULL,
               *fixedsequence = NULL,
               **seqchains = NULL,
               **outchains = NULL,
               **atomchains = NULL;
            int nAtomChains, len1;
            
            pdb = wpdb->pdb;
            
            if((seqchains = (char **)blArray2D(sizeof(char), MAXCHAINS, blMAXCHAINLABEL))==NULL)
            {
               fprintf(stderr,"Error: No memory for seqchains array\n");
               return(1);
            }

            /* Read MODRES and SEQRES records                                 */
            modres = blGetModresWholePDB(wpdb);
            seqres = blGetSeqresAsStringWholePDB(wpdb, seqchains, modres, TRUE);

            /* Extract sequence from PDB linked list                             */
            if((atomchains = blGetPDBChainLabels(pdb, &nAtomChains)) == NULL)
            {
               fprintf(stderr,"Error: No memory for atom chain labels\n");
               return(1);
            }
            
            /* Convert PDB linked list to a sequence                             */
            if((sequence = blPDB2SeqX(pdb))==NULL)
            {
               fprintf(stderr,"Error: No memory for sequence data\n");
               return(1);
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

            /* Fiddle with sequences to combine information from SEQRES and ATOM
               records
            */
            if((fixedsequence = FixSequence(seqres,sequence,seqchains,
                                            atomchains,outchains,FALSE,
                                            nAtomChains))
               ==NULL)
               return(1);


            
            fixedPDB = RepairPDB(pdb, fixedsequence);
            blWritePDB(out, fixedPDB);
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \param[out]     *outfile    Output filename (or blank string)
   \return                     Success

   Parse the command line

-  29.10.21 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
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

-  29.10.21 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\npdbrepair V1.0 (c) 2021 Prof. Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"\nUsage: pdbrepair [in.pdb [out.pdb]]\n");
   
   fprintf(stderr,"If files are not specified, stdin and stdout are \
used.\n");
   fprintf(stderr,"Currently just adds missing ATOM records for atoms \
based on the\n");
   fprintf(stderr,"SEQRES records. Coordinates are set to 9999.999\n\n");
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
-  27.07.21 Returns the sequence if SEQRES not read
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


   if((seqres == NULL) || (sequence == NULL))
   {
      for(i=0; i<nAtomChains; i++)
      {
         strcpy(outchains[i], atomchains[i]);
      }
      return(strdup(sequence));
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

            strcat(outseq,seqs[0][i]);
            strcat(outseq,"*");
            strcpy(outchains[NOutChain++], seqchains[i]);
            
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
         outseq[i] = safetoupper(align1[i]);
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




PDB *RepairPDB(PDB *pdb, char *fixedseq)
{
   return(pdb);
}

