/************************************************************************/
/**


TODO:

-t option to trim sequences (I.e. ignore lower case letters at the ends of the alignments)

Maybe need a local sequence alignment and/or option to accept the ATOM residues rather than SEQRES

Code to rewrite SEQRES based on coordinates.

Implement the fix sequence code by calling routines from mutmodel

Need to fill in missing atoms

Check with other modified residues


   \file       pdbrepair.c
   
   \version    V0.1
   \date       29.10.21
   \brief      Add missing ATOM records based on SEQRES
   
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

typedef struct _restype
{
   char resnam[8],
        aa;
   int  hetatm;
   char atnams[MAXATINAA+1][8];
}  RESTYPE;
   

/************************************************************************/
/* Globals
*/

RESTYPE gResTypes[] =
{
   {"ALA",'A',0,{"N   ","CA  ","C   ","O   ","CB  ",""}},
   {"CYS",'C',0,{"N   ","CA  ","C   ","O   ","CB  ","SG  ",""}},
   {"ASP",'D',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","OD1 ","OD2 ",""}},
   {"GLU",'E',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","OE1 ","OE2 ",""}},
   {"PHE",'F',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ","CE1 ","CE2 ","CZ  ",""}},
   {"GLY",'G',0,{"N   ","CA  ","C   ","O   ",""}},
   {"HIS",'H',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","ND1 ","CD2 ","CE1 ","NE2 ",""}},
   {"ILE",'I',0,{"N   ","CA  ","C   ","O   ","CB  ","CG1 ","CG2 ","CD1 ",""}},
   {"LYS",'K',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","CE  ","NZ  ",""}},
   {"LEU",'L',0,{"N   ","CA  ","C   ","O   ","CB  ","CD1 ","CD2 ",""}},
   {"MET",'M',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","SD  ","CE  ",""}},
   {"ASN",'N',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","OD1 ","ND2 ",""}},
   {"PRO",'P',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ",""}},
   {"GLN",'Q',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","OE1 ","NE2 ",""}},
   {"ARG",'R',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","NE  ","CZ  ","NH1 ","NH2 ",""}},
   {"SER",'S',0,{"N   ","CA  ","C   ","O   ","CB  ","OG  ",""}},
   {"THR",'T',0,{"N   ","CA  ","C   ","O   ","CB  ","OG1 ","CG2 ",""}},
   {"VAL",'V',0,{"N   ","CA  ","C   ","O   ","CB  ","CG1 ","CG2 ",""}},
   {"TRP",'W',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ","NE1 ","CE2 ","CE3 ","CZ2 ","CZ3 ","CH2",""}},
   {"TYR",'Y',0,{"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ","CE1 ","CE2 ","CZ  ","OH  ",""}},
   {"PCA",'E',1,{"N   ","CA  ","CB  ","CG  ","CD  ","OE  ","C   ","O   ",""}},
   {""  ,'\0',0,{""}}
};

   

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);
PDB *RepairPDB(PDB *pdb, char *fixedSequence, MODRES *modres);
int GetPDBChains(PDB *pdb, char *chains);
char *FixSequence(char *seqresSequence, char *atomSequence, char **seqresChains, 
                  char **atomChains, char **outchains, BOOL IgnoreSEQRES,
                  int nAtomChains);
char *CombineSequence(char *align1, char *align2, int align_len);
void PrintNumbering(FILE *out, PDB *pdb, MODRES *modres);
void AppendThisResidue(PDB **ppdbOut, PDB **ppOut,
                       PDB *resIn,    PDB *nextResIn);
void AppendNewResidue(PDB **ppdbOut, PDB **ppOut,
                       char restype);
void AppendFixedResidue(PDB **ppdbOut, PDB **ppOut,
                        PDB *resIn,    PDB *nextResIn,
                        char restype);
void AppendRemainingAtoms(PDB **ppdbOut, PDB **ppOut,
                          PDB *pdbIn);



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
            char *seqresSequence = NULL,
               *atomSequence = NULL,
               *fixedSequence = NULL,
               **seqresChains = NULL,
               **outchains = NULL,
               **atomChains = NULL;
            int nAtomChains, len1;
            
            pdb = wpdb->pdb;
            
            if((outchains = (char **)blArray2D(sizeof(char),
                                               MAXCHAINS,
                                               blMAXCHAINLABEL))==NULL)
            {
               fprintf(stderr,"Error: No memory for outchains array\n");
               return(1);
            }
            if((seqresChains = (char **)blArray2D(sizeof(char),
                                               MAXCHAINS,
                                               blMAXCHAINLABEL))==NULL)
            {
               fprintf(stderr,"Error: No memory for seqresChains array\n");
               return(1);
            }

            /* Read MODRES and SEQRES records                           */
            modres = blGetModresWholePDB(wpdb);
            seqresSequence = blGetSeqresAsStringWholePDB(wpdb, seqresChains,
                                                 modres, TRUE);

            /* Get list of chains from the PDB linked list              */
            if((atomChains = blGetPDBChainLabels(pdb, &nAtomChains))
               == NULL)
            {
               fprintf(stderr,"Error: No memory for atom chain labels\n");
               return(1);
            }
            
            /* Convert PDB linked list to a sequence                    */
            if((atomSequence = blPDB2SeqX(pdb))==NULL)
            {
               fprintf(stderr,"Error: No memory for sequence data\n");
               return(1);
            }
            /* Append a * since blPDB2Seq() doesn't do this; note that
               this will have to change if we fix blPDB2Seq() in future
            */
            len1 = strlen(atomSequence);
            if((atomSequence=(char *)realloc((void *)atomSequence,
                                         (len1+2)*sizeof(char)))==NULL)
            {
               fprintf(stderr,"Error: No memory to expand sequence \
data\n");
               return(1);
            }
            strcat(atomSequence,"*");

#ifdef DEBUG
            fprintf(stderr,"\nSEQRES sequence:\n");
            fprintf(stderr,"%s\n", seqresSequence);
            fprintf(stderr,"\nATOM sequence:\n");
            fprintf(stderr,"%s\n", atomSequence);
#endif
            
            /* Fiddle with sequences to combine information from SEQRES
               and ATOM records
            */
            if((fixedSequence = FixSequence(seqresSequence,atomSequence,
                                            seqresChains,
                                            atomChains,outchains,FALSE,
                                            nAtomChains))
               ==NULL)
               return(1);
            
#ifdef DEBUG
            fprintf(stderr, "Fixed sequence:\n");
            fprintf(stderr, "%s\n", fixedSequence);
#endif
            fixedPDB = RepairPDB(pdb, fixedSequence, modres);
            wpdb->pdb = fixedPDB;
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
/*>char *FixSequence(char *seqresSequence, char *atomSequence, char **seqresChains, 
                     char **atomChains, char **outchains, 
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
-  27.07.21 Returns the atomSequence if SEQRES not read
*/
char *FixSequence(char *seqresSequence, char *atomSequence, char **seqresChains, 
                  char **atomChains, char **outchains, BOOL IgnoreSEQRES,
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
   if(!strcmp(seqresSequence,atomSequence) && !strcmp(seqresChains[0],atomChains[0])) /* FIXME! */
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
   nchain[0] = blCountchar(seqresSequence,   '*');        /* SEQRES sequence    */
   if(len1 && seqresSequence[len1-1] != '*')
      nchain[0]++;
   nchain[1] = blCountchar(atomSequence, '*');        /* ATOM sequence      */
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
            LOWER(seqs[0][i]);
            
            strcat(outseq,seqs[0][i]);
            strcat(outseq,"*");
            strcpy(outchains[NOutChain++], seqresChains[i]);
            
         }
      }
   }

#ifdef DEBUG
   fprintf(stderr, "\n=================================================\n\n");
#endif

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
   fprintf(stderr,"\n-----------------------------------------------\n");
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




void  AppendNewResidue(PDB **ppdbOut, PDB **ppOut,
                       char restype)
{
   PDB *pOut = *ppOut;
   int resOffset,
       atomOffset;
   char prevChain[8];

   for(resOffset=0; gResTypes[resOffset].aa != '\0'; resOffset++)
   {
      if(gResTypes[resOffset].aa == restype)
         break;
   }
   
   for(atomOffset=0; gResTypes[resOffset].atnams[atomOffset][0] != '\0'; atomOffset++)
   {
      if(*ppdbOut == NULL)
      {
         INIT((*ppdbOut), PDB);
         pOut = *ppdbOut;
         strcpy(prevChain, " ");
      }
      else
      {
         strcpy(prevChain, pOut->chain);
         ALLOCNEXT(pOut, PDB);
      }
         
      if(pOut == NULL)
      {
         FREELIST(*ppdbOut, PDB);
         *ppdbOut = NULL;
         *ppOut   = NULL;
         
         return;
      }

      /* Initialize this PDB item  */
      pOut->x              = 9999.999;
      pOut->y              = 9999.999;
      pOut->z              = 9999.999;
      pOut->occ            = 0.0;
      pOut->bval           = 99.0;
      pOut->access         = 0.0;
      pOut->radius         = 0.0;
      pOut->partial_charge = 0.0;

      pOut->atnum          = 0;
      pOut->resnum         = 0;
      pOut->formal_charge  = 0;
      pOut->nConect        = 0;
      pOut->entity_id      = 0;
      pOut->atomtype       = 0;
      
      pOut->element[0]     = gResTypes[resOffset].atnams[atomOffset][0];
      pOut->element[1]     = '\0';
      
      pOut->altpos         = ' ';
      pOut->secstr         = ' ';
      
      strcpy(pOut->record_type, gResTypes[resOffset].hetatm?"HETATM":"ATOM  ");
      strcpy(pOut->record_type, "INSERT");
      strcpy(pOut->atnam,       gResTypes[resOffset].atnams[atomOffset]);

      pOut->atnam_raw[0]   = ' ';
      strcpy(pOut->atnam_raw+1, gResTypes[resOffset].atnams[atomOffset]);
      pOut->atnam_raw[4]   = '\0';
      
      strcpy(pOut->resnam,      gResTypes[resOffset].resnam);
      strcpy(pOut->insert,      " ");
      strcpy(pOut->chain,       prevChain);
      strcpy(pOut->segid,       " ");

   }

   *ppOut = pOut;
}


void AppendThisResidue(PDB **ppdbOut, PDB **ppOut,
                       PDB *resIn,    PDB *nextResIn)
{
   PDB *pOut = *ppOut;
   PDB *pIn  = NULL;
   
   for(pIn = resIn; pIn != nextResIn; NEXT(pIn))
   {
      if(*ppdbOut == NULL)
      {
         INIT((*ppdbOut), PDB);
         pOut = *ppdbOut;
      }
      else
      {
         ALLOCNEXT(pOut, PDB);
      }
         
      if(pOut == NULL)
      {
         FREELIST(*ppdbOut, PDB);
         *ppdbOut = NULL;
         *ppOut   = NULL;
         
         return;
      }

      blCopyPDB(pOut, pIn);
   }
   *ppOut = pOut;
}


void AppendFixedResidue(PDB **ppdbOut, PDB **ppOut,
                        PDB *resIn,    PDB *nextResIn,
                        char restype)
{
   PDB *pOut = *ppOut;
   PDB *pIn  = NULL;
   
   for(pIn = resIn; pIn != nextResIn; NEXT(pIn))
   {
      if(*ppdbOut == NULL)
      {
         INIT((*ppdbOut), PDB);
         pOut = *ppdbOut;
      }
      else
      {
         ALLOCNEXT(pOut, PDB);
      }
         
      if(pOut == NULL)
      {
         FREELIST(*ppdbOut, PDB);
         *ppdbOut = NULL;
         *ppOut   = NULL;
         
         return;
      }

      blCopyPDB(pOut, pIn);
      strcpy(pOut->record_type, "FIXED ");
      
   }
   *ppOut = pOut;
}


void AppendRemainingAtoms(PDB **ppdbOut, PDB **ppOut,
                          PDB *pdbIn)
{
   PDB *pOut = *ppOut;
   PDB *pIn  = NULL;
   
   for(pIn = pdbIn; pIn != NULL; NEXT(pIn))
   {
      if(*ppdbOut == NULL)
      {
         INIT((*ppdbOut), PDB);
         pOut = *ppdbOut;
      }
      else
      {
         ALLOCNEXT(pOut, PDB);
      }
         
      if(pOut == NULL)
      {
         FREELIST(*ppdbOut, PDB);
         *ppdbOut = NULL;
         *ppOut   = NULL;
         
         return;
      }

      blCopyPDB(pOut, pIn);
   }
   *ppOut = pOut;
}


PDB *RepairPDB(PDB *pdbIn, char *fixedSequence, MODRES *modres)
{
   PDB  *pdbOut    = NULL,
        *resIn     = NULL,
        *pOut      = NULL,
        *nextResIn = NULL;
   char *fixedRes  = fixedSequence,
        pdbRes;

   /* Initialize pointer to PDB residues                                */
   resIn=pdbIn;
   
   /* Step through the fixedSequence                                    */
   for(fixedRes=fixedSequence; *fixedRes!='\0'; fixedRes++)
   {
      /* Skip chain break indications in fixedSequence                  */
      while(*fixedRes == '*') fixedRes++;
      /* Catch skipping beyond the final *                              */
      if(*fixedRes == '\0') break;

      /* Find next residue                                              */
      nextResIn = blFindNextResidue(resIn);

      /* Find the amino acid type for this residue                      */
      pdbRes = blThrone(resIn->resnam);
      if(pdbRes == 'X')
      {
         char tmpthree[8];
         blFindOriginalResType(resIn->resnam, tmpthree, modres);
         pdbRes = blThrone(tmpthree);
      }

      /* If the residue is in lower case, it's an insertion
      */
      if(islower(*fixedRes))
      {
#ifdef DEBUG
         fprintf(stderr,"INSERT: PDB %c, SEQ %c\n",
                 pdbRes, *fixedRes);
#endif
         AppendNewResidue(&pdbOut, &pOut, toupper(*fixedRes));
         if(pdbOut == NULL)
            return(NULL);
      }
      else
      {
         /* If they match, append this residue */
         if(pdbRes == *fixedRes)
         {
#if (DEBUG > 1)
         fprintf(stderr,"APPEND: PDB %c, SEQ %c\n",
                 pdbRes, *fixedRes);
#endif
            AppendThisResidue(&pdbOut, &pOut, resIn, nextResIn);
         }
         else
         {
#ifdef DEBUG
            fprintf(stderr,"FIXING: PDB %c, SEQ %c\n",
                    pdbRes, *fixedRes);
#endif
            AppendFixedResidue(&pdbOut, &pOut, resIn, nextResIn,
                               toupper(*fixedRes));
         }
         if(pdbOut == NULL)
            return(NULL);

         /* Move on to the next PDB residue */
         resIn = nextResIn;
      }
   }

   AppendRemainingAtoms(&pdbOut, &pOut, resIn);
   
   return(pdbOut);
}

