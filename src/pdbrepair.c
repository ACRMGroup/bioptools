/************************************************************************/
/**
TODO:

Implement the fix sequence code by calling routines from mutmodel

Need to fill in missing atoms

Check with other modified residues

Maybe need a local sequence alignment and/or option to accept the ATOM
residues rather than SEQRES

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
#define MAXBUFF        160
#define MAXCHAINS      240
#define MAXATINRESBASE  40
#define GAPPEN           2
#define safetoupper(x) ((islower(x))?toupper(x):(x))
#define safetolower(x) ((isupper(x))?tolower(x):(x))
#define CONECT_TOL 0.2
#define RESTYPE_PROTEIN 1
#define RESTYPE_DNA     2
#define RESTYPE_RNA     4
#define RESTYPE_HET     8
#define RESTYPE_WAT     16
#define MAXRESTYPE      16

typedef struct _restype
{
   char resnam[8],
        aa;
   int  type;
   char atnams[MAXATINRESBASE+1][8];
}  RESTYPE;

typedef struct _chaintype
{
   char chain[8];
   int  type;
   struct _chaintype *next;
}  CHAINTYPE;


/************************************************************************/
/* Globals
*/
RESTYPE gResTypes[] =
{
   {"ALA",'A',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ",""}},
   {"CYS",'C',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","SG  ",""}},
   {"ASP",'D',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","OD1 ","OD2 ",""}},
   {"GLU",'E',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","OE1 ","OE2 ",""}},
   {"PHE",'F',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ","CE1 ",
     "CE2 ","CZ  ",""}},
   {"GLY",'G',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ",""}},
   {"HIS",'H',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","ND1 ","CD2 ","CE1 ",
     "NE2 ",""}},
   {"ILE",'I',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG1 ","CG2 ","CD1 ",""}},
   {"LYS",'K',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","CE  ","NZ  ",""}},
   {"LEU",'L',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CD1 ","CD2 ",""}},
   {"MET",'M',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","SD  ","CE  ",""}},
   {"ASN",'N',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","OD1 ","ND2 ",""}},
   {"PRO",'P',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ",""}},
   {"GLN",'Q',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","OE1 ","NE2 ",""}},
   {"ARG",'R',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ","NE  ","CZ  ",
     "NH1 ","NH2 ",""}},
   {"SER",'S',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","OG  ",""}},
   {"THR",'T',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","OG1 ","CG2 ",""}},
   {"VAL",'V',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG1 ","CG2 ",""}},
   {"TRP",'W',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ","NE1 ",
     "CE2 ","CE3 ","CZ2 ","CZ3 ","CH2",""}},
   {"TYR",'Y',RESTYPE_PROTEIN,
    {"N   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ","CD2 ","CE1 ",
     "CE2 ","CZ  ","OH  ",""}},
   {"PCA",'E',RESTYPE_HET,
    {"N   ","CA  ","CB  ","CG  ","CD  ","OE  ","C   ","O   ",""}},
   {"  U",'U',RESTYPE_RNA,
    {"P   ","OP1 ","OP2 ","O5' ","C5' ","C4' ","O4' ","C3' ","O3' ",
     "C2' ","O2' ","C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","O4  ",
     "C5  ","C6  ",""}},
   {"  A",'A',RESTYPE_RNA,
    {"P   ","OP1 ","OP2 ","O5' ","C5' ","C4' ","O4' ","C3' ","O3' ",
     "C2' ","O2' ","C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","O4  ",
     "C5  ","C6  ",""}},
   {"  C",'C',RESTYPE_RNA,
    {"P   ","OP1 ","OP2 ","O5' ","C5' ","C4' ","O4' ","C3' ","O3' ",
     "C2' ","O2' ","C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","N4  ",
     "C5  ","C6  ",""}},
   {"  G",'G',RESTYPE_RNA,
    {"P   ","OP1 ","OP2 ","O5' ","C5' ","C4' ","O4' ","C3' ","O3' ",
     "C2' ","O2' ","C1' ","N9  ","C8  ","N7  ","C5  ","C6  ","O6  ",
     "N1  ","C2  ","N2  ","N3  ","C4  ",""}},
   {" DG",'G',RESTYPE_DNA,
    {"O5' ","C5' ","C4' ","O4' ","C3' ","O3' ","C2' ","C1' ","N9  ",
     "C8  ","N7  ","C5  ","C6  ","O6  ","N1  ","C2  ","N2  ","N3  ",
     "C4  ",""}},
   {" DT",'T',RESTYPE_DNA,
    {"P   ","OP1 ","OP2 ","O5' ","C5' ","C4' ","O4' ","C3' ","O3' ",
     "C2' ","C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","O4  ","C5  ",
     "C7  ","C6  ",""}},
   {" DC",'C',RESTYPE_DNA,
    {"P   ","OP1 ","OP2 ","O5' ","C5' ","C4' ","O4' ","C3' ","O3' ",
     "C2' ","C1' ","N1  ","C2  ","O2  ","N3  ","C4  ","N4  ","C5  ",
     "C6  ",""}},
   {" DA",'A',RESTYPE_DNA,
    {"P   ","OP1 ","OP2 ","O5' ","C5' ","C4' ","O4' ","C3' ","O3' ",
     "C2' ","C1' ","N9  ","C8  ","N7  ","C5  ","C6  ","N6  ","N1  ",
     "C2  ","N3  ","C4  ",""}},
   {"HOH",'O',RESTYPE_WAT,
    {"O   ",""}},
   {""  ,'\0',0,{""}}
};

   

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *trimSeqs);
void Usage(void);
PDB *RepairPDB(PDB *pdb, char *fixedSequence, MODRES *modres,
               BOOL *repaired);
void AppendThisResidue(PDB **ppdbOut, PDB **ppOut,
                       PDB *resIn,    PDB *nextResIn);
void AppendNewResidue(PDB **ppdbOut, PDB **ppOut,
                       char restype);
void AppendFixedResidue(PDB **ppdbOut, PDB **ppOut,
                        PDB *resIn,    PDB *nextResIn,
                        char restype);
BOOL AppendRemainingAtomRecords(PDB **ppdbOut, PDB **ppOut,
                                PDB *pdbIn);
char *TrimSequence(char *sequence);
int FindResType(PDB *p);
CHAINTYPE *GetChainTypes(PDB *pdb);
STRINGLIST *CreateSEQRES(PDB *pdb, CHAINTYPE *chaintypes);
char *OnethrRNA(char one);
char *OnethrDNA(char one);
STRINGLIST *myCreateSEQRES(PDB *fixedPDB);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for filling in missing records 

-  29.10.21 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   BOOL trimSequences = FALSE;
   
        
   if(ParseCmdLine(argc, argv, infile, outfile, &trimSequences))
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
            PDB        *fixedPDB,
                       *pdb;
            MODRES     *modres         = NULL;
            char       *seqresSequence = NULL,
                       *atomSequence   = NULL,
                       *fixedSequence  = NULL,
                       **seqresChains  = NULL,
                       **outchains     = NULL,
                       **atomChains    = NULL;
            int        nAtomChains,
                       len1;
            BOOL       repaired        = FALSE;
            CHAINTYPE  *chainTypes     = NULL;
#ifdef DEBUG
            CHAINTYPE  *ct             = NULL;
#endif
            
            pdb = wpdb->pdb;

            chainTypes = GetChainTypes(pdb);
#ifdef DEBUG
            for(ct=chainTypes; ct!=NULL; NEXT(ct))
            {
               printf("Chain %s Type %d\n", ct->chain, ct->type);
            }
            return(0);
#endif
            
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
               fprintf(stderr,"Error: No memory for seqresChains \
array\n");
               return(1);
            }

            /* Read MODRES and SEQRES records                           */
            modres = blGetModresWholePDB(wpdb);
            seqresSequence = blGetSeqresAsStringWholePDB(wpdb,
                                                         seqresChains,
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
            if((fixedSequence = blFixSequence(seqresSequence,atomSequence,
                                              seqresChains,
                                              atomChains,outchains,FALSE,
                                              nAtomChains,
                                              FALSE,FALSE,NULL))
               ==NULL)
               return(1);

#ifdef DEBUG
            fprintf(stderr, "Fixed sequence:\n");
            fprintf(stderr, "%s\n", fixedSequence);
#endif
            if(trimSequences)
            {
               if((fixedSequence = TrimSequence(fixedSequence))
                  ==NULL)
                  return(1);
#ifdef DEBUG
               fprintf(stderr, "Trimmed sequence:\n");
               fprintf(stderr, "%s\n", fixedSequence);
#endif
            }
            
            fixedPDB = RepairPDB(pdb, fixedSequence, modres, &repaired);
            wpdb->pdb = fixedPDB;
            if(trimSequences)
            {
               STRINGLIST *seqres         = NULL;
               seqres = myCreateSEQRES(fixedPDB);
               blReplacePDBHeader(wpdb, "SEQRES", seqres);
            }

            if(repaired)
               blBuildConectData(wpdb->pdb, CONECT_TOL);
            
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
                     BOOL *trimSeqs)
   ----------------------------------------------------------------------
*//**
   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \param[out]     *outfile    Output filename (or blank string)
   \param[out]     *trimSeqs   Trim missing residues from the end of the
                               SEQRES records
   \return                     Success

   Parse the command line

-  29.10.21 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *trimSeqs)
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
         case 't':
            *trimSeqs = TRUE;
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
   fprintf(stderr,"\nUsage: pdbrepair [-t] [in.pdb [out.pdb]]\n");
   fprintf(stderr,"       -t Trim SEQRES data for missing resiudes at \
the start or end\n");
   fprintf(stderr,"          of a chain\n");
   
   
   fprintf(stderr,"\nIf files are not specified, stdin and stdout are \
used.\n");
   fprintf(stderr,"Currently just adds missing ATOM records for atoms \
based on the\n");
   fprintf(stderr,"SEQRES records. Coordinates are set to 9999.999\n\n");
}



/************************************************************************/
void AppendNewResidue(PDB **ppdbOut, PDB **ppOut, char restype)
{
   PDB  *pOut = *ppOut;
   int  resOffset,
        atomOffset;
   char prevChain[8];
   static int newResnum = -1;

   for(resOffset=0; gResTypes[resOffset].aa != '\0'; resOffset++)
   {
      if(gResTypes[resOffset].aa == restype)
         break;
   }
   
   for(atomOffset=0;
       gResTypes[resOffset].atnams[atomOffset][0] != '\0';
       atomOffset++)
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

      /* Initialize this PDB item                                       */
      pOut->x              = 9999.999;
      pOut->y              = 9999.999;
      pOut->z              = 9999.999;
      pOut->occ            = 0.0;
      pOut->bval           = 99.0;
      pOut->access         = 0.0;
      pOut->radius         = 0.0;
      pOut->partial_charge = 0.0;

      pOut->atnum          = 0;
      pOut->resnum         = newResnum;
      pOut->formal_charge  = 0;
      pOut->nConect        = 0;
      pOut->entity_id      = 0;
      pOut->atomtype       = 0;
      
      pOut->element[0]     = gResTypes[resOffset].atnams[atomOffset][0];
      pOut->element[1]     = '\0';
      
      pOut->altpos         = ' ';
      pOut->secstr         = ' ';
      
      strcpy(pOut->record_type,
             gResTypes[resOffset].type==RESTYPE_HET?"HETATM":"ATOM  ");
#ifdef BREAKSSEQRES
      strcpy(pOut->record_type, "INSERT");
#endif
      strcpy(pOut->atnam,
             gResTypes[resOffset].atnams[atomOffset]);

      pOut->atnam_raw[0]   = ' ';
      strcpy(pOut->atnam_raw+1, gResTypes[resOffset].atnams[atomOffset]);
      pOut->atnam_raw[4]   = '\0';
      
      strcpy(pOut->resnam,      gResTypes[resOffset].resnam);
      strcpy(pOut->insert,      " ");
      strcpy(pOut->chain,       prevChain);
      strcpy(pOut->segid,       " ");

   }
   newResnum--;

   *ppOut = pOut;
}


/************************************************************************/
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


/************************************************************************/
/* Just behaves the same as AppendThisResidue
 */
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
#ifdef DEBUG
      strcpy(pOut->record_type, "FIXED ");
#endif
      
   }
   *ppOut = pOut;
}


/************************************************************************/
BOOL AppendRemainingAtomRecords(PDB **ppdbOut, PDB **ppOut, PDB *pdbIn)
{
   PDB  *pOut = *ppOut,
        *pIn  = NULL;
   BOOL repaired = FALSE;
   
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
         
         return(FALSE);
      }

      blCopyPDB(pOut, pIn);
      repaired = TRUE;
   }
   *ppOut = pOut;

   return(repaired);
}


/************************************************************************/
PDB *RepairPDB(PDB *pdbIn, char *fixedSequence, MODRES *modres,
               BOOL *repaired)
{
   PDB  *pdbOut    = NULL,
        *resIn     = NULL,
        *pOut      = NULL,
        *nextResIn = NULL;
   char *fixedRes  = fixedSequence,
        pdbRes;

   /* Assume nothing was done                                           */
   *repaired = FALSE;

   /* Initialize pointer to PDB residues                                */
   resIn=pdbIn;
   
   /* Step through the fixedSequence                                    */
   for(fixedRes=fixedSequence; *fixedRes!='\0'; fixedRes++)
   {
      /* Skip chain break indications in fixedSequence                  */
      while(*fixedRes == '*') fixedRes++;
      /* Catch skipping beyond the final *                              */
      if(*fixedRes == '\0') break;

      if(resIn!=NULL)
      {
         /* Find next residue                                           */
         nextResIn = blFindNextResidue(resIn);

         /* Find the amino acid type for this residue                   */
         pdbRes = blThrone(resIn->resnam);
      }
      else
      {
         pdbRes = '-';
      }
      
      if(pdbRes == 'X')
      {
         char tmpthree[8];
         blFindOriginalResType(resIn->resnam, tmpthree, modres);
         pdbRes = blThrone(tmpthree);
      }

      /* If the residue is in lower case, it's an insertion             */
      if(islower(*fixedRes))
      {
#ifdef DEBUG
         fprintf(stderr,"INSERT: PDB %c, SEQ %c\n",
                 pdbRes, *fixedRes);
#endif
         AppendNewResidue(&pdbOut, &pOut, toupper(*fixedRes));
         if(pdbOut == NULL)
            return(NULL);
         *repaired = TRUE;
      }
      else
      {
         /* If they match, append this residue                          */
         if(pdbRes == *fixedRes)
         {
#if (DEBUG > 2)
            fprintf(stderr,"APPEND: PDB %c, SEQ %c\n",
                    pdbRes, *fixedRes);
#endif
            AppendThisResidue(&pdbOut, &pOut, resIn, nextResIn);
         }
         else
         {
#if (DEBUG > 1)
            fprintf(stderr,"FIXING: PDB %c, SEQ %c\n",
                    pdbRes, *fixedRes);
#endif
            AppendFixedResidue(&pdbOut, &pOut, resIn, nextResIn,
                               toupper(*fixedRes));
            *repaired = TRUE;
         }
         if(pdbOut == NULL)
            return(NULL);

         /* Move on to the next PDB residue                             */
         resIn = nextResIn;
      }
   }

   if(AppendRemainingAtomRecords(&pdbOut, &pOut, resIn))
      *repaired = TRUE;
      
   if(*repaired)
   {
      blRenumAtomsPDB(pdbOut,    1);
      blRenumResiduesPDB(pdbOut, 1);
   }
   
   return(pdbOut);
}


/************************************************************************/
char *TrimSequence(char *inSeq)
{
   char *startPtr,
      *outSeq,
      *outPtr,
      *boundaryPtr,
      *endPtr;

   /* Allocate space for trimmed sequence                               */
   if((outSeq = (char *)malloc((1+strlen(inSeq))*sizeof(char))) == NULL)
      return(NULL);
   outPtr = outSeq;

   startPtr = inSeq;
   /* While we have a sequence left to process                          */
   while((*startPtr != '\0') && (boundaryPtr = strchr(startPtr, '*'))
         != NULL)
   {
      /* Skip over leading lower-case characters                        */
      while(islower(*startPtr)) startPtr++;
      /* Skip back over trailing lower-case characters                  */
      endPtr = boundaryPtr-1;
      while(islower(*endPtr) && (endPtr > startPtr)) endPtr--;

#if (DEBUG > 10)
      fprintf(stderr, "Trimmed start: %c; Trimmed end: %c\n",
              *startPtr, *endPtr);
#endif
      
      /* Copy to the output                                             */
      do {
         *outPtr++ = *startPtr++;
      } while(startPtr <= endPtr);
      *outPtr++ = '*';

      /* Skip over the *                                                */
      startPtr = boundaryPtr+1;
   }
   *outPtr = '\0';

   free(inSeq);
   return(outSeq);
}


/************************************************************************/
CHAINTYPE *GetChainTypes(PDB *pdb)
{
   CHAINTYPE *chaintypes = NULL,
             *ct         = NULL;
   PDB       *startChain = NULL,
             *nextChain  = NULL,
             *startRes   = NULL,
             *nextRes    = NULL;
   int       type        = 0;
   
   for(startChain=pdb; startChain!=NULL; startChain=nextChain)
   {
      type      = 0;
      nextChain = blFindNextChain(startChain);
      
      for(startRes=startChain; startRes!=nextChain; startRes=nextRes)
      {
         nextRes = blFindNextResidue(startRes);
         type |= FindResType(startRes);
      }

      if(ct==NULL)
      {
         INIT(chaintypes, CHAINTYPE);
         ct = chaintypes;
      }
      else
      {
         ALLOCNEXT(ct, CHAINTYPE);
      }
      if(ct==NULL)
      {
         fprintf(stderr, "Error: No memory for chaintype information\n");
         return(NULL);
      }

      ct->type = type;
      strncpy(ct->chain, startChain->chain, 8);
   }

   return(chaintypes);
}

/************************************************************************/
int FindResType(PDB *p)
{
   int i;
   for(i=0; gResTypes[i].aa != '\0'; i++)
   {
      if(!strncmp(p->resnam, gResTypes[i].resnam, 3))
      {
         return(gResTypes[i].type);
      }
   }
   return(RESTYPE_HET);
}

/************************************************************************/
/*>STRINGLIST *CreateSEQRES(PDB *pdb, CHAINTYPE *chaintypes)
   ------------------------------------
*//**
   \param[in]   *pdb   Start of PDB linked list
   \return             STRINGLIST linked list of SEQRES data

   Creates a linked list of strings representing the SEQRES data

   17.11.21  Original   By: ACRM
*/
STRINGLIST *CreateSEQRES(PDB *pdb, CHAINTYPE *chaintypes)
{
   HASHTABLE  *seqByChain = NULL;
   STRINGLIST *seqres     = NULL;
   int        nChains     = 0;
   char       **chains    = NULL,
              buffer[MAXBUFF],
              aa[8];
   CHAINTYPE  *ct         = chaintypes;
   
   if((seqByChain = blPDB2SeqXByChain(pdb))!=NULL)
   {
      if((chains = blGetPDBChainLabels(pdb, &nChains))==NULL)
      {
         return(NULL);
      }
      else
      {
         int i;

         /* Prevent blOnethr() from using the nucleic acid info         */
         gBioplibSeqNucleicAcid = 0;
         
         
         for(i=0; i<nChains; i++)
         {
            int  chainLen,
                 lineNum   = 1,
                 resNum    = 0;
            char *sequence = NULL;
            BOOL Stored    = TRUE;

            /* Skip any non ATOM chains                                 */
            while(!(ct->type | RESTYPE_PROTEIN) &&
                  !(ct->type | RESTYPE_DNA) &&
                  !(ct->type | RESTYPE_RNA))
            {
               NEXT(ct);
            }
            
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
                  if((ct->type & RESTYPE_PROTEIN))
                  {
                     sprintf(aa, "%-4s", blOnethr(sequence[resNum]));
                  }
                  else if((ct->type & RESTYPE_RNA))
                  {
                     sprintf(aa, "%-4s", OnethrRNA(sequence[resNum]));
                  }
                  else
                  {
                     sprintf(aa, "%-4s", OnethrDNA(sequence[resNum]));
                  }
                  
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
            NEXT(ct);  /* Next chaintype                                */
         }
      }
   }
   
   if(seqByChain != NULL)
      blFreeHash(seqByChain);
   
   return(seqres);
}


char *OnethrRNA(char one)
{
   int i;
   for(i=0; gResTypes[i].aa != '\0'; i++)
   {
      if((gResTypes[i].type == RESTYPE_RNA) &&
         (gResTypes[i].aa   == one))
      {
         return(gResTypes[i].resnam);
      }
   }
   return("UNK");
}

char *OnethrDNA(char one)
{
   int i;
   for(i=0; gResTypes[i].aa != '\0'; i++)
   {
      if((gResTypes[i].type == RESTYPE_DNA) &&
         (gResTypes[i].aa   == one))
      {
         return(gResTypes[i].resnam);
      }
   }
   return("UNK");
}

STRINGLIST *myCreateSEQRES(PDB *pdb)
{
   STRINGLIST *seqres = NULL;
   PDB  *res,
        *nextRes;
   char buffer[MAXBUFF],
        aa[8],
        lastChain[8];
   BOOL Stored     = TRUE;
   int  lineNum    = 1,
        chainCount = 0,
        resCount   = 0,
        chainLengths[MAXCHAINS];

   strcpy(lastChain, pdb->chain);
   resCount = 0;
   chainCount = 0;

   /* Run through the chains to find their lengths                      */
   for(res=pdb; res!=NULL; res=nextRes)
   {
      nextRes = blFindNextResidue(res);
      if(!CHAINMATCH(res->chain, lastChain))
      {
         chainLengths[chainCount++] = resCount;
         resCount = 0;
         strcpy(lastChain, res->chain);
      }
      resCount++;
   }
   chainLengths[chainCount++] = resCount;
   
   /* Now create the seqres records                                     */
   lastChain[0] = '\0';
   chainCount   = -1;
   resCount     = 0;
   for(res=pdb; res!=NULL; res=nextRes)
   {
      nextRes = blFindNextResidue(res);
      if(!CHAINMATCH(res->chain, lastChain))
      {
         /* Start of a new chain                                        */
         if(!Stored)
         {
            strcat(buffer, "\n");
            seqres = blStoreString(seqres, buffer);
            Stored = TRUE;
         }

         /* Reset chain info                                            */
         strcpy(lastChain, res->chain);
         resCount = 0;
         lineNum  = 1;
         chainCount++;

         /* Initialize SEQRES record                                    */
         sprintf(buffer, "SEQRES%4d %c%5d  ",
                 lineNum++,
                 res->chain[0],
                 chainLengths[chainCount]);
      }

      /* Save the current SEQRES if we have 13aa                        */
      if(!(resCount%13))
      {
         if(!Stored)
         {
            strcat(buffer, "\n");
            seqres = blStoreString(seqres, buffer);
            Stored = TRUE;

            /* Initialize SEQRES record                                 */
            sprintf(buffer, "SEQRES%4d %c%5d  ",
                    lineNum++,
                    res->chain[0],
                    chainLengths[chainCount]);
         }
      }

      /* Add the amino acid                                             */
      sprintf(aa, "%-4s", res->resnam);
      strcat(buffer, aa);
      Stored = FALSE;
      resCount++;
   }

   /* Store the last one                                                */
   if(!Stored)
   {
      strcat(buffer, "\n");
      seqres = blStoreString(seqres, buffer);
   }
   return(seqres);
}

