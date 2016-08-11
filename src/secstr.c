#define TESTCODE 1
/************************************************************************/
/**

   \File       secstruc.c
   
   \version    V1.1
   \date       10.07.15
   \brief      Secondary structure calculation
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 1988-2015
   \author     Dr. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College London,
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

   V1.0   19.05.99 Original, written while at Inpharmatica   By: ACRM
   V1.1   10.07.15 Modified for BiopLib
          17.02.16 xyz for coordinates now count from zero
                   CalcDihedral() also now uses arrays that count from
                   zero
          18.02.16 NUM_MC_ATOM_TYPES now counts from zero
                   MAX_NUM_ANGLES now counts from zero

*************************************************************************/
/* Doxygen
   -------
   #GROUP    Handling PDB Data
   #SUBGROUP Analyzing structures
   #FUNCTION blCalcSecStrucPDB()
   Calculate secondary structure populating the secstr field of the PDB
   structure.
*/
/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/pdb.h"
#include "bioplib/array.h"
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/angle.h"

#include "secstr.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define COORD_DIM          3    /* number of coord dimensions: x, y, z  */
#define NUM_STRAND_CHARS  26    /* number of available strand characters
                                   (the 26 letters of the alphabet)     */
#define MAX_NUM_HBOND      4    /* number of allowed H-bonds            */
#define MAX_NUM_CHN       200   /* number of allowed 'chains'; the way
                                   the code is used, it should be handled
                                   an individual chain, but it can also
                                   cope with (up to 50) multiple chains
                                   and has to allow this to deal with 
                                   apparent chain breaks resulting from 
                                   missing residues                     */

/* Used by MakeTurnsAndBridges()                                        */
#define NRULES                3
#define PARLEL                1

#define NUM_BRIDGE            4 /* Max bridge index                     */
#define NUM_STRAND_PAIR       2 /* Paired strands                       */
#define NUM_BRIDGE_PAIR       2 /* Paired strands for bridge            */

#define NUM_DIHED_DATA        3 /* number of dihedral data points       */
#define RESTYPE_PROLINE      15 /* residue type number for Proline - an
                                   offset into the KnownResidueTypes
                                   array                                */

#define NUM_MC_ATOM_TYPES     6 /* number of M/C atom types             */
#define ATOM_N                0 /* atom indexes                         */
#define ATOM_CA               1
#define ATOM_C                2
#define ATOM_O                3
#define ATOM_H                4
#define ATOM_CB               5

#define MAX_NUM_ANGLES           7 /* number of angles calculated       */
#define NUM_STANDARD_DIHEDRALS   4 /* number of standard dihedrals      */
#define DIHED_PHI                0 /* dihedral indexes                  */
#define DIHED_PSI                1
#define DIHED_OMEGA              2
#define DIHED_CHIRAL             3
#define DIHED_IMPLD              4
#define DIHED_KAPPA              5
#define DIHED_TCO                6
  
#define NUM_STRUC_SYMS         (11) /* number of structure symbols      */
#define SECSTR_IDX_ALPHAH        0  /* Symbol indexes                   */
#define SECSTR_IDX_SHEET1        1
#define SECSTR_IDX_SHEET         2
#define SECSTR_IDX_BRIDGE        3
#define SECSTR_IDX_BRIDG1        4
#define SECSTR_IDX_BRIDG2        5
#define SECSTR_IDX_THRTNH        6
#define SECSTR_IDX_PIH           7
#define SECSTR_IDX_TURN          8
#define SECSTR_IDX_BEND          9
#define SECSTR_IDX_CHISGN       10

#define NUM_HELIX_TYPE        3 /* number of helix types                */
#define HELIX_CHUNK_SIZE      4 /* Helices are made in chunks of this   */
#define THREE_TEN_CHUNK_SIZE  3 /* 3_10 helices in chunks of this       */
#define PIH_CHUNK_SIZE        5 /* Pi helices in chunks of this         */

#define BRIDGE_SEP            2 /* bridge separation                    */
#define BULGE_SIZE_L          5
#define BULGE_SIZE_S          2

#define CST                   0 /* Offsets to store HBond energies      */
#define NST                   2

#define NUM_TURN_TYPE         3 /* number of turn types                 */

#define RADIAN (180.0/3.141592) /* constant to convert RADs to degrees  */
#define NULLVAL           999.9 /* used as a NULL value                 */

#define MAX_PEPTIDE_BOND    2.5 /* maximum peptide bond length          */
#define MAX_CA_DISTANCE     5.0 /* maximum CA-CA distance               */

#define BEND_SIZE          70.0 /* mainchain atoms bent by more than this
                                   angle are flagged as "bend"          */

/* Constants used by MakeHBonds() to create hydrogen bonds              */
#define MAX_HBOND_ENERGY    -0.5
#define HBOND_ENERGY_LIMIT  -9.9 
#define HBOND_MAX_CA_DIST    8.0
#define HBOND_Q1             0.42 
#define HBOND_Q2             0.2 
#define HBOND_F            332 

/* Macros for approximate equality and inequality                       */
#define ACCURACY 0.0001
#define APPROXEQ(x, y) (ABS((x)-(y)) < ACCURACY)
#define APPROXNE(x, y) (ABS((x)-(y)) >= ACCURACY)

/* Macro to calculate atom distance                                     */
#define ATDIST(a, b)                                                     \
   sqrt(((a)[0]-(b)[0])*((a)[0]-(b)[0]) +                                \
        ((a)[1]-(b)[1])*((a)[1]-(b)[1]) +                                \
        ((a)[2]-(b)[2])*((a)[2]-(b)[2]))

/* Return the nearest integer (cast to a REAL)                          */
#define ANINT(x) ((REAL)(int)((x) + (((x) >= 0.0)?0.5:(-0.5))))

/* Free memory in blCalcSecStrPDB()                                     */
#define FREE_SECSTR_MEMORY                                               \
   do {                                                                  \
   if(gotAtom != NULL)                                                   \
      blFreeArray2D((char **)gotAtom, NUM_MC_ATOM_TYPES, seqlenMalloc);  \
   if(mcCoords != NULL)                                                  \
      blFreeArray3D((char ***)mcCoords, NUM_MC_ATOM_TYPES, seqlenMalloc, \
                    COORD_DIM);                                          \
   if(hbondEnergy != NULL)                                               \
      blFreeArray2D((char **)hbondEnergy, seqlenMalloc, MAX_NUM_HBOND);  \
   if(mcAngles != NULL)                                                  \
      blFreeArray2D((char **)mcAngles, MAX_NUM_ANGLES, seqlenMalloc);    \
   if(detailSS != NULL) free(detailSS);                                  \
   if(finalSS != NULL) free(finalSS);                                    \
   if(breakSymbol != NULL) free(breakSymbol);                            \
   if(ssTable != NULL)                                                   \
      blFreeArray2D((char **)ssTable, NUM_STRUC_SYMS, seqlenMalloc);     \
   if(residueID != NULL)                                                 \
      blFreeArray2D((char **)residueID, seqlenMalloc, 16);               \
   if(residueTypes != NULL) free(residueTypes);                          \
   if(hbond != NULL)                                                     \
      blFreeArray2D((char **)hbond, seqlenMalloc, MAX_NUM_HBOND);        \
   if(bridgePoints != NULL)                                              \
      blFreeArray2D((char **)bridgePoints, NUM_BRIDGE_PAIR,              \
                    seqlenMalloc);                                       \
   } while(0)


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
static void ExtractPDBData(PDB  *pdbStart, PDB *pdbStop, 
                           REAL ***mcCoords, BOOL **gotAtom, 
                           char **seqbcd, BOOL *caOnly, int *seqlen,
                           int *residueTypes, char *KnownResidueIndex);
static void SetPDBSecStruc(PDB *pdbStart, PDB *pdbStop, char *finalSS);
static int CountResidues(PDB *pdbStart, PDB *pdbStop);
static int FindResidueCode(char *resnam, char *KnownResidueIndex);
static void MarkBends(char **ssTable, char *detailSS, char *finalSS,
                      int seqlen);
static void SetChainEnds(char *breakSymbol, int *chainSize,
                         int *numChains, int *chainEnd,
                         int numberInChain, char **residueID,
                         int resIndex, int seqlen, BOOL verbose);
static REAL CalcDihedral(int angnum, REAL *atoma, REAL *atomb,
                         REAL *atomc, REAL *atomd);
static int FindHBondOffset(REAL energy, REAL bonde1, REAL bonde2);
static void SetHBond(REAL **hbondEnergy, int **hbond, REAL energy, 
                     int donorResIndex, int acceptorResIndex, 
                     int *bondCount, BOOL verbose);
static void FindNextPartialSheet(int startPoint, int sheetNum, 
                                 char **ssTable, int *sheetCode,
                                 int **bridgePoints, BOOL *found,
                                 int seqlen);
static void FindFirstUnsetSheet(int startPoint, char **ssTable,
                                int *sheetCode, int *startOfSheet, 
                                int *endOfSheet, int seqlen);
static void MarkHelicesInSummary(char *summ, char helixChar, 
                              char altHelixChar, char turnChar,
                              char altTurnChar, int helixSize, 
                              int seqlen);
static void MarkHelices(char **ssTable, char *detailSS, char *finalSS,
                       int helixPos, int helixSize, char helixCharacter, 
                       char altHelixCharacter, int seqlen);
static void MarkSheetsAndBridges(int seqlen, char **ssTable, 
                                 char *detailSS, char *finalSS, 
                                 char *ssSymbols, char *altSSSymbols);
static void RejectHBond(int resIndex, int *bond, REAL *bonde);
static void SetSheet(int bridgePoint, int sheetNum, int *sheetCode,
                     char **ssTable, int seqlen);
static void LabelSheets(char **ssTable, int **bridgePoints,
                        int *sheetCode, int seqlen, BOOL verbose);
static int FindNextCodeOccurrence(int strandStart, int strandEnd,
                                  int code, int **strandCode);
static void MarkTurns(char **ssTable, char *detailSS, char *finalSS,
                      char turnChar, char altTurnChar, int seqlen);
static void FindNextStrand(int **strandCode, int sheetStart, 
                           int sheetEnd, int *strandCount, 
                           int *startIndex, int lastStrand,
                           int *bestvalue);
static void FindChainBreaks(REAL ***mcCoords, BOOL **gotAtom,
                            char **residueID, char *breakSymbol,
                            int *chainSize, int *numChains, int *chainEnd,
                            BOOL caOnly, int seqlen, BOOL verbose);
static void AddHydrogens(REAL ***mcCoords, BOOL **gotAtom, int *chainSize,
                         int numChains, BOOL verbose);
static void MakeHBonds(REAL ***mcCoords, BOOL **gotAtom, int **hbond,
                       REAL **hbondEnergy, int *residueTypes,
                       int *chainEnd, int seqlen, BOOL verbose);
static void CalcMCAngles(REAL ***mcCoords, REAL **mcAngles,
                         BOOL **gotAtom, int *chainSize, int numChains,
                         BOOL caOnly, int seqlen);
static BOOL MakeTurnsAndBridges(int **hbond, char **ssTable,
                                REAL **mcAngles, int **bridgePoints,
                                int *chainEnd, int seqlen,
                                BOOL verbose);
static void MakeSummary(char **ssTable, char *detailSS, char *finalSS,
                        int seqlen);
static void SetBendResidues(REAL **mcAngles, char **ssTable, int seqlen);


/************************************************************************/
/*>int blCalcSecStrucPDB(PDB *pdbStart, PDB *pdbStop, BOOL verbose)
   -----------------------------------------------------------------
*//**
   \param[in]  *pdbStart   Start of PDB linked list
   \param[in]  *pdbStop    End of PDB linked list (NULL or pointer
                           to start of next chain)
   \param[in]  verbose     Provide informational messages
   \return                 0   - success
                           -ve - error (See SECSTR_ERR_xxxxx)
 
   Calculate secondary structure populating the ss field of the PDB
   structure.
      
-  19.05.99 Original   By: ACRM
-  27.05.99 Standard format for messages
-  09.06.99 Warning only issued if not quiet
-  10.07.15 Modified for BiopLib
-  09.03.16 Zero-basing
-  10.08.16 Completed zero-basing
*/
int blCalcSecStrucPDB(PDB *pdbStart, PDB *pdbStop, BOOL verbose)
{
   char **ssTable    = NULL,
        *detailSS    = NULL,
        *finalSS     = NULL,
        *breakSymbol = NULL,     /* Chain break symbols                 */
        **residueID  = NULL;     /* Residue identifiers (including chain
                                    ID and insert code)                 */
   int  seqlen, 
        numChains, 
        resCount,
        seqlenMalloc,
        chainSize[MAX_NUM_CHN],  /* number of residues in chain         */
        chainEnd[MAX_NUM_CHN],
        retval         = SECSTR_ERR_NOERR,
        **hbond        = NULL,
        **bridgePoints = NULL,
        *residueTypes  = NULL;   /* amino acid sequence stored as offsets
                                    into the KnownResidueIndex[] array  */

   REAL ***mcCoords    = NULL,   /* Array of mainchain coords           */
        **hbondEnergy  = NULL,
        **mcAngles     = NULL;

   BOOL **gotAtom      = NULL,   /* Flags for mainchain atoms found     */
        caOnly;

   static char KnownResidueIndex[] = "ALAASXCYSASPGLUPHEGLYHISILEXXX\
LYSLEUMETASNXXXPROGLNARGSERTHRXXXVALTRPXXXTYRGLXUNKPCAINI";

   seqlenMalloc = CountResidues(pdbStart, pdbStop);

   /* Allocate memory for arrays                                        */
   detailSS     = (char *)malloc(sizeof(char) * seqlenMalloc);
   finalSS      = (char *)malloc(sizeof(char) * seqlenMalloc);
   breakSymbol  = (char *)malloc(sizeof(char) * seqlenMalloc);
   residueTypes = (int  *)malloc(sizeof(int)  * seqlenMalloc);

   gotAtom      = (BOOL **)blArray2D(sizeof(BOOL), NUM_MC_ATOM_TYPES, 
                                     seqlenMalloc);
   hbondEnergy  = (REAL **)blArray2D(sizeof(REAL), seqlenMalloc, 
                                     MAX_NUM_HBOND);
   mcAngles     = (REAL **)blArray2D(sizeof(REAL), MAX_NUM_ANGLES, 
                                     seqlenMalloc);
   ssTable      = (char **)blArray2D(sizeof(char), NUM_STRUC_SYMS, 
                                     seqlenMalloc);
   residueID    = (char **)blArray2D(sizeof(char), seqlenMalloc, 
                                     16);
   hbond        = (int  **)blArray2D(sizeof(int), seqlenMalloc, 
                                     MAX_NUM_HBOND);
   bridgePoints = (int  **)blArray2D(sizeof(int), NUM_BRIDGE_PAIR,
                                     seqlenMalloc);
   mcCoords     = (REAL ***)blArray3D(sizeof(REAL), NUM_MC_ATOM_TYPES, 
                                      seqlenMalloc, COORD_DIM);

   /* Check all allocations succeeded                                   */
   if((detailSS     != NULL) &&
      (finalSS      != NULL) &&
      (breakSymbol  != NULL) &&
      (residueTypes != NULL) &&
      (gotAtom      != NULL) &&
      (hbondEnergy  != NULL) &&
      (mcAngles     != NULL) &&
      (ssTable      != NULL) &&
      (residueID    != NULL) &&
      (hbond        != NULL) &&
      (bridgePoints != NULL) &&
      (mcCoords     != NULL))
   {
      /* Extract the required data from the PDB linked list             */
      ExtractPDBData(pdbStart, pdbStop, mcCoords, gotAtom, residueID,
                     &caOnly, &seqlen, residueTypes, KnownResidueIndex);

      if(caOnly)  /* Secondary structure undefined - just insert '?'    */
      {
         for(resCount=0; resCount<seqlen; resCount++)
         {
            finalSS[resCount] = '?';
         }

         if(verbose)
         {
            fprintf(stderr,"Sec Struc: (warning) protein chain %c is \
CA-only. Secondary structure undefined.\n", pdbStart->chain[0]);
         }
      }
      else
      {
         /* Sets breakSymbol[], chainSize[], numChains, chainEnd        */
         FindChainBreaks(mcCoords, gotAtom, residueID, breakSymbol, 
                         chainSize, &numChains, chainEnd, caOnly, seqlen,
                         verbose);
      
         AddHydrogens(mcCoords, gotAtom, chainSize, numChains, verbose);
      
         /* Sets hbond, hbondEnergy                                     */
         MakeHBonds(mcCoords, gotAtom, hbond, hbondEnergy, residueTypes, 
                    chainEnd, seqlen, verbose);
      
         /* Sets mcAngles[]                                             */
         CalcMCAngles(mcCoords, mcAngles, gotAtom, chainSize, numChains, 
                      caOnly, seqlen);

         /* Sets ssTable[][], bridgePoints[][]                          */
         if(!MakeTurnsAndBridges(hbond, ssTable, mcAngles, bridgePoints, 
                                 chainEnd, seqlen, verbose))
         {
            FREE_SECSTR_MEMORY;
            return(SECSTR_ERR_NOMEM);
         }
      
         /* Updates ssTable[][]                                         */
         SetBendResidues(mcAngles, ssTable, seqlen);
      
         /* Sets detailSS[], finalSS[]                                  */
         MakeSummary(ssTable, detailSS, finalSS, seqlen);
      }

      /* Set the results back into the PDB linked list                  */
      SetPDBSecStruc(pdbStart, pdbStop, finalSS);
   }
   else
   {
      retval = SECSTR_ERR_NOMEM;
   }
   
   FREE_SECSTR_MEMORY;
   return(retval);
}


/************************************************************************/
/*>static void SetPDBSecStruc(PDB *pdbStart, PDB *pdbStop, char *finalSS)
   ----------------------------------------------------------------------
*//**
   \param[in,out]  *pdbStart  Start of PDB linked list
   \param[in]      *pdbStop   End of PDB linked list
   \param[in]      *finalSS   Array of secondary structure assignments

   Copies secondary structure assignment data into the PDB linked list.
   If no backbone nitrogen is present, then the secondary structure field
   is set to a '?'.
   
-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
-  09.03.16 Completed zero-basing
*/
static void SetPDBSecStruc(PDB *pdbStart, PDB *pdbStop, char *finalSS)
{
   PDB  *p, *r,
        *nextRes = NULL;
   int  i = 0;
   char ssChar;
   BOOL gotN;

   for(p=pdbStart; p!=pdbStop; p=nextRes)
   {
      nextRes = blFindNextResidue(p);
      
      /* Check we have a backbone N                                     */
      gotN = FALSE;
      for(r=p; r!=nextRes; NEXT(r))
      {
         if(!strncmp(r->atnam, "N   ", 4))
         {
            gotN = TRUE;
            break;
         }
      }

      /* Extract the secondary structure value or use '?' if the backbone
         N was not found
      */
      ssChar = ((gotN) ? finalSS[i++] : '?');

      /* Set the values into the PDB linked list                        */
      for(r=p; r!=nextRes; NEXT(r))
      {
         r->secstr = ssChar;
      }
   }
}


/************************************************************************/
/*>static int CountResidues(PDB *pdbStart, PDB *pdbStop)
   -----------------------------------------------------
*//**
   \param[in]   *pdbStart   Start of PDB linked list
   \param[in]   *pdbStop    End of PDB linked list
   \return                  Number of residues

   Counts how many residue are present in the PDB linked list

-  10.07.15 Original   By: ACRM
-  13.07.15 Modified for BiopLib
-  09.03.16 Completed zero-basing
*/
static int CountResidues(PDB *pdbStart, PDB *pdbStop)
{
   int  resCount;
   PDB  *p, *nextRes;
   
   resCount = 0;
   
   for(p=pdbStart; p!=pdbStop; p = nextRes)
   {
      resCount++;
      nextRes = blFindNextResidue(p);
   }

   return(resCount);
}


/************************************************************************/
/*>static int FindResidueCode(char *resnam, char *KnownResidueIndex)
   ------------------------------------------------------------------
*//**
   \param[in]  *resnam            Residue name
   \param[in]  *KnownResidueIndex Array of known residue names
   \return                        Numeric encoding of the residue

   Steps through the KnownResidueIndex array to find a numeric encoding
   for a given amino acid type

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
-  09.03.16 Completed zero-basing
*/
static int FindResidueCode(char *resnam, char *KnownResidueIndex)
{
   int i, codLen;
   
   codLen = strlen(KnownResidueIndex);
   
   for(i=0; i<codLen; i+=3)
   {
      if(!strncmp(resnam, KnownResidueIndex+i, 3))
         return((int)(i/3));
   }

   return(-1);
}


/************************************************************************/
/*>static void ExtractPDBData(PDB *pdbStart, PDB *pdbStop, 
                              REAL ***mcCoords, BOOL **gotAtom,
                              char **residueID, BOOL *caOnly, 
                              int *seqlen, int *residueTypes, 
                              char *KnownResidueIndex)
   -------------------------------------------------------------------
*//**
   \param[in]  *pdbStart          Start of PDB linked list
   \param[in]  *pdbStop           End of PDB linked list
   \param[out] ***mcCoords        Mainchain coordinates
   \param[out] **gotAtom          Flags to indicate whether each backbone
                                  atom type is present
   \param[out] **residueID        Residue labels - used for error messages
   \param[out] *caOnly            Is the chain CA only?
   \param[out] *seqlen            How many residues
   \param[out] *residueTypes      amino acid sequence stored as offsets
                                  into the KnownResidueIndex[] array
   \param[in]  *KnownResidueIndex Index array of residue types for 
                                  numeric encoding

   Extracts the required data from the PDB linked list

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
-  17.02.16 x,y,z for coordinates now count from 0 instead of 1
-  09.03.16 Completed zero-basing
*/
static void ExtractPDBData(PDB  *pdbStart, PDB  *pdbStop, 
                           REAL ***mcCoords, BOOL **gotAtom,
                           char **residueID, BOOL *caOnly, 
                           int *seqlen,  int *residueTypes, 
                           char *KnownResidueIndex)
{
   char buffer[MAXBUFF],
        prevInsert[8];
   int  i,
        prevResnum = (-99999999),
        atomIdx,
        caCount  = 0,
        nCount   = 0,
        resCount = (-1);
   PDB  *p;

   strncpy(prevInsert, "-", 8);

   /* Assume not C-alpha only - we check the N and CA counts and reset 
      this later if there were too few Ns for the number of CAs
   */
   *caOnly = FALSE;
   
   for(p=pdbStart; p!=pdbStop; NEXT(p))
   {
      /* If it's the start of a new residue, record the residue level
         information
      */
      if((p->resnum != prevResnum) || !INSERTMATCH(p->insert, prevInsert))
      {
         prevResnum = p->resnum;
         strcpy(prevInsert, p->insert);
         resCount++;

         for(i=0; i<NUM_MC_ATOM_TYPES; i++)
         {
            gotAtom[i][resCount] = FALSE;
         }

         sprintf(buffer,"%s.%d%s", p->chain, p->resnum, p->insert);
         strncpy(residueID[resCount], buffer, 16);
         residueTypes[resCount] = FindResidueCode(p->resnam, 
                                                  KnownResidueIndex);
      }

      /* See if it's backbone/C-beta                                    */
      if(!strncmp(p->atnam, "N   ",4))
      {
         atomIdx = ATOM_N;
         nCount++;
      }
      else if(!strncmp(p->atnam, "CA  ",4))
      {
         atomIdx = ATOM_CA;
         caCount++;
      }
      else if(!strncmp(p->atnam, "C   ",4))
      {
         atomIdx = ATOM_C;
      }
      else if(!strncmp(p->atnam, "O   ",4))
      {
         atomIdx = ATOM_O;
      }
      else if(!strncmp(p->atnam, "CB  ",4))
      {
         atomIdx = ATOM_CB;
      }
      else if(!strncmp(p->atnam, "H   ",4))
      {
         atomIdx = ATOM_H;
      }
      else
      {
         atomIdx = -1;
      }
      
      /* If it was, then we store it                                    */
      if(atomIdx != (-1))
      {
         gotAtom[atomIdx][resCount]     = TRUE;
         mcCoords[atomIdx][resCount][0] = p->x;
         mcCoords[atomIdx][resCount][1] = p->y;
         mcCoords[atomIdx][resCount][2] = p->z;
      }
   }

   if(nCount < caCount/2)
   {
      *caOnly = TRUE;
   }
   
   *seqlen = resCount+1;
}


/************************************************************************/
/*>static void MarkBends(char **ssTable, char *detailSS, char *finalSS,
                         int seqlen)
   --------------------------------------------------------------------
*//**
   \param[in]   ssTable    The sec struc table [struc symbol][res]
   \param[out]  detailSS   detailed ss assignment array
   \param[out]  finalSS    final ss assignment array
   \param[in]   seqlen     sequence length

   Enter bend symbols into summary information arrays

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
-  09.03.16 Completed zero-basing
*/
static void MarkBends(char **ssTable, char *detailSS, char *finalSS,
                      int seqlen)
{
   int resCount;
   
   for(resCount=0; resCount<seqlen; resCount++)
   {
      if(ssTable[SECSTR_IDX_BEND][resCount] != SECSTR_COIL)
      {
         if(detailSS[resCount] == SECSTR_COIL) 
            detailSS[resCount] = ssTable[SECSTR_IDX_BEND][resCount];

         if(finalSS[resCount] == SECSTR_COIL) 
            finalSS[resCount] = ssTable[SECSTR_IDX_BEND][resCount];
      }
   }
}


/************************************************************************/
/*>static void FindChainBreaks(REAL ***mcCoords, BOOL **gotAtom,
                               char **residueID, char *breakSymbol,
                               int *chainSize, int *numChains, 
                               int *chainEnd, BOOL caOnly, int seqlen,
                               BOOL verbose)
   ----------------------------------------------------------------------
*//**
   \param[in]  ***mcCoords    Mainchain coordinates
   \param[in]  **gotAtom      Flags to indicate whether each backbone
                              atom type is present
   \param[in]  **residueID    Residue ID            
   \param[in]  caOnly         Input flag            
   \param[in]  seqlen         Input seq length      
   \param[in]  verbose        Print messages
   \param[out] *breakSymbol   Array indexed by residue number indicating
                              chain breaks
   \param[out] *chainSize     Array indexed by chain number indicating
                              the number of residues in each chain
   \param[out] *numChains     Chain count
   \param[out] *chainEnd      Array indexed by chain number indicating
                              the end of each chain

   Look for chain breaks. We introduce these wherever the peptide bond
   is longer than MAX_PEPTIDE_BOND or either the C or N is missing. If 
   the chain is CA only, we introduce a break wherever the CA distance is
   greater than MAX_CA_DISTANCE

   We should be passing in individual chains, but this allows us not to
   (providing the MAX_NUM_CHN isn't exceeded) and also accounts for 
   pseudo-chain breaks as a result of missing residues.

   The code initialies the breakSymbol[] array so that all residues are
   set to coil. Where it identifies chain breaks, it calls SetChainEnds()
   which updates numChain and chainEnd, and sets the breakSymbol[] either
   side of the break to SECSTR_BREAK ('!')

-  19.05.99 Original   By: ACRM
-  27.05.99 Standard format for messages
-  13.07.15 Modified for BiopLib
-  09.03.16 Completed zero-basing
*/
static void FindChainBreaks(REAL ***mcCoords, BOOL **gotAtom,
                            char **residueID, char *breakSymbol,
                            int *chainSize, int *numChains, int *chainEnd,
                            BOOL caOnly, int seqlen, BOOL verbose)
{
   int resIndex, 
       numberInChain;

   for(resIndex=0; resIndex<seqlen; resIndex++)
   {
      breakSymbol[resIndex] = SECSTR_COIL;
   }
   numberInChain = 1;
   *numChains    = 1;
   
   for(resIndex=1; resIndex<seqlen; resIndex++)
   {
      if(caOnly)
      {
         if(ATDIST(mcCoords[ATOM_CA][resIndex-1],
                   mcCoords[ATOM_CA][resIndex]) > MAX_CA_DISTANCE)
         {
            SetChainEnds(breakSymbol, chainSize, numChains, chainEnd, 
                         numberInChain, residueID, resIndex, seqlen,
                         verbose);
            numberInChain = 1;
         }
         else
         {
            numberInChain++;
         }
      }
      else
      {
         if(gotAtom[ATOM_C][resIndex-1] && gotAtom[ATOM_N][resIndex])
         {
            if(ATDIST(mcCoords[ATOM_C][resIndex-1],
                      mcCoords[ATOM_N][resIndex]) > MAX_PEPTIDE_BOND)
            {
               SetChainEnds(breakSymbol, chainSize, numChains, chainEnd, 
                            numberInChain, residueID, resIndex, seqlen,
                            verbose);
               numberInChain = 1;
            }
            else
            {
               numberInChain++;
            }
         }
         else
         {
            SetChainEnds(breakSymbol, chainSize, numChains, chainEnd, 
                         numberInChain, residueID, resIndex, seqlen,
                         verbose);
            numberInChain = 1;
         }
      }
   }
   SetChainEnds(breakSymbol, chainSize, numChains, chainEnd, 
                numberInChain, residueID, seqlen, seqlen, verbose);
}

      
/************************************************************************/
/*>static void SetChainEnds(char *breakSymbol, int *chainSize,
                            int *numChains, int *chainEnd,
                            int numberInChain, char **residueID, 
                            int resIndex, int seqlen, BOOL verbose) 
   -------------------------------------------------------------
*//**
   \param[out]    *breakSymbol   Array indexed by residue number 
                                 indicating chain breaks
   \param[in]     *chainSize     Array indexed by chain number indicating
                                 the number of residues in each chain
   \param[in,out] *numChains     Chain count
   \param[out]    *chainEnd      Array indexed by chain number indicating
                                 the end of each chain
   \param[in]     numberInChain  Number of residues in the current chain
   \param[in]     **residueID    Residue labels - used for error messages
   \param[in]     resIndex       Current position in the list of residues
   \param[in]     seqlen         Total number of residues
   \param[in]     verbose        Give error messages

   Stores the information on where there are chain ends (or breaks in 
   the chain owing to missing residues)

   Updates numChain and chainEnd, and sets the breakSymbol[] either
   side of the break to SECSTR_BREAK ('!')

-  19.05.99 Original   By: ACRM
-  27.05.99 Standard format for messages
-  13.07.15 Modified for BiopLib
*/
static void SetChainEnds(char *breakSymbol, int *chainSize,
                         int *numChains, int *chainEnd,
                         int numberInChain, char **residueID,
                         int resIndex, int seqlen, BOOL verbose) 
{
   chainSize[*numChains-1] = numberInChain;

   if(*numChains == 1) 
   {
      chainEnd[*numChains] = numberInChain;
   }
   else
   {
      chainEnd[*numChains] = chainEnd[*numChains-1] + numberInChain;
   }

   if(resIndex < seqlen) 
   {
      breakSymbol[resIndex-1] = SECSTR_BREAK;
      breakSymbol[resIndex]   = SECSTR_BREAK;
      
      if(verbose)
      {
         fprintf(stderr,"Sec Struc: (info) Chain break between %d(%s) \
and %d(%s)\n",
                 resIndex-1, residueID[resIndex-1], resIndex, 
                 residueID[resIndex]);
      }
      
      (*numChains)++;
      if(*numChains >=  MAX_NUM_CHN) 
      {
         *numChains = MAX_NUM_CHN-1;
         if(verbose)
         {
            fprintf(stderr,"Sec Struc: (warning) Maximum number of \
chain-breaks exceeded (MAX_NUM_CHN = %d)\n", MAX_NUM_CHN);
         }
      }
   }
}

/************************************************************************/
/*>static void AddHydrogens(REAL ***mcCoords, BOOL **gotAtom, 
                            int *chainSize, int numChains, BOOL verbose)
   ---------------------------------------------------------------------
*//**
   \param[in,out] ***mcCoords  Mainchain coordinates
   \param[in]     **gotAtom    Flags to indicate whether each backbone
                               atom type is present
   \param[in]     *chainSize   Array indexed by chain number indicating
                               the number of residues in each chain
   \param[in]     numChains    Chain count
   \param[in]     verbose      Give error messages

   Create backbone hydrogens. These are placed 1A from the N such that
   N-H is parallel to the preceding C=0

-  19.05.99 Original   By: ACRM
-  27.05.99 Warning in standard format and checks quiet flag
-  13.07.15 Modified for BiopLib
-  17.02.16 x,y,z for coordinates now count from 0 instead of 1
*/
static void AddHydrogens(REAL ***mcCoords, BOOL **gotAtom, int *chainSize,
                         int numChains, BOOL verbose)
{
   int  resCount, 
        i, 
        chainStart, 
        chainNumber;
   REAL atvec[COORD_DIM], colen;
   
   chainStart = 0;
   /* For each chain                                                    */
   for(chainNumber=0; chainNumber<numChains; chainNumber++)
   {
      for(resCount=chainStart+1; 
          resCount<chainStart+chainSize[chainNumber]; 
          resCount++)
      {  /* Check we have the N, C and O                                */
         if(gotAtom[ATOM_N][resCount]   && 
            gotAtom[ATOM_C][resCount-1] &&
            gotAtom[ATOM_O][resCount-1])
         {
            colen = 0.0;
            for(i=0; i<COORD_DIM; i++)
            {
               atvec[i] = mcCoords[ATOM_C][resCount-1][i] - 
                          mcCoords[ATOM_O][resCount-1][i];
               colen += atvec[i] * atvec[i];
            }
            colen = sqrt(colen);

            for(i=0; i<COORD_DIM; i++)
            {
               atvec[i] /= colen;
            }

            for(i=0; i<COORD_DIM; i++)
            {
               mcCoords[ATOM_H][resCount][i] = 
                  mcCoords[ATOM_N][resCount][i] + atvec[i];
            }

            gotAtom[ATOM_H][resCount] = TRUE;
         }
         else
         {
            if(verbose)
            {
               fprintf(stderr,"Sec Struc: (warning) Hydrogen atom not \
generated for residue %d\n",
                       resCount+1);
            }
            gotAtom[ATOM_H][resCount] = FALSE;
         }
      }
      chainStart += chainSize[chainNumber];
   }
}



/************************************************************************/
/*>static REAL CalcDihedral(int  angnum, REAL *atoma, REAL *atomb,
                            REAL *atomc, REAL *atomd)
   ---------------------------------------------------------------
*//**
   \param[in] angnum  Which angle is being calculating
   \param[in] *atoma  First atom
   \param[in] *atomb  Second atom
   \param[in] *atomc  Third atom
   \param[in] *atomd  Fourth atom
   \return            Dihedral

   Calculate the dihedral described by 4 atoms. angnum specifies which
   dihedral is to be calculated (see the DIHED_XXXX defines) as the
   'special' ones (KAPPA (CA-CA-CA-CA) and TCO (C=O, C=O)) need special 
   treatment. 

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
-  17.02.16 x,y,z for coordinates now count from 0 instead of 1
            dotproduct[][] and dihatm[] now count from 0
*/
static REAL CalcDihedral(int  angnum, 
                         REAL *atoma,
                         REAL *atomb,
                         REAL *atomc,
                         REAL *atomd)
{
   REAL *dihatm[NUM_DIHED_DATA], 
        codist[COORD_DIM][NUM_DIHED_DATA],
        atomDistance[NUM_DIHED_DATA], 
        dotProduct[NUM_DIHED_DATA][NUM_DIHED_DATA],
        ang, 
        cosAngle; 
   int  i, j, k;

   /* It's a standard dihedral, call our normal bioplib routine         */
   if((angnum <= NUM_STANDARD_DIHEDRALS) || (angnum >= MAX_NUM_ANGLES))
   {
/*HERE*/
/*
fprintf(stderr,"%.2f %.2f %.2f\n",
        atoma[0], atoma[1], atoma[2]);
fprintf(stderr,"%.2f %.2f %.2f\n",
        atomb[0], atomb[1], atomb[2]);
fprintf(stderr,"%.2f %.2f %.2f\n",
        atomc[0], atomc[1], atomc[2]);
fprintf(stderr,"%.2f %.2f %.2f\n",
        atomd[0], atomd[1], atomd[2]);
*/
      return(RADIAN * blPhi(atoma[0], atoma[1], atoma[2],
                            atomb[0], atomb[1], atomb[2],
                            atomc[0], atomc[1], atomc[2],
                            atomd[0], atomd[1], atomd[2]));
   }

   /* A non-standard dihedral, use the special code                     */
   dihatm[0] = atoma;
   dihatm[1] = atomb;
   dihatm[2] = atomc;
   dihatm[3] = atomd;
   
   for(i=0; i<NUM_DIHED_DATA; i++)
   {
      atomDistance[i] = ATDIST(dihatm[i],dihatm[i+1]);
      for(j=0; j<COORD_DIM; j++)
      {
         codist[j][i] = dihatm[i+1][j] - dihatm[i][j];
      }
   }
   
   for(i=0; i<NUM_DIHED_DATA; i++)
   {
      for(j=0; j<NUM_DIHED_DATA; j++)
      {
         dotProduct[i][j] = 0.0;
      }
   }
   
   for(i=0; i<NUM_DIHED_DATA-1; i++)
   {
      for(j=i+1; j<NUM_DIHED_DATA; j++)
      {
         for(k=0; k<COORD_DIM; k++)
         {
            dotProduct[i][j] += (codist[k][i] * codist[k][j]);
         }
      }
   }
   
   for(i=0; i<NUM_DIHED_DATA-1; i++)
   {
      for(j=i+1; j<NUM_DIHED_DATA; j++)
      {
         if(APPROXEQ(atomDistance[i], 0.0) || 
            APPROXEQ(atomDistance[j], 0.0))
         {
            dotProduct[i][j] = 1.0;
         }
         else
         {
            dotProduct[i][j] = (dotProduct[i][j] / atomDistance[i]) / 
               atomDistance[j];
         }
      }
   }
   
   cosAngle = dotProduct[0][2];

   if(fabs(cosAngle) > 1.0) 
   {
      REAL tmp = cosAngle/fabs(cosAngle);

      cosAngle = ANINT(tmp); 
   }
   
   ang = acos(cosAngle) * RADIAN;
   
   return(ang);
}


/************************************************************************/
/*>static int FindHBondOffset(REAL energy, REAL bonde1, REAL bonde2)
   -----------------------------------------------------------------
*//**
   \param[in] energy   Energy of the new HBond
   \param[in] bonde1   Energy of 2nd best HBond
   \param[in] bonde2   Energy of best HBond
   \return             Best HBond 

   Finds which of three HBond energies is the best. Note that positive
   scores are good.

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static int FindHBondOffset(REAL energy, REAL bonde1, REAL bonde2)
{
   if(energy < bonde1)
   {
      return(1);
   }
   else if(energy < bonde2) 
   {
      return(2);
   }

   return(0);
}


/************************************************************************/
/*>static void SetHBond(REAL **hbondEnergy, int **hbond, REAL energy, 
                        int donorResIndex, int acceptorResIndex, 
                        int *bondCount, BOOL verbose)
   ------------------------------------------------------------------
*//**
   \param[out]    **hbondEnergy    HBond Energy [resnum][hbondIndex]
                                   Array of energies for the found
                                   HBonds
   \param[out]    **hbond          HBond [resnum][hbondIndex]
                                   Array of the HBonds we've found
   \param[in]     energy           Energy of a new HBond that we want
                                   to store
   \param[in]     donorResIndex    Index of the new donor
   \param[in]     acceptorResIndex Index of the new acceptor
   \param[in,out] *bondCount       Number of HBonds
   \param[in]     verbose          Print messages

   We only keep the best 2 HBonds, so find whether this one is better than
   any previous ones and slot it in. Throw away the third HBond if we
   have too many. Note that positive energies are good.

   Note that we allow 2 hbonds from C=O and 2 from N-H

-  19.05.99 Original   By: ACRM
-  27.05.99 Standard format for messages
-  13.07.15 Modified for BiopLib
*/
static void SetHBond(REAL **hbondEnergy,
                     int  **hbond,
                     REAL energy, 
                     int  donorResIndex, 
                     int  acceptorResIndex, 
                     int  *bondCount,
                     BOOL verbose)
{
   int acceptorOffset, 
       donorOffset, 
       rejectOffset;

   donorOffset    = FindHBondOffset(energy, 
                                    hbondEnergy[donorResIndex][NST], 
                                    hbondEnergy[donorResIndex][NST+1]);
   acceptorOffset = FindHBondOffset(energy, 
                                    hbondEnergy[acceptorResIndex][CST], 
                                    hbondEnergy[acceptorResIndex][CST+1]);

   if(acceptorOffset == 0 || donorOffset == 0)
   {
      if(verbose)
      {
         fprintf(stderr,"Sec Struc: (info) Third H-bond skipped. \
Donor index: %4d;  Acceptor Index: %4d; Energy: %6.2f\n", 
                 donorResIndex+1, acceptorResIndex+1, energy);
      }
   }
   else
   {
      if(hbond[donorResIndex][NST+1] != 0)
      {
         if(verbose)
         {
            fprintf(stderr,"Sec Struc: (info) Third H-bond skipped. \
Donor index: %4d; Acceptor index: %4d; Energy %6.2f\n", 
                    donorResIndex+1, 
                    hbond[donorResIndex][NST+1], 
                    hbondEnergy[donorResIndex][NST+1]);
         }
         rejectOffset = hbond[donorResIndex][NST+1];
         RejectHBond(donorResIndex+1,
                     &(hbond[rejectOffset-1][CST]),
                     &(hbondEnergy[rejectOffset-1][CST]));
         (*bondCount)--;
      }

      if(hbond[acceptorResIndex][CST+1] != 0)
      {
         if(verbose)
         {
            fprintf(stderr,"Sec Struc: (info) Third H-bond skipped. \
Donor index: %4d; Acceptor index: %4d; Energy %6.2f\n", 
                    hbond[acceptorResIndex][CST+1], 
                    acceptorResIndex+1, 
                    hbondEnergy[acceptorResIndex][CST+1]);
         }
         rejectOffset = hbond[acceptorResIndex][CST+1];
         RejectHBond(acceptorResIndex+1,
                     &(hbond[rejectOffset-1][NST]),
                     &(hbondEnergy[rejectOffset-1][NST]));
         (*bondCount)--;
      }

      if(acceptorOffset == 1)
      {
         hbond[acceptorResIndex][CST+1] = hbond[acceptorResIndex][CST];

         hbondEnergy[acceptorResIndex][CST+1] = 
            hbondEnergy[acceptorResIndex][CST];
      }

      if(donorOffset == 1)
      {
         hbond[donorResIndex][NST+1] = hbond[donorResIndex][NST];

         hbondEnergy[donorResIndex][NST+1] = 
            hbondEnergy[donorResIndex][NST];
      }

      hbond[acceptorResIndex][CST+acceptorOffset-1]       = 
         donorResIndex+1;
      hbondEnergy[acceptorResIndex][CST+acceptorOffset-1] = 
         energy;
      hbond[donorResIndex][NST+donorOffset-1]             = 
         acceptorResIndex+1;
      hbondEnergy[donorResIndex][NST+donorOffset-1]       = 
         energy;

      (*bondCount)++;
   }
}


/************************************************************************/
/*>static void FindNextPartialSheet(int startPoint, int sheetNum, 
                                    char **ssTable, int *sheetCode,
                                    int **bridgePoints, BOOL *found, 
                                    int seqlen)
   -------------------------------------------------------------------
*//**
   \param[in]  startPoint     Start of a sheet (index into residues)
   \param[in]  sheetNum       Sheet number
   \param[in]  **ssTable      The sec struc table [struc symbol][res]
   \param[out] *sheetCode     The sheet number for each residue
   \param[in]  **bridgePoints Bridges between strands
   \param[out] *found         Have we found this sheet
   \param[in]  seqlen         Length of the sequence

   Find the next partially created sheet. Set any negative sheet numbers
   to fully set and then set the partners to partially set if they are
   not already fully set.

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void FindNextPartialSheet(int  startPoint, 
                                 int  sheetNum, 
                                 char **ssTable,
                                 int  *sheetCode,
                                 int  **bridgePoints,
                                 BOOL *found,
                                 int  seqlen)
{      
   int resCount, 
       startOfSheet, 
       endOfSheet, 
       i;
   
   *found   = FALSE;
   resCount = startPoint - 1;
   
   while (sheetCode[resCount] != (-sheetNum))
   {
      resCount++;
      if(resCount >= seqlen)
         return;
   }
   
   startOfSheet = resCount;
   *found = TRUE;
   
   while(sheetCode[resCount] == (-sheetNum))
   { 
      resCount++;
      if(resCount >= seqlen) 
         break;
   }
   
   endOfSheet = resCount-1;
   
   for(resCount = startOfSheet; resCount <= endOfSheet; resCount++)
   {
      sheetCode[resCount] = sheetNum;
      for(i=0; i<NUM_BRIDGE_PAIR; i++)
      {
         if(bridgePoints[i][resCount] != 0)
         {
            if(sheetCode[bridgePoints[i][resCount]-1] == 0)
            {
               SetSheet(bridgePoints[i][resCount], sheetNum, sheetCode, 
                        ssTable, seqlen);
            }
         }
      }
   }
}


/************************************************************************/
/*>static void FindFirstUnsetSheet(int startPoint, char **ssTable,
                                   int *sheetCode, int *startOfSheet, 
                                   int *endOfSheet, int seqlen)
   ---------------------------------------------------------------
*//**
   \param[in]  startPoint    Start of sheet
   \param[in]  **ssTable     The sec struc table [struc symbol][res]
   \param[in]  *sheetCode    The sheet number for each residue
   \param[out] *startOfSheet Start of a strand
   \param[out] *endOfSheet   End of a strand
   \param[in]  seqlen        Length of the sequence

   Work through the residues from startPoint to find the first unset 
   sheet

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void FindFirstUnsetSheet(int startPoint, char **ssTable,
                                int *sheetCode, int *startOfSheet, 
                                int *endOfSheet, int seqlen)
{
   int   resCount;
   
   *startOfSheet = 0;
   resCount      = startPoint-1;
   
   while((ssTable[SECSTR_IDX_SHEET][resCount] == SECSTR_COIL) || 
         (sheetCode[resCount] != 0)) 
   {
      resCount++;
      if(resCount >= seqlen) 
         return;
   }
   *startOfSheet = resCount+1;
   
   do
   {
      resCount++;
   }  while((resCount <= seqlen-1) && 
            (ssTable[SECSTR_IDX_SHEET][resCount] != SECSTR_COIL));

   *endOfSheet = resCount;
}


      
/************************************************************************/
/*>static void MarkHelicesInSummary(char *summ, char helixChar, 
                                 char altHelixChar, char turnChar, 
                                 char altTurnChar, int helixSize, 
                                 int seqlen)
   -----------------------------------------------------------------------
*//**
   \param[in,out] *summ        Secondary structure string
   \param[in]     helixChar    Helix character (G or I)
   \param[in]     altHelixChar Helix character (g or i)
   \param[in]     turnChar     Turn character (T)
   \param[in]     altTurnChar  Turn character (t)
   \param[in]     helixSize    Length of helix
   \param[in]     seqlen       Length of the sequence

   Mark 3_{10} and Pi helix residues in the summary list

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void MarkHelicesInSummary(char *summ, char helixChar, 
                                 char altHelixChar, char turnChar,
                                 char altTurnChar, int helixSize, 
                                 int seqlen)
{
   int resCount, 
       posInHelix,
       startPos = 0;
      
   for(resCount=1; resCount<seqlen; resCount++)
   {
      if((summ[resCount] == helixChar) || 
         (summ[resCount] == altHelixChar))
      {
         if(startPos == 0) 
            startPos = resCount;
      }
      else
      {
         if(startPos != 0)
         {
            if((resCount - startPos) < helixSize)
            {
               for(posInHelix = startPos; 
                   posInHelix <= (resCount - 1); 
                   posInHelix++)
               {
                  if(summ[posInHelix] == helixChar)
                  {
                     summ[posInHelix] = turnChar;
                  }
                  else
                  {
                     summ[posInHelix] = altTurnChar;
                  }
               }
            }
            startPos = 0;
         }
      }
   }
   
   if(startPos != 0)
   {
      if((seqlen - startPos + 1) < helixSize)
      {
         for(posInHelix=startPos; posInHelix <= seqlen; posInHelix++)
         {
            if(summ[resCount] == helixChar)
            {
               summ[posInHelix] = turnChar;
            }
            else
            {
               summ[posInHelix] = altTurnChar;
            }
         }
      }
   }
}


/************************************************************************/
/*>static void MarkHelices(char **ssTable, char *detailSS, char *finalSS,
                          int helixPos, int helixSize, 
                          char helixCharacter, char altHelixCharacter, 
                          int seqlen)
   ---------------------------------------------------------------------
*//**
   \param[in]  **ssTable         The sec struc table [struc symbol][res]
   \param[out] *detailSS         detailed ss assignment array
   \param[out] *finalSS          final ss assignment array
   \param[in]  helixPos          start of a helix
   \param[in]  helixSize         length of a helix
   \param[in]  helixCharacter    helix character (H)
   \param[in]  altHelixCharacter alt helix character (h)
   \param[in]  seqlen            Sequence length

   Examine the secondary structure bonding table to set any alpha 
   helices in the summary secondary structure arrays

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void MarkHelices(char **ssTable, char *detailSS, char *finalSS,
                       int helixPos, int helixSize, char helixCharacter, 
                       char altHelixCharacter, int seqlen)
{
   int resCount, 
       posInHelix;

   for(resCount=1; resCount < (seqlen - helixSize); resCount++)
   {
      if(ssTable[helixPos][resCount] == SECSTR_BEND_START ||
         ssTable[helixPos][resCount] == SECSTR_BEND_BOTH)
      {
         if(ssTable[helixPos][resCount-1] == SECSTR_BEND_START ||
            ssTable[helixPos][resCount-1] == SECSTR_BEND_BOTH)
         {
            for(posInHelix = 0; posInHelix < helixSize; posInHelix++)
            {
               if(detailSS[resCount+posInHelix] == SECSTR_COIL)
               {
                  detailSS[resCount+posInHelix] = helixCharacter;
               }

               if(finalSS[resCount+posInHelix] == SECSTR_COIL || 
                  finalSS[resCount+posInHelix] == altHelixCharacter)
               {
                  finalSS[resCount+posInHelix] = helixCharacter;
               }
               
            }

            if(finalSS[resCount-1] == SECSTR_COIL)
            {
               finalSS[resCount-1] = altHelixCharacter;
            }

            if(finalSS[resCount+helixSize] == SECSTR_COIL)
            {
               finalSS[resCount+helixSize] = altHelixCharacter;
            }
         }
      }
   }
}


/************************************************************************/
/*>static void CalcMCAngles(REAL ***mcCoords, REAL **mcAngles, 
                            BOOL **gotAtom, int *chainSize, int numChains,
                            BOOL caOnly, int seqlen)
   -----------------------------------------------------------------------
*//**
   \param[in]  ***mcCoords  Mainchain coordinates
   \param[out] **mcAngles   Mainchain torsion angles
   \param[in]  **gotAtom    Flags to indicate whether each backbone
                            atom type is present
   \param[in]  *chainSize   Array of chain (segment) sizes
   \param[in]  numChains    Number of chains (segments)
   \param[in]  caOnly       Sequence is C-alpha only
   \param[in]  seqlen       Sequence length

   Calculate main chain dihedral angles if all the atoms are present

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void CalcMCAngles(REAL ***mcCoords, REAL **mcAngles,
                         BOOL **gotAtom, int *chainSize, int numChains,
                         BOOL caOnly, int seqlen)
{
   static int angleIndices[MAX_NUM_ANGLES][3] = 
      {
         {2,  0, DIHED_PHI},
         {1, -1, DIHED_PSI}, 
         {1, -1, DIHED_OMEGA}, 
         {2, -2, DIHED_CHIRAL}, 
         {1,  0, DIHED_IMPLD}, 
         {3, -2, DIHED_KAPPA}, 
         {2,  0, DIHED_TCO}
      };
   static int angleAtoms[MAX_NUM_ANGLES][4] = 
      {
         {ATOM_C,  ATOM_N,  ATOM_CA, ATOM_C},
         {ATOM_N,  ATOM_CA, ATOM_C,  ATOM_N}, 
         {ATOM_CA, ATOM_C,  ATOM_N,  ATOM_CA}, 
         {ATOM_CA, ATOM_CA, ATOM_CA, ATOM_CA}, 
         {ATOM_CA, ATOM_N,  ATOM_C,  ATOM_CB}, 
         {ATOM_CA, ATOM_CA, ATOM_CA, ATOM_CA}, 
         {ATOM_C,  ATOM_O,  ATOM_C,  ATOM_O}
      };
   static int angleOffsets[MAX_NUM_ANGLES][4] =
      {
         {-2, -1, -1, -1},
         {-1, -1, -1,  0}, 
         {-1, -1,  0,  0}, 
         {-2, -1,  0,  1}, 
         {-1, -1, -1, -1}, 
         {-3, -1, -1,  1}, 
         {-1, -1, -2, -2}
      };
   static int  numAtomsRequired[MAX_NUM_ANGLES] =
      {2, 2, 2, 4, 1, 5, 2};
   
   int resCount, 
       angleIndex, 
       chainStart, 
       chainCount, 
       firstAngle,
       finalAngle;
   
   firstAngle = 0;
   finalAngle = MAX_NUM_ANGLES-1;

   if(caOnly)
   {
      firstAngle = DIHED_CHIRAL;
      finalAngle = DIHED_KAPPA;
   }
   
   for(angleIndex=firstAngle; angleIndex<=finalAngle; angleIndex++)
   {
      chainStart = 0;

      for(chainCount=0; chainCount<numChains; chainCount++)
      {
         /* If too few residues in chain, create extra NULL coordinates */
         if(chainSize[chainCount] < numAtomsRequired[angleIndex])
         {
            for(resCount =  chainStart; 
                resCount < (chainStart + chainSize[chainCount]); 
                resCount++)
            {
               mcAngles[angleIndices[angleIndex][2]][resCount] = NULLVAL;
            }
         }
         else
         {
            for(resCount =  (chainStart + angleIndices[angleIndex][0]); 
                resCount <= (chainStart + chainSize[chainCount] + 
                             angleIndices[angleIndex][1]); 
                resCount++)
            {
               if(gotAtom[angleAtoms[angleIndex][0]]
                          [resCount+angleOffsets[angleIndex][0]] &&
                  gotAtom[angleAtoms[angleIndex][1]]
                          [resCount+angleOffsets[angleIndex][1]] &&
                  gotAtom[angleAtoms[angleIndex][2]]
                          [resCount+angleOffsets[angleIndex][2]] &&
                  gotAtom[angleAtoms[angleIndex][3]]
                          [resCount+angleOffsets[angleIndex][3]])
               {



                  mcAngles[angleIndices[angleIndex][2]][resCount-1] = 
                     CalcDihedral(angleIndex,
                                  mcCoords[angleAtoms[angleIndex][0]]
                                          [resCount+
                                           angleOffsets[angleIndex][0]],
                                  mcCoords[angleAtoms[angleIndex][1]]
                                          [resCount+
                                           angleOffsets[angleIndex][1]],
                                  mcCoords[angleAtoms[angleIndex][2]]
                                          [resCount+
                                           angleOffsets[angleIndex][2]],
                                  mcCoords[angleAtoms[angleIndex][3]]
                                          [resCount+
                                           angleOffsets[angleIndex][3]])
                     ;
               }
               else
               {
                  mcAngles[angleIndices[angleIndex][2]][resCount-1] = 
                     NULLVAL;
               }
            }

            for(resCount = chainStart; 
                resCount < chainStart + angleIndices[angleIndex][0] - 1; 
                resCount++)
            {
               mcAngles[angleIndices[angleIndex][2]][resCount] = NULLVAL;
            }

            for(resCount = chainStart + 
                           chainSize[chainCount] + 
                           angleIndices[angleIndex][1];
                resCount < (chainStart + chainSize[chainCount]);
                resCount++)
            {
               mcAngles[angleIndices[angleIndex][2]][resCount] = 
                  NULLVAL;
            }
         }
         chainStart += chainSize[chainCount];
      }
   }
}


/************************************************************************/
/*>static void SetBendResidues(REAL **mcAngles, char **ssTable, 
                               int seqlen)
   ------------------------------------------------------------
*//**
   \param[in]  **mcAngles Mainchain torsion angles
   \param[out] **ssTable  The sec struc table [struc symbol][res]
   \param[in]  seqlen     The sequence length

   Find any residues where the mainchain is bent by more than BEND_SIZE

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void SetBendResidues(REAL **mcAngles, char **ssTable, int seqlen)
{
   int  resCount;

   for(resCount=0; resCount<seqlen; resCount++)
   {
      if((mcAngles[DIHED_KAPPA][resCount] > BEND_SIZE) &&
          APPROXNE(mcAngles[DIHED_KAPPA][resCount], NULLVAL)) 
      {
         ssTable[SECSTR_IDX_BEND][resCount] = SECSTR_BEND;
      }
   }
}


/************************************************************************/
/*>static void MakeHBonds(REAL ***mcCoords, BOOL **gotAtom, int **hbond,
                          REAL **hbondEnergy, int *residueTypes, 
                          int *chainEnd, int seqlen, BOOL verbose)
   ---------------------------------------------------------------------
*//**
   \param[in]  ***mcCoords   Mainchain coordinates
   \param[in]  **gotAtom     Flags to indicate whether each backbone
                             atom type is present
   \param[out] **hbond       HBond [resnum][hbondIndex]
                             Array of the HBonds we've found
   \param[out] **hbondEnergy HBond Energy [resnum][hbondIndex]
                             Array of energies for the found
                             HBonds
   \param[in]  *residueTypes amino acid sequence stored as offsets
                             into the KnownResidueIndex[] array
   \param[in]  *chainEnd     Array indexed by chain number indicating
                             the end of each chain
   \param[in]  seqlen        Sequence length
   \param[in]  verbose       Print messages

   Identify mainchain hydrogen bonds. We allow each residue to make 2
   HBonds from C=O and 2 from N-H. We use the Kabsch and Sander energy 
   formula to define the presence of an HBond which corresponds to a
   maximum distance of 5.2A for linearly arranged atoms with a maximum 
   angle of 63 degrees:

                     /                       \
                     |  1     1     1     1  |
   E = q  * q  * f * | --- + --- - --- - --- |
        1    2       | d     d     d     d   |
                     \  ON    CH    OH    CN /

   The calculations are optimised by making a distance check and residues
   must have at least 1 intervening residue. We also skip prolines as
   donors!

-  19.05.99 Original   By: ACRM
-  27.05.99 Standard format for messages
-  13.07.15 Modified for BiopLib
*/
static void MakeHBonds(REAL ***mcCoords, BOOL **gotAtom, int **hbond,
                       REAL **hbondEnergy, int *residueTypes,
                       int *chainEnd, int seqlen, BOOL verbose)
{
   int  otherRes, 
        i, 
        nbonds, 
        firstChain, 
        otherChain, 
        resCount;
   REAL distON, 
        distOH, 
        distCH, 
        distCN, 
        energy,
        caDist;
   

   for(resCount=0; resCount<seqlen; resCount++)
   {
      for(i=0; i<MAX_NUM_HBOND; i++)
      {
         hbond[resCount][i] = 0;
         hbondEnergy[resCount][i] = 0.0;
      }
   }
   nbonds     = 0;
   firstChain = 1;

   for(resCount=0; resCount<seqlen; resCount++)
   {
      if(resCount > chainEnd[firstChain]) 
         firstChain++;
      
      /* Check both N and H are present                                 */
      if(gotAtom[ATOM_N][resCount] && 
         gotAtom[ATOM_H][resCount] &&
         residueTypes[resCount] != RESTYPE_PROLINE)
      {
         otherChain = 1;
         for(otherRes=0; otherRes<seqlen; otherRes++)
         {
            if(otherRes >= chainEnd[otherChain]) 
               otherChain++;

            if(((abs(resCount - otherRes) == 1) && 
                (firstChain != otherChain))     ||  
               abs(resCount - otherRes) >= 2) 
            {
               if(gotAtom[ATOM_CA][resCount] && 
                  gotAtom[ATOM_CA][otherRes]) 
               {
                  caDist = ATDIST(mcCoords[ATOM_CA][otherRes],
                                  mcCoords[ATOM_CA][resCount]);
                  if(caDist < HBOND_MAX_CA_DIST) 
                  {
                     if(gotAtom[ATOM_C][otherRes] && 
                        gotAtom[ATOM_O][otherRes]) 
                     {
                        distON = ATDIST(mcCoords[ATOM_O][otherRes],
                                        mcCoords[ATOM_N][resCount]);
                        distOH = ATDIST(mcCoords[ATOM_O][otherRes],
                                        mcCoords[ATOM_H][resCount]);
                        distCH = ATDIST(mcCoords[ATOM_C][otherRes],
                                        mcCoords[ATOM_H][resCount]);
                        distCN = ATDIST(mcCoords[ATOM_C][otherRes],
                                        mcCoords[ATOM_N][resCount]);

                        if(APPROXEQ(distON,0.0) || 
                           APPROXEQ(distOH,0.0) ||
                           APPROXEQ(distCH,0.0) || 
                           APPROXEQ(distCN,0.0)) 
                        {
                           if(verbose)
                           {
                              fprintf(stderr,"Sec Struc: (warning) \
Coincident atoms in hydrogen bonding, donor %4d acceptor %4d\n", 
                                      resCount+1, otherRes+1);
                           }
                        }
                        else
                        {
                           energy = HBOND_Q1 * HBOND_Q2 * HBOND_F * 
                              (1.0/distON + 
                               1.0/distCH - 
                               1.0/distOH - 
                               1.0/distCN);

                           if(energy < HBOND_ENERGY_LIMIT) 
                           {
                              if(verbose)
                              {
                                 fprintf(stderr,"Sec Struc: (warning) \
Atom indices %d and %d too close O-N: %8.3f C-H: %8.3f O-H: %8.3f \
C-N: %8.3f\n", 
                                         resCount+1, otherRes+1, 
                                         distON, distCH, distOH, distCN);
                              }
                              energy = HBOND_ENERGY_LIMIT;
                           }

                           if(energy < MAX_HBOND_ENERGY) 
                           {
                              SetHBond(hbondEnergy, hbond, energy, 
                                       resCount, otherRes, &nbonds,
                                       verbose);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   if(verbose)
   {
      fprintf(stderr,"Sec Struc: (info) Total Number of H-bonds: %5d\n",
              nbonds);
   }
}


/************************************************************************/
/*>static void MakeSummary(char **ssTable, char *detailSS, char *finalSS,
                           int  seqlen)
   ----------------------------------------------------------------------
*//**
   \param[in]  **ssTable The sec struc table [struc symbol][res]
   \param[out] *detailSS detailed ss assignment array
   \param[out] *finalSS  final ss assignment array
   \param[in]  seqlen    the sequence length

   Create the secodary structure summary from the secondary structure
   bonding table using Kabsch/Sander priority rules

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void MakeSummary(char **ssTable,
                        char *detailSS,
                        char *finalSS,
                        int  seqlen)
{
   static char ssSymbols[]    = "H EB  GITS+";
   static char altSSSymbols[] = "h eb  gits+";

   int resCount;
   
   for(resCount = 0; resCount < seqlen; resCount++)
   {
      detailSS[resCount] = SECSTR_COIL;
      finalSS[resCount]  = SECSTR_COIL;
   }

   /* Alpha helices                                                     */
   MarkHelices(ssTable, detailSS, finalSS, SECSTR_IDX_ALPHAH, 
              HELIX_CHUNK_SIZE, ssSymbols[SECSTR_IDX_ALPHAH],
              altSSSymbols[SECSTR_IDX_ALPHAH], seqlen);

   /* Sheets and bridges                                                */
   MarkSheetsAndBridges(seqlen, ssTable, detailSS, finalSS, ssSymbols, 
                        altSSSymbols);

   /* 3_10 helices                                                      */
   MarkHelices(ssTable, detailSS, finalSS, SECSTR_IDX_THRTNH, 
              THREE_TEN_CHUNK_SIZE, ssSymbols[SECSTR_IDX_THRTNH], 
              altSSSymbols[SECSTR_IDX_THRTNH], seqlen);

   /* Pi helices                                                        */
   MarkHelices(ssTable, detailSS, finalSS, SECSTR_IDX_PIH, PIH_CHUNK_SIZE,
               ssSymbols[SECSTR_IDX_PIH], altSSSymbols[SECSTR_IDX_PIH], 
               seqlen);

   /* Move to summary information                                       */
   MarkHelicesInSummary(detailSS, ssSymbols[SECSTR_IDX_THRTNH], 
                     altSSSymbols[SECSTR_IDX_THRTNH], 
                     ssSymbols[SECSTR_IDX_TURN], 
                     altSSSymbols[SECSTR_IDX_TURN],
                     THREE_TEN_CHUNK_SIZE, seqlen);
   MarkHelicesInSummary(finalSS, ssSymbols[SECSTR_IDX_THRTNH], 
                     altSSSymbols[SECSTR_IDX_THRTNH], 
                     ssSymbols[SECSTR_IDX_TURN], 
                     altSSSymbols[SECSTR_IDX_TURN],
                     THREE_TEN_CHUNK_SIZE, seqlen);
   MarkHelicesInSummary(detailSS, ssSymbols[SECSTR_IDX_PIH], 
                     altSSSymbols[SECSTR_IDX_PIH], 
                     ssSymbols[SECSTR_IDX_TURN], 
                     altSSSymbols[SECSTR_IDX_TURN], PIH_CHUNK_SIZE, 
                     seqlen);
   MarkHelicesInSummary(finalSS, ssSymbols[SECSTR_IDX_PIH], 
                     altSSSymbols[SECSTR_IDX_PIH], 
                     ssSymbols[SECSTR_IDX_TURN], 
                     altSSSymbols[SECSTR_IDX_TURN], PIH_CHUNK_SIZE, 
                     seqlen);

   /* Turns                                                             */
   MarkTurns(ssTable, detailSS, finalSS, ssSymbols[SECSTR_IDX_TURN], 
             altSSSymbols[SECSTR_IDX_TURN], seqlen);

   /* Bends                                                             */
   MarkBends(ssTable, detailSS, finalSS, seqlen);
}


/************************************************************************/
/*>static BOOL MakeTurnsAndBridges(int **hbond, char **ssTable, 
                                   REAL **mcAngles, int **bridgePoints,
                                   int *chainEnd, int seqlen,BOOL verbose)
   -----------------------------------------------------------------------
*//**
   \param[in]  **hbond        HBond [resnum][hbondIndex]
                              Array of the HBonds we've found
   \param[out] **ssTable      The sec struc table [struc symbol][res] 
   \param[in]  **mcAngles     Mainchain torsion angles
   \param[out] **bridgePoints Bridges between strands
   \param[in]  *chainEnd      Array indexed by chain number indicating
                              the end of each chain
   \param[in]  seqlen         The sequence length
   \param[in]  verbose        Print messages
   \return                    Memory allocation OK

   Identify turns and bridges

   TODO: This code needs better commenting and breaking down into smaller
   subroutines

-  19.05.99 Original   By: ACRM
-  27.05.99 Standard format for messages
-  16.06.99 Initialise lastStrand to 0
-  13.07.15 Modified for BiopLib
-  09.08.16 Zero-basing
*/
static BOOL MakeTurnsAndBridges(int **hbond, char **ssTable,
                                REAL **mcAngles, int **bridgePoints,
                                int *chainEnd, int seqlen,
                                BOOL verbose)
{
   static int  turnSize[NUM_TURN_TYPE]  = {3, 4, 5};
   static char turnChar[NUM_TURN_TYPE]  = {'3', '4', '5'};
   static int  turnIndex[NUM_TURN_TYPE] = {SECSTR_IDX_THRTNH, 
                                           SECSTR_IDX_ALPHAH, 
                                           SECSTR_IDX_PIH};
   static int  baseOfSt[NRULES][4] =
      {{-1,  0,  0,  1}, 
       { 0,  0,  0,  0}, 
       {-1, -1, -2,  1}
      };
   static int  baseOffsets[NRULES][2] =
      {{2, -1}, 
       {1,  0}, 
       {2, -1}
      };
   static int  bridgeOffsets[NRULES] = 
      {1, -1, -1};

   static char upperCaseLetters[] = 
      " ABCDEFGHIJKLMNOPQRSTUVWXYZ";
   static char lowerCaseLetters[] = 
      " abcdefghijklmnopqrstuvwxyz";
   static char nstchp[] = 
      " ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890";

   static int bridgeRules[2][NRULES] = 
      {{-1,  1, -1}, 
       { 1,  1, -1}
      };
   
   int  resCount, 
        bondCount, 
        turnCount, 
        i, 
        ruleCount, 
        bondPointer, 
        bridgeIndex, 
        bridgePointer, 
        theBridge, 
        bridgeDirection, 
        theNewBridge, 
        resCountInSS, 
        bridgePoint, 
        newBridgePoint, 
        strandNumber, 
        startOfSheet, 
        endOfSheet, 
        strandCount, 
        charCode, 
        strandCharOffset, 
        thisStrandOffset,
        nextStrandOffset, 
        bulgeCount, 
        bridgeCount, 
        strandOffset, 
        theBridgePoint, 
        bridgeArrayIndex,
        bestValue, 
        firstChain,
        seqlenMalloc,
        newBridgeIndex = 0,
        lastStrand = 0, 
        **bridge     = NULL,
        **strandCode = NULL,
        *sheetCode   = NULL;
   char ch, 
        strandChar, 
        bridgeChar;

   /* Allocate memory for arrays                                        */
   seqlenMalloc = seqlen;

   bridge     = (int **)blArray2D(sizeof(int), NUM_BRIDGE, seqlenMalloc);
   strandCode = (int **)blArray2D(sizeof(int), NUM_STRAND_PAIR, 
                                  seqlenMalloc);
   sheetCode  = (int *)malloc(sizeof(int) * seqlenMalloc);

   if((bridge == NULL) || (strandCode == NULL) || (sheetCode == NULL))
   {
      if(bridge     != NULL)
         blFreeArray2D((char **)bridge,     NUM_BRIDGE, seqlenMalloc);
      if(strandCode != NULL)
         blFreeArray2D((char **)strandCode, NUM_STRAND_PAIR, 
                       seqlenMalloc);
      if(sheetCode  != NULL)
         free(sheetCode);

      return(FALSE);
   }
      
   /* Initialize everything as being coil                               */
   for(resCount = 0; resCount < seqlen; resCount++)
   {
      for(i = 0; i<NUM_STRUC_SYMS; i++)
      {
         ssTable[i][resCount] = SECSTR_COIL;
      }
   }
   
   firstChain = 1;
   
   for(resCount=0; resCount<seqlen; resCount++)
   {
      if(resCount >= chainEnd[firstChain]) 
         firstChain++;
      
      for(bondCount=0; bondCount<2; bondCount++)
      {
         for(turnCount=0; turnCount<NUM_TURN_TYPE; turnCount++)
         {
            if(hbond[resCount][bondCount] != 0) 
            {
               if((hbond[resCount][bondCount]-resCount-1 == 
                    turnSize[turnCount]) && 
                   (hbond[resCount][bondCount] <= chainEnd[firstChain]))
               {
                  if(ssTable[turnIndex[turnCount]][resCount] == 
                     SECSTR_BEND_END)
                  {
                     ssTable[turnIndex[turnCount]][resCount] = 
                        SECSTR_BEND_BOTH;
                  }
                  else
                  {
                     ssTable[turnIndex[turnCount]][resCount] = 
                        SECSTR_BEND_START;
                  }
                  ssTable[turnIndex[turnCount]]
                         [hbond[resCount][bondCount]-1] = SECSTR_BEND_END;
                  
                  for(i = resCount + 1; 
                      i < resCount + turnSize[turnCount]; 
                      i++)
                  {
                     if(ssTable[turnIndex[turnCount]][i] == SECSTR_COIL)
                     {
                        ssTable[turnIndex[turnCount]][i] = 
                           turnChar[turnCount];
                     }
                     
                  }
                  for(i = resCount; 
                      i<= resCount + turnSize[turnCount]; 
                      i++)
                  {
                     ssTable[SECSTR_IDX_TURN][i] = SECSTR_TURN;
                  }
               }
            }
         }
      }
   }
   
   /* Initialize bridge counts to zero                                  */
   for(resCount = 0; resCount < seqlen; resCount++)
   {
      for(bondCount = 0; bondCount < 4; bondCount++)
      {
         bridge[bondCount][resCount] = 0;
      }
   }
   
   for(ruleCount = 0; ruleCount<NRULES; ruleCount++)
   {
      for(resCount = baseOffsets[ruleCount][0]; 
          resCount <= (seqlen + baseOffsets[ruleCount][1]); 
          resCount++)
      {
         int bond_i;
         
         for(bond_i = 0; bond_i < 2; bond_i++)
         {
            bondPointer =
               hbond[resCount+baseOfSt[ruleCount][0]-1][bond_i];

            if((ruleCount <= PARLEL) || (bondPointer > resCount))
            {
               bridgePointer = bondPointer + baseOfSt[ruleCount][1];

               if((bondPointer > -baseOfSt[ruleCount][2]) && 
                  (abs(resCount - bridgePointer) > BRIDGE_SEP))
               {
                  int bond_j;
                  
                  for(bond_j = 0; bond_j < 2; bond_j++)
                  {
                     if(hbond[bondPointer+baseOfSt[ruleCount][2]-1]
                             [bond_j] == 
                        resCount+baseOfSt[ruleCount][3])
                     {
                        int altBridgeIndex;

                        bridgeIndex = 0;

                        if(bridge[bridgeIndex][resCount-1] != 0) 
                           bridgeIndex = 1;

                        bridge[bridgeIndex][resCount-1] = 
                           bridgeOffsets[ruleCount] * bridgePointer;

                        bridge[bridgeIndex+2][resCount-1] = 
                           bridgeRules[0][ruleCount];

                        altBridgeIndex = 1;

                        if(bridge[altBridgeIndex-1][bridgePointer-1] !=0) 
                           altBridgeIndex = 2;

                        bridge[altBridgeIndex-1][bridgePointer-1] = 
                           bridgeOffsets[ruleCount] * resCount;

                        bridge[altBridgeIndex+2-1][bridgePointer-1] = 
                           bridgeRules[1][ruleCount];
                     }
                  }
               }
            }
         }
      }
   }

   for(resCount = 0; resCount < seqlen; resCount++)
   {
      if(APPROXNE(mcAngles[DIHED_CHIRAL][resCount], NULLVAL)) 
      {
         if(mcAngles[DIHED_CHIRAL][resCount] < 0.0)
         {
            ssTable[SECSTR_IDX_CHISGN][resCount] = '-';
         }
         else
         {
            ssTable[SECSTR_IDX_CHISGN][resCount] = '+';
         }
      }
   }
   
   for(resCount = 0; resCount < seqlen; resCount++)
   {
      sheetCode[resCount]       = 0;
      strandCode[0][resCount]   = 0;
      strandCode[1][resCount]   = 0;
      bridgePoints[0][resCount] = 0;
      bridgePoints[1][resCount] = 0;
   }
   strandNumber = 0;

   for(resCount=0; resCount<seqlen; resCount++)
   {
      for(bridgeIndex=0; bridgeIndex<2; bridgeIndex++)
      {
         if((bridge[bridgeIndex][resCount] != 0) && 
            (abs(bridge[bridgeIndex][resCount]) >= resCount))
         {
            /* bridge[][] array contains the bridge number +1           */
            theBridge = abs(bridge[bridgeIndex][resCount]) - 1; 
            bridgeDirection = bridge[bridgeIndex][resCount] / 
               (theBridge+1);
            bridgePoint = 1;

            if(abs(bridge[0][theBridge]) == resCount+1) 
               bridgePoint = 0;

            for(resCountInSS = (resCount + 2); 
                resCountInSS <= MIN(seqlen,resCount+BULGE_SIZE_L); 
                resCountInSS++)
            {
               for(newBridgeIndex = 0; 
                   newBridgeIndex < NUM_STRAND_PAIR; 
                   newBridgeIndex++)
               {
                  if((bridge[newBridgeIndex][resCountInSS-1] * 
                      bridgeDirection) > 0)
                  {
                     /* bridge[][] array contains the bridge number +1 */
                     theNewBridge = abs(bridge[newBridgeIndex]
                                              [resCountInSS-1]) - 1;

                     if(abs(theNewBridge-theBridge)!=0)
                     {
                        if((theNewBridge-theBridge) /
                           abs(theNewBridge-theBridge) == 
                            bridgeDirection)
                        {
                           if((resCountInSS-resCount-1 >  BULGE_SIZE_S && 
                               abs(theNewBridge-theBridge) <= 
                               BULGE_SIZE_S)  ||
                              (resCountInSS-resCount-1 <= BULGE_SIZE_S && 
                               abs(theNewBridge-theBridge) <= 
                               BULGE_SIZE_L)) 
                           {
                              newBridgePoint = 1;
                              if(abs(bridge[0][theNewBridge]) == 
                                 resCountInSS) 
                              {
                                 newBridgePoint = 0;
                              }
                              
                              for(i = resCount; i< resCountInSS; i++)
                              {
                                 if(ssTable[SECSTR_IDX_SHEET][i] == 
                                    SECSTR_COIL)
                                 {
                                    ssTable[SECSTR_IDX_SHEET][i] = 
                                       SECSTR_SHEET;
                                 }
                              }

                              if(bridge[bridgeIndex+2][resCount] < 0)
                              {
                                 ssTable[SECSTR_IDX_SHEET][resCount] = 
                                    SECSTR_SHEET_SMALL;
                              }
                              
                              if(bridge[newBridgeIndex+2][resCountInSS-1] 
                                 < 0)
                              {
                                 ssTable[SECSTR_IDX_SHEET][resCountInSS-1] =
                                    SECSTR_SHEET_SMALL;
                              }
                              
                              for(i = theBridge; 
                                  (bridgeDirection < 0) ? 
                                     i>=theNewBridge : 
                                     i<=theNewBridge; 
                                  i+=bridgeDirection)
                              {
                                 if(ssTable[SECSTR_IDX_SHEET][i] == 
                                    SECSTR_COIL)
                                 {
                                    ssTable[SECSTR_IDX_SHEET][i] = 
                                       SECSTR_SHEET;
                                 }
                              }

                              if(bridge[bridgePoint+2][theBridge] < 0)
                              {
                                 ssTable[SECSTR_IDX_SHEET][theBridge] = 
                                    SECSTR_SHEET_SMALL;
                              }

                              if(bridge[newBridgePoint+2][theNewBridge] <
                                 0)
                              {
                                 ssTable[SECSTR_IDX_SHEET][theNewBridge] =
                                    SECSTR_SHEET_SMALL;
                              }
                              
                              if(strandCode[bridgeIndex][resCount] == 0)
                              {
                                 strandNumber++;
                                 strandCode[bridgeIndex][resCount] = 
                                    strandNumber * bridgeDirection;
                                 strandCode[bridgePoint][theBridge] = 
                                    strandNumber * bridgeDirection;
                              }

                              strandCode[newBridgeIndex][resCountInSS-1] = 
                                 strandCode[bridgeIndex][resCount];
                              strandCode[newBridgePoint][theNewBridge] = 
                                 strandCode[bridgePoint][theBridge];
                              
                              /* Force break from both newBridgeIndex and 
                                 resCountInSS loops 
                              */
                              resCountInSS = MAX(seqlen,
                                                 resCount+1+BULGE_SIZE_L);
                              break;
                           }
                        }
                     }
                  }
               }
            }
            
            ch = SECSTR_BRIDGE_FWD;
            if(bridgeDirection < 0) 
               ch = SECSTR_BRIDGE_BACKWD;
            ssTable[SECSTR_IDX_BRIDGE][resCount]  = ch;
            ssTable[SECSTR_IDX_BRIDGE][theBridge] = ch;
         }
      }
   }


   resCount = 0;
   
   do{
      while(ssTable[SECSTR_IDX_SHEET][resCount] == SECSTR_COIL)
      {
         resCount++;
         if(resCount >= seqlen) 
            goto RES_COUNT_ERROR;
      }
      
      startOfSheet = resCount;
      resCount = startOfSheet + 1;
      
      while(ssTable[SECSTR_IDX_SHEET][resCount] != SECSTR_COIL)
      {
         resCount++;
         if(resCount >= seqlen)
            break;
      }
      endOfSheet = (resCount - 1);
      lastStrand = 0;
      
      FindNextStrand(strandCode, startOfSheet, endOfSheet, &strandCount, 
                     &newBridgeIndex, lastStrand, &bestValue);
      
      while(bestValue != 0)
      {
         lastStrand = fabs(strandCode[newBridgeIndex][strandCount]);
         charCode   = lastStrand%NUM_STRAND_CHARS;

         if(charCode == 0) 
            charCode = NUM_STRAND_CHARS;

         if(strandCode[newBridgeIndex][strandCount] < 0)
         {
            strandChar = upperCaseLetters[charCode];
         }
         else
         {
            strandChar = lowerCaseLetters[charCode];
         }

         strandCharOffset = SECSTR_IDX_BRIDG1;
         
         if(ssTable[SECSTR_IDX_BRIDG1][strandCount] != SECSTR_COIL) 
            strandCharOffset = SECSTR_IDX_BRIDG2;
         
         if(strandCharOffset != SECSTR_IDX_BRIDG2)
         {
            thisStrandOffset = strandCount+1;

            do
            {
               nextStrandOffset = 
                  FindNextCodeOccurrence(thisStrandOffset, endOfSheet,
                                         strandCode[newBridgeIndex]
                                                    [strandCount],
                                         strandCode);
               if(nextStrandOffset <= 0)
                  break;

               if(ssTable[SECSTR_IDX_BRIDG1][nextStrandOffset-1] != 
                  SECSTR_COIL) 
               {
                  strandCharOffset = SECSTR_IDX_BRIDG2;
               }

               thisStrandOffset = nextStrandOffset;
            } while(strandCharOffset != SECSTR_IDX_BRIDG2);
         }
         
         ssTable[strandCharOffset][strandCount] = strandChar;

         bridgePoints[strandCharOffset-SECSTR_IDX_BRIDG1][strandCount] =
            fabs(bridge[newBridgeIndex][strandCount]);

         thisStrandOffset = strandCount+1;
         
         for(;;)
         {
            nextStrandOffset = 
               FindNextCodeOccurrence(thisStrandOffset,endOfSheet,
                                      strandCode[newBridgeIndex]
                                                [strandCount],
                                      strandCode);

            if(nextStrandOffset <= 0)
               break;

            for(bulgeCount =  thisStrandOffset; 
                bulgeCount < (nextStrandOffset - 1); 
                bulgeCount++)
            {
               ssTable[strandCharOffset][bulgeCount] = SECSTR_BULGE;
            }

            ssTable[strandCharOffset][nextStrandOffset-1] = strandChar;
            strandOffset = 0;

            if(strandCode[0][nextStrandOffset-1] != 
                strandCode[newBridgeIndex][strandCount]) 
            {
               strandOffset = 1;
            }

            bridgePoints[strandCharOffset-SECSTR_IDX_BRIDG1]
                        [nextStrandOffset-1] = 
               fabs(bridge[strandOffset][nextStrandOffset-1]);

            thisStrandOffset = nextStrandOffset;
         }
         
         FindNextStrand(strandCode, startOfSheet, endOfSheet, 
                        &strandCount, &newBridgeIndex, lastStrand, 
                        &bestValue);
      }
      resCount = endOfSheet + 1;
   } while (resCount < seqlen);
   
   
RES_COUNT_ERROR: 
   if(verbose && (lastStrand > NUM_STRAND_CHARS))
   {
      fprintf(stderr,"Sec Struc: (warning) There are %3d strands. Labels \
have restarted\n",
                 lastStrand);
   }
   
   LabelSheets(ssTable, bridgePoints, sheetCode, seqlen, verbose);
   
   for(resCount=0; resCount<seqlen; resCount++)
   {
      if(sheetCode[resCount] != 0)
      {
         charCode                             = sheetCode[resCount];
         ssTable[SECSTR_IDX_SHEET1][resCount] = nstchp[charCode];
      }
   }
   
   bridgeCount = 0;

   for(resCount=0; resCount<seqlen; resCount++)
   {
      for(newBridgeIndex=0; 
          newBridgeIndex<NUM_STRAND_PAIR; 
          newBridgeIndex++)
      {
         if((bridge[newBridgeIndex][resCount]      != 0) && 
            (strandCode[newBridgeIndex][resCount]  == 0) && 
            (fabs(bridge[newBridgeIndex][resCount]) >= resCount))
         {
            bridgeCount++;
            theBridgePoint = fabs(bridge[newBridgeIndex][resCount]);
            charCode       = bridgeCount%NUM_STRAND_CHARS;

            if(charCode == 0) 
               charCode = NUM_STRAND_CHARS;

            if(bridge[newBridgeIndex][resCount] < 0)
            {
               bridgeChar = upperCaseLetters[charCode];
            }
            else
            {
               bridgeChar = lowerCaseLetters[charCode];
            }
            
            bridgeArrayIndex = SECSTR_IDX_BRIDG1;
            
            if(ssTable[SECSTR_IDX_BRIDG1][resCount] != SECSTR_COIL) 
               bridgeArrayIndex = SECSTR_IDX_BRIDG2;

            ssTable[bridgeArrayIndex][resCount] = bridgeChar;

            bridgePoints[bridgeArrayIndex-SECSTR_IDX_BRIDG1][resCount] =
               theBridgePoint;

            bridgeArrayIndex = SECSTR_IDX_BRIDG1;

            if(ssTable[SECSTR_IDX_BRIDG1][theBridgePoint-1] != 
               SECSTR_COIL) 
               bridgeArrayIndex = SECSTR_IDX_BRIDG2;

            ssTable[bridgeArrayIndex][theBridgePoint-1] = bridgeChar;

            bridgePoints[bridgeArrayIndex-SECSTR_IDX_BRIDG1]
                        [theBridgePoint-1] = resCount+1;
         }
      }
   }
   
   if(verbose && (bridgeCount > NUM_STRAND_CHARS))
   {
      fprintf(stderr,"Sec Struc: (warning) There are %3d bridges. Labels \
have restarted\n", 
              bridgeCount);
   }

   if(bridge     != NULL)
      blFreeArray2D((char **)bridge,     NUM_BRIDGE, seqlenMalloc);
   if(strandCode != NULL)
      blFreeArray2D((char **)strandCode, NUM_STRAND_PAIR, seqlenMalloc);
   if(sheetCode  != NULL)
      free(sheetCode);

   return(TRUE);
}


/************************************************************************/
/*>static void RejectHBond(int resIndex, int *bond, REAL *bonde)
   -------------------------------------------------------------
*//**
   \param[in]     resIndex Current position in the list of residues
   \param[in,out] *bond    Bond array
   \param[in,out] *bonde   Bond energy array

   Reject a hydrogen bond because others have better energy

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void RejectHBond(int resIndex, int *bond, REAL *bonde)
{
   int i=1;

   if(bond[0] != resIndex) 
      i++;

   if(i==1)
   {
      bond[0]  = bond[1];
      bonde[0] = bonde[1];
   }

   bond[1]  = 0;
   bonde[1] = 0.0;
}


/************************************************************************/
/*>static void SetSheet(int bridgePoint, int sheetNum, int *sheetCode, 
                        char **ssTable, int seqlen)
   -------------------------------------------------------------------
*//**
   \param[in]  bridgePoint A bridge point residue
   \param[in]  sheetNum    The sheet number
   \param[out] *sheetCode  The sheet number for each residue
   \param[in]  **ssTable   The sec struc table [struc symbol][res]
   \param[in]  seqlen      The sequence length

   Find the bounds of the parent sheet containing bridgePoint and flag
   this as partially set as sheet

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void SetSheet(int bridgePoint, 
                     int sheetNum, 
                     int *sheetCode, 
                     char **ssTable,
                     int seqlen)
{
   int resCount, 
       sheetStart, 
       sheetEnd;
      

   resCount = bridgePoint-1;
   
   while(ssTable[SECSTR_IDX_SHEET][resCount] != SECSTR_COIL)
   {
      resCount++;
      if(resCount >= seqlen) 
         break;
   }
   
   sheetEnd = resCount - 1;
   resCount = bridgePoint-1;
   
   while(ssTable[SECSTR_IDX_SHEET][resCount] != SECSTR_COIL)
   {
      resCount--;
      if(resCount < 1)
         break;
   }
   
   sheetStart = resCount + 1;
   
   for(resCount = sheetStart; resCount <= sheetEnd; resCount++)
   {
      sheetCode[resCount] = (-sheetNum);
   }
}

      
/************************************************************************/
/*>static void LabelSheets(char **ssTable, int **bridgePoints,
                           int *sheetCode, int seqlen, BOOL verbose)
   -----------------------------------------------------------------
*//**
   \param[in]  **ssTable      The sec struc table [struc symbol][res]
   \param[in]  **bridgePoints Bridges between strands
   \param[out] *sheetCode     The sheet number for each residue
   \param[in]  seqlen         The sequence length
   \param[in]  verbose        Print messages

   Label the sheets. At this stage we use numeric labels. These are turned
   into alphabetic symbols later

-  19.05.99 Original   By: ACRM
-  27.05.99 Standard format for messages
-  13.07.15 Modified for BiopLib
*/
static void LabelSheets(char **ssTable, int **bridgePoints,
                        int *sheetCode, int seqlen, BOOL verbose)
{
   int  resCount, 
        sheetStart, 
        sheetEnd, 
        sheetNum, 
        resCountInSheet, 
        bridgeIndex;
   BOOL found;

   for(resCount=0; resCount<seqlen; resCount++)
   {
      sheetCode[resCount] = 0;
   }

   sheetNum = 0;
   resCount = 1;
   
   for(;;)
   {
      FindFirstUnsetSheet(resCount, ssTable, sheetCode, &sheetStart, 
                          &sheetEnd, seqlen);
      if(sheetStart != 0)
      {
         sheetNum++;
         for(resCountInSheet = sheetStart; 
             resCountInSheet <= sheetEnd; 
             resCountInSheet++)
         {
            sheetCode[resCountInSheet-1] = sheetNum;

            for(bridgeIndex = 0; 
                bridgeIndex < NUM_BRIDGE_PAIR; 
                bridgeIndex++)
            {
               if(bridgePoints[bridgeIndex][resCountInSheet-1] != 0)
               {
                  if(sheetCode[bridgePoints[bridgeIndex]
                                           [resCountInSheet-1]-1] == 0)
                  {
                     SetSheet(bridgePoints[bridgeIndex][resCountInSheet-1], 
                              sheetNum, sheetCode, ssTable, seqlen);
                  }
               }
            }
         }
         
         do
         {  
            FindNextPartialSheet(sheetEnd+1, sheetNum, ssTable, 
                                 sheetCode, bridgePoints, &found, seqlen);
         } while(found);
         
         resCount = sheetEnd + 1;

         if(resCount > seqlen)
            break;
      }
      else
      {
         break;
      }
   }
   
   if(verbose && (sheetNum > NUM_STRAND_CHARS))
   {
      fprintf(stderr,"Sec Struc: (warning) There are %3d sheets. Labels \
have restarted\n",
              sheetNum);
   }
}


/************************************************************************/
/*>static int FindNextCodeOccurrence(int strandStart, int strandEnd,
                                     int code, int **strandCode)
   -----------------------------------------------------------------
*//**
   \param[in] strandStart   Start of a strand
   \param[in] strandEnd     End of a strand
   \param[in] code          Current strand code
   \param[in] **strandCode  Array of strand codes[bridge][strand]
   \return                  Offset into strandCode array

   Search for the next occurrence of 'code' in the strandCode array

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static int FindNextCodeOccurrence(int strandStart, int strandEnd,
                                  int code,        int **strandCode)
{
   int nxpl;
   
   for(nxpl=strandStart; nxpl<strandEnd; nxpl++)
   {
      if((strandCode[0][nxpl] == code) ||
         (strandCode[1][nxpl] == code))
      {
         return(nxpl+1);
      }
   }
   
   return(0);
}

/************************************************************************/
/*>static void MarkTurns(char **ssTable, char *detailSS, char *finalSS,
                         char turnChar, char altTurnChar, int seqlen)
   --------------------------------------------------------------------
*//**
   \param[in]  **ssTable   The sec struc table [struc symbol][res]
   \param[out] *detailSS   detailed ss assignment array
   \param[out] *finalSS    final ss assignment array
   \param[in]  turnChar    Character for a turn (T)
   \param[in]  altTurnChar Alternate character for a turn (t)
   \param[in]  seqlen      The sequence length

   Search for isolated turn type bonds. If found set them to the 'turn'
   symbol

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void MarkTurns(char **ssTable, char *detailSS,   char *finalSS,
                      char turnChar,  char altTurnChar, int seqlen)
{
   int  resCount, 
        symbolCount, 
        helixCount, 
        helixPri[NUM_HELIX_TYPE], 
        helixSize[NUM_HELIX_TYPE];

   helixPri[0]  = SECSTR_IDX_ALPHAH;
   helixPri[1]  = SECSTR_IDX_THRTNH;
   helixPri[2]  = SECSTR_IDX_PIH;
   helixSize[0] = HELIX_CHUNK_SIZE;
   helixSize[1] = THREE_TEN_CHUNK_SIZE;
   helixSize[2] = PIH_CHUNK_SIZE;
   
   for(helixCount = 0; helixCount < NUM_HELIX_TYPE; helixCount++)
   {
      if(ssTable[helixPri[helixCount]][0] == SECSTR_BEND_START &&
         ssTable[helixPri[helixCount]][1] != SECSTR_BEND_START) 
      {
         for(symbolCount = 1; 
             symbolCount < helixSize[helixCount]; 
             symbolCount++)
         {
            if(detailSS[symbolCount] == SECSTR_COIL) 
            {
               detailSS[symbolCount] = turnChar;
            }

            if(finalSS[symbolCount] == SECSTR_COIL || 
               finalSS[symbolCount] == altTurnChar) 
            {
               finalSS[symbolCount] = turnChar;
            }
         }
         
         if(finalSS[0] == SECSTR_COIL) 
            finalSS[0] = altTurnChar;

         if(finalSS[helixSize[helixCount]] == SECSTR_COIL)
            finalSS[helixSize[helixCount]] = altTurnChar;
      }
      
      for(resCount = 1; 
          resCount < seqlen - helixSize[helixCount]; 
          resCount++)
      {
         if((ssTable[helixPri[helixCount]][resCount] == 
             SECSTR_BEND_START ||
             ssTable[helixPri[helixCount]][resCount] == 
             SECSTR_BEND_BOTH) &&
            (ssTable[helixPri[helixCount]][resCount-1] != 
             SECSTR_BEND_START &&
             ssTable[helixPri[helixCount]][resCount-1] != 
             SECSTR_BEND_BOTH) &&
            (ssTable[helixPri[helixCount]][resCount+1] != 
             SECSTR_BEND_START &&
             ssTable[helixPri[helixCount]][resCount+1] != 
             SECSTR_BEND_BOTH)) 
         {
            for(symbolCount = resCount + 1; 
                symbolCount < (resCount + helixSize[helixCount]); 
                symbolCount++)
            {
               if(detailSS[symbolCount] == SECSTR_COIL) 
               {
                  detailSS[symbolCount] = turnChar;
               }

               if(finalSS[symbolCount] == SECSTR_COIL || 
                  finalSS[symbolCount] == altTurnChar) 
               {
                  finalSS[symbolCount] = turnChar;
               }
            }

            if(finalSS[resCount] == SECSTR_COIL) 
               finalSS[resCount] = altTurnChar;

            if(finalSS[resCount+helixSize[helixCount]] == SECSTR_COIL)
               finalSS[resCount+helixSize[helixCount]] = altTurnChar;
         }
      }
   }
}


/************************************************************************/
/*>static void MarkSheetsAndBridges(int seqlen, char **ssTable, 
                                    char *detailSS, char *finalSS, 
                                    char *ssSymbols, char *altSSSymbols)
   ---------------------------------------------------------------------
*//**
   \param seqlen         The sequence length
   \param **ssTable      The sec struc table [struc symbol][res]
   \param *detailSS      detailed ss assignment array
   \param *finalSS       final ss assignment array
   \param *ssSymbols     array of symbols for secondary structures
   \param *altSSSymbols  array of alternative symbols for sec strucs

   Put in symbols for the secondary structure elements for sheets and
   bridges.

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void MarkSheetsAndBridges(int seqlen, char **ssTable, 
                                 char *detailSS, char *finalSS, 
                                 char *ssSymbols, char *altSSSymbols)
{
   int resCount;
   
   for(resCount = 0; resCount < seqlen; resCount++)
   {
      if(ssTable[SECSTR_IDX_SHEET][resCount] != SECSTR_COIL)
      {
         if(detailSS[resCount] == SECSTR_COIL) 
            detailSS[resCount] = ssSymbols[SECSTR_IDX_SHEET];

         if(finalSS[resCount] == SECSTR_COIL || 
            finalSS[resCount] == altSSSymbols[SECSTR_IDX_ALPHAH] || 
            finalSS[resCount] == altSSSymbols[SECSTR_IDX_SHEET])
         {
            finalSS[resCount] = ssSymbols[SECSTR_IDX_SHEET];
         }
         
         if((ssTable[SECSTR_IDX_SHEET][resCount] == SECSTR_SHEET_SMALL) &&
            (resCount > 0))
         {
            if(finalSS[resCount-1] == SECSTR_COIL || 
               finalSS[resCount-1] == ssSymbols[SECSTR_IDX_BRIDGE]) 
            {
               finalSS[resCount-1] = altSSSymbols[SECSTR_IDX_SHEET];
            }

            if(resCount < seqlen-1)
            {
               if(finalSS[resCount+1] == SECSTR_COIL || 
                  finalSS[resCount-1] == ssSymbols[SECSTR_IDX_BRIDGE]) 
               {
                  finalSS[resCount+1] = altSSSymbols[SECSTR_IDX_SHEET];
               }
            }
         }
      }

      if(ssTable[SECSTR_IDX_BRIDGE][resCount] != SECSTR_COIL)
      {
         if(detailSS[resCount] == SECSTR_COIL) 
            detailSS[resCount] = ssSymbols[SECSTR_IDX_BRIDGE];

         if(finalSS[resCount] == SECSTR_COIL) 
            finalSS[resCount] = ssSymbols[SECSTR_IDX_BRIDGE];
      }
   }
}


/************************************************************************/
/*>static void FindNextStrand(int **strandCode, int sheetStart, 
                              int sheetEnd, int *strandCount, 
                              int *startIndex, int lastStrand,
                              int *bestValue)
   ------------------------------------------------------------
*//**
   \param[in]  **strandCode Array of strand codes[bridge][strand]
   \param[in]  sheetStart   Start of a sheet
   \param[in]  sheetEnd     End of a sheet
   \param[out] *strandCount The strand count
   \param[out] *startIndex  The new bridge 
   \param[in]  lastStrand   The last strand
   \param[out] *bestValue   Best strand

   Find where the next beta strand starts

-  19.05.99 Original   By: ACRM
-  13.07.15 Modified for BiopLib
*/
static void FindNextStrand(int **strandCode, int sheetStart, int sheetEnd,
                           int *strandCount, int *startIndex, 
                           int lastStrand, int *bestValue)
{
   int  nxstr, nxpl;
   
   *bestValue = 0;
   
   for(nxstr=sheetStart; nxstr<=sheetEnd; nxstr++)
   {
      for(nxpl=0; nxpl<NUM_STRAND_PAIR; nxpl++)
      {
         if(abs(strandCode[nxpl][nxstr]) > lastStrand)
         {
            if(*bestValue == 0)
            {
               *bestValue = abs(strandCode[nxpl][nxstr]);
               *startIndex = nxpl;
               *strandCount = nxstr;
            }
            else
            {
               if(abs(strandCode[nxpl][nxstr]) < *bestValue)
               {
                  *bestValue = abs(strandCode[nxpl][nxstr]);
                  *startIndex = nxpl;
                  *strandCount = nxstr;
               }
            }
         }
      }
   }
}



