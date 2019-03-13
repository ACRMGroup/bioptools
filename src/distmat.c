/*************************************************************************

   Program:    distmat
   File:       distmat.c
   
   Version:    V2.1
   Date:       13.03.19
   Function:   Calculate inter-CA distances on a set of common-labelled
               PDB files
   
   Copyright:  (c) UCL, Dr. Andrew C. R. Martin 2009-2019
   Author:     Dr. Andrew C. R. Martin
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

   This is a simple program for calculating means and standard
   deviations for inter residue distances from a set of PDB files
   having common numbering (e.g. antibodies). Input is either a single
   PDB file or a 'file of files' - i.e. a file containing a list of
   PDB files to be processed.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0   01.04.09  Original
   V1.1   06.04.09  Added -n and -m options
   V1.2   30.11.16  Minor cleanup for bioptools - Added option to take
                    a single PDB file as input rather than a set (-p)
   V2.0   01.12.16  Rewrite to use hashes. Removed -n and -m since these
                    are't needed any more - everything is dynamically
                    allocated.
   V2.1   13.03.19  Increased some buffer sizes

*************************************************************************/
/* #define DEBUG 1 */

/************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "bioplib/MathUtil.h"
#include "bioplib/hash.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXCHAINLABEL 8
#define MAXBUFF       160
#define DEF_MAXRES    300   /* Approximate number of residues in a file.
                               This is simply used to initialize the hash
                               table. Increasing it will increase 
                               efficiency for large structures, but waste
                               memory for small ones
                            */
#define MAXLABEL      16    /* ResidueLabel size                        */
#define ATOMS_CA      0     /* Selection types                          */
#define ATOMS_ALL     1
#define ATOMS_SC      2

typedef struct respair
{
   REAL sx;
   REAL sxsq;
   int  nval;
}  RESPAIR;
   

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *singleFile, int *atomTypes, char *chains);
void HandleInput(FILE *in, FILE *out, BOOL singleFile, 
                 HASHTABLE *hashTable, int atomTypes, char *chains);
void ProcessFile(FILE *fp, HASHTABLE *hashTable, int atomTypes,
                 char **chainList);
void ProcessPDB(PDB *pdb, HASHTABLE *hashTable);
void StoreData(HASHTABLE *hashTable, PDB *res1, PDB *res2, REAL dist);
PDB *ReduceAtomList(PDB *pdb, int atomTypes);
void DisplayResults(FILE *out, HASHTABLE *hashTable);
PDB *FindEndOfChain(PDB *chain);
BOOL ValidChain(PDB *pdb, char **chains);
PDB *SelectPDBChains(PDB *pdb, char **chains);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program for distance analysis

-  01.04.09 Original   By: ACRM
-  06.04.09 Added -n and -m parameters
*/
int main(int argc, char **argv)
{
   char  infile[MAXBUFF],
         outfile[MAXBUFF];
   FILE  *in  = stdin,
         *out = stdout;
   BOOL  singleFile  = FALSE;
   ULONG hashSize    = DEF_MAXRES * DEF_MAXRES;
   int   atomTypes   = ATOMS_CA;
   HASHTABLE *hashTable = NULL;
   char      chains[MAXBUFF];

   chains[0] = '\0';

   if(ParseCmdLine(argc, argv, infile, outfile, &singleFile, &atomTypes,
                   chains))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((hashTable = blInitializeHash(hashSize))!=NULL)
         {
            HandleInput(in, out, singleFile, hashTable, atomTypes, chains);
            DisplayResults(out, hashTable);
         }
         else
         {
            fprintf(stderr,"ERROR: Unable to initialize hash table.\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"ERROR: Unable to open files.\n");
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
                     BOOL *singleFile, int *atomTypes, char *chains)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            BOOL   *singleFile  Input is a single PDB file instead of
                                a list
            char   *chains      Comma-separated list of chains to keep
   Returns: BOOL                Success?

   Parse the command line

-  01.04.09 Original    By: ACRM
-  06.04.09 Added -n and -m and their parameters
-  30.11.16 Added -p
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *singleFile, int *atomTypes, char *chains)
{
   argc--;
   argv++;

   infile[0]    = outfile[0] = '\0';
   *singleFile  = FALSE;
   
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'p':
            *singleFile = TRUE;
            break;
         case 'c':
            argc--;
            argv++;
            if(!argc) return(FALSE);
            strncpy(chains, argv[0], MAXBUFF);
            break;
         case 'a':
            *atomTypes = ATOMS_ALL;
            break;
         case 's':
            *atomTypes = ATOMS_SC;
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
/*>void HandleInput(FILE *in, FILE *out, BOOL singleFile, 
                    HASHTABLE *hashTable, int atomTypes, char *chains)
   -------------------------------------------------------------------
*//**
   \input[in]     in          Input file pointer
   \input[in]     out         Output file pointer
   \input[in]     singleFile  Input is a PDB file rather than a file of
                              files
   \input[in,out] hashTable   Analysis data for each residue pair
   \input[in]     atomTypes   Atom types to include
   \input[in]     chains      Comma separated list of chain names 
                              (or blank)

   Handle the input file - extract the PDB filenames and process each 
   in turn, or just the one file if singleFile is set.

-  01.04.09 Original   By: ACRM
-  06.04.09 Handles maxchain
-  30.11.16 Added singleFile
-  01.12.16 Major rewrite
*/
void HandleInput(FILE *in, FILE *out, BOOL singleFile, 
                 HASHTABLE *hashTable, int atomTypes, char *chains)
{
   char filename[MAXBUFF];
   char **chainList = NULL;
   
   if(chains[0])
   {
      if((chainList = blSplitStringOnCommas(chains, MAXCHAINLABEL))==NULL)
      {
         fprintf(stderr,"No memory for storing chain labels: %s\n",
                 chains);
         exit(1);
      }
   }

   if(singleFile)
   {
      ProcessFile(in, hashTable, atomTypes, chainList);
   }
   else
   {
      while(fgets(filename,MAXBUFF,in))
      {
         FILE *fp;
         
         TERMINATE(filename);
         
         if((fp = fopen(filename,"r"))!=NULL)
         {
            fprintf(stderr,"INFO: Processing file: %s\n",filename);
            ProcessFile(fp, hashTable, atomTypes, chainList);
            fclose(fp);
         }
         else
         {
            fprintf(stderr,"WARNING: Unable to read file: %s\n", 
                    filename);
         }
      }
   }
}


/************************************************************************/
/*>void ProcessFile(FILE *fp, HASHTABLE *hashTable, int atomTypes, 
                    char **chainList)
   ---------------------------------------------------------------
*//**
   \param[in]      fp          File pointer for input file
   \input[in,out]  hashTable   Analysis data for each residue pair
   \param[in]      atomTypes   Atom types to keep
   \param[in]      chainList   List of chains to keep (Keep all if NULL)

   Processes an individual PDB file, selecting required atoms and
   chains if necessary 

-  01.12.16 Original - Complete new version   By: ACRM  
*/
void ProcessFile(FILE *fp, HASHTABLE *hashTable, int atomTypes, 
                 char **chainList)
{
   PDB *pdb;
   int natoms;
   
   if((pdb = blReadPDBAtoms(fp, &natoms))!=NULL)
   {
      pdb = ReduceAtomList(pdb, atomTypes);
      if(chainList)
      {
         pdb = SelectPDBChains(pdb, chainList);
      }
      if(pdb!=NULL)
         ProcessPDB(pdb, hashTable);
      FREELIST(pdb, PDB);
   }
   else
   {
      fprintf(stderr,"Error: No memory to read PDB file\n");
      exit(1);
   }
}


/************************************************************************/
/*>void ProcessPDB(PDB *pdb, HASHTABLE *hashTable)
   -----------------------------------------------
*//**
   \input[in]      pdb        PDB linked list
   \input[in,out]  hashTable  Analysis data for each residue pair

   Does the actual analysis of a PDB linked list

-  01.12.16 Original - Complete new version   By: ACRM  
*/
void ProcessPDB(PDB *pdb, HASHTABLE *hashTable)
{
   PDB *atom1, *atom2, 
       *res1,  *res1Next, 
       *res2,  *res2Next;
   REAL minDistSq = (REAL)0.0;

   /* Step through a residue at a time                                  */
   for(res1=pdb; res1!=NULL; res1=res1Next)
   {
      /* Find the start of the next residue                             */
      res1Next = blFindNextResidue(res1);
      
      /* Step through again a residue at a time                         */
      for(res2=pdb; res2!=NULL; res2=res2Next)
      {
         /* Find the start of the next residue                          */
         res2Next = blFindNextResidue(res2);

         /* Initialize minimum distance between the residues            */
         minDistSq = DISTSQ(res1, res2);

         /* Step through atoms in first residue                         */
         for(atom1=res1;
             ((atom1!=NULL)&&(atom1!=res1Next));
             NEXT(atom1))
         {
            /* Step through atoms in second residue to find 
               the minimum distance between the two residues
            */
            for(atom2=res2; 
                ((atom2!=NULL)&&(atom2!=res2Next));
                 NEXT(atom2))
            {
               REAL dSq = DISTSQ(atom1, atom2);
               if(dSq < minDistSq)
               {
                  minDistSq = dSq;
               }
            }
         }

         StoreData(hashTable, res1, res2, sqrt(minDistSq));
      }
   }
}


/************************************************************************/
/*>void StoreData(HASHTABLE *hashTable, PDB *res1, PDB *res2, REAL dist)
   ---------------------------------------------------------------------
*//**
   \input[in,out]  hashTable   Analysis data for each residue pair
   \input[in]      res1        First residue PDB pointer
   \input[in]      res2        Second residue PDB pointer
   \input[in]      dist        Distance between residues

   Creates a hash key from the two residues, allocates memory for data
   for this residue pair if needed and store the data in the hash keyed
   by the residue pair.

-  01.12.16 Original - Complete new version   By: ACRM  
-  13.03.19 Added 1 to resID1, resID2 and resPair sizes  
*/
void StoreData(HASHTABLE *hashTable, PDB *res1, PDB *res2, REAL dist)
{
   RESPAIR *rp = NULL;
   char resID1[MAXLABEL+1], resID2[MAXLABEL+1], resPair[(MAXLABEL+1)*2];
   
   MAKERESID(resID1, res1);
   MAKERESID(resID2, res2);
   sprintf(resPair, "%s-%s", resID1, resID2);
   
   /* This key not in the hash so create it */
   if(!blHashKeyDefined(hashTable, resPair))
   {
      if((rp = (RESPAIR *)malloc(sizeof(RESPAIR))) == NULL)
      {
         fprintf(stderr, "Error: no memory for RESPAIR data\n");
         exit(1);
      }

      rp->sx   = 0.0;
      rp->sxsq = 0.0;
      rp->nval = 0;

      blSetHashValuePointer(hashTable, resPair, (BPTR)rp);
   }
   else
   {
      if((rp=(RESPAIR *)blGetHashValuePointer(hashTable, resPair))==NULL)
      {
         fprintf(stderr, "Error: internal Hash confused!\n");
         exit(1);
      }
   }
   
   blCalcExtSD(dist, 0, &(rp->sx), &(rp->sxsq), &(rp->nval), NULL, NULL);
}


/************************************************************************/
/*>PDB *ReduceAtomList(PDB *pdb, int atomTypes)
   --------------------------------------------
*//**
   \input[in]   pdb        PDB linked list
   \input[in]   atomTypes  Atom types to keep
   \return                 New PDB linked list

   Removes atoms not being analyzed from the linked list. atomTypes is 
   one of ATOMS_ALL, ATOMS_CA or ATOMS_SC

-  01.12.16 Original - Complete new version   By: ACRM  
*/
PDB *ReduceAtomList(PDB *pdb, int atomTypes)
{
   PDB  *reduced;
   char *sel[40];
   int  natoms;
   
   switch(atomTypes)
   {
   case ATOMS_ALL:
      return(pdb);
      break;
   case ATOMS_CA:
      SELECT(sel[0],"CA  ");
      if(sel[0] == NULL) return(NULL);
      reduced = blSelectAtomsPDBAsCopy(pdb, 1, sel, &natoms);
      FREELIST(pdb, PDB);
      return(reduced);
      break;
   case ATOMS_SC:
      SELECT(sel[0], "CB  ");
      SELECT(sel[1], "CD  ");
      SELECT(sel[2], "CD1 ");
      SELECT(sel[3], "CD2 ");
      SELECT(sel[4], "CE  ");
      SELECT(sel[5], "CE1 ");
      SELECT(sel[6], "CE2 ");
      SELECT(sel[7], "CE3 ");
      SELECT(sel[8], "CG  ");
      SELECT(sel[9], "CG1 ");
      SELECT(sel[10],"CG2 ");
      SELECT(sel[11],"CH2 ");
      SELECT(sel[12],"CZ  ");
      SELECT(sel[13],"CZ2 ");
      SELECT(sel[14],"CZ3 ");
      SELECT(sel[15],"ND1 ");
      SELECT(sel[16],"ND2 ");
      SELECT(sel[17],"NE  ");
      SELECT(sel[18],"NE1 ");
      SELECT(sel[19],"NE2 ");
      SELECT(sel[20],"NH1 ");
      SELECT(sel[21],"NH2 ");
      SELECT(sel[22],"NZ  ");
      SELECT(sel[23],"OD1 ");
      SELECT(sel[24],"OD2 ");
      SELECT(sel[25],"OE1 ");
      SELECT(sel[26],"OE2 ");
      SELECT(sel[27],"OG  ");
      SELECT(sel[28],"OG1 ");
      SELECT(sel[29],"OH  ");
      SELECT(sel[30],"SD  ");
      SELECT(sel[31],"SG  ");
      if(sel[31] == NULL) return(NULL);
      reduced = blSelectAtomsPDBAsCopy(pdb, 32, sel, &natoms);
      FREELIST(pdb, PDB);
      return(reduced);
      break;
   }
   return(pdb);
}


/************************************************************************/
/*>void DisplayResults(FILE *out, HASHTABLE *hashTable)
   -----------------------------------------------------
   \input[in]      out         Output file pointer
   \input[in,out]  hashTable   Analysis data for each residue pair

   Display the results. Run through each indexed location and calculate 
   the mean and standard deviation then print the residue IDs with these
   values.

-  01.04.09 Original   By: ACRM
-  01.12.16 Major rewrite
-  13.03.19 Added 1 to res1 and res2 sizes and terminate string
*/
void DisplayResults(FILE *out, HASHTABLE *hashTable)
{
   char **keys = NULL;
   int  i;
   REAL mean, 
        sd;


   if((keys = blGetHashKeyList(hashTable))!=NULL)
   {
      for(i=0; keys[i] != NULL; i++)
      {
         RESPAIR *rp;
         char    res1[MAXLABEL+1],
                 res2[MAXLABEL+1],
                 *chp;
      
         if((rp=(RESPAIR *)blGetHashValuePointer(hashTable,keys[i]))==NULL)
         {
            fprintf(stderr, "Error: internal Hash confused!\n");
            exit(1);
         }

         strncpy(res1, keys[i], MAXLABEL);
         TERMAT(res1, '-');
         chp = strchr(keys[i], '-');
         strncpy(res2, chp+1, MAXLABEL);
         res2[MAXLABEL]='\0';
      
         blCalcExtSD((REAL)0.0, 1, 
                     &(rp->sx), &(rp->sxsq), &(rp->nval), &mean, &sd);
         fprintf(out,"%s %s %6.3f %6.3f\n", 
                 res1, res2, mean, sd);
      }
      blFreeHashKeyList(keys);
   }            
}



/************************************************************************/
/*>PDB *SelectPDBChains(PDB *pdb, char **chains)
   ---------------------------------------------
*//**
   \input[in]   pdb      PDB linked list
   \input[in]   chain    Array of chains to keep
   \return               Reduced PDB linked list

   Selects the specified chains from the PDB linked list

-  01.12.16 Original - based on code in pdbgetchains   By: ACRM  
*/
PDB *SelectPDBChains(PDB *pdb, char **chains)
{
   PDB  *chainStart    = NULL,
        *endOfChain    = NULL,
        *nextChain     = NULL,
        *keptChainsEnd = NULL;
   
   if(chains == NULL)
      return(pdb);

   for(chainStart=pdb; chainStart!=NULL; chainStart=nextChain)
   {
      endOfChain       = FindEndOfChain(chainStart);
      nextChain        = endOfChain->next;
      endOfChain->next = NULL;

      if(ValidChain(chainStart, chains))
      {
         if(keptChainsEnd!=NULL)
            keptChainsEnd->next = chainStart;
         keptChainsEnd = endOfChain;
      }
      else
      {
         if(chainStart == pdb)
            pdb = nextChain;

         FREELIST(chainStart, PDB);
      }
   }

   return(pdb);
}

/************************************************************************/
/*>BOOL ValidChain(PDB *pdb, char **chains)
   ----------------------------------------
*//**
   \input[in]   pdb     PDB item
   \input[in]   chains  Array of chain labels
   \return              Is it a valid chain?

   Tests whether the chain specified in pdb appears in the list of chains

-  01.12.16 Original - based on code in pdbgetchains   By: ACRM  
*/
BOOL ValidChain(PDB *pdb, char **chains)
{
   int i;

   for(i=0; ((chains[i] != NULL) && (chains[i][0] != '\0')); i++)
   {
      if(CHAINMATCH(pdb->chain, chains[i]))
         return(TRUE);
   }
   
   return(FALSE);
}

/************************************************************************/
/*>PDB *FindEndOfChain(PDB *chain)
   -------------------------------
*//**
   \input[in]   chain    PDB linked list
   \return               Pointer to the last item in the current chain

   Finds the last atom in the current chain (not the start of the next
   chain!)

-  01.12.16 Original - based on code in pdbgetchains   By: ACRM  
*/
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

-  01.04.09 Original   By: ACRM
-  06.04.09 V1.1
-  30.11.16 V1.2
-  01.12.16 V2.0
-  13.03.19 V2.1
*/
void Usage(void)
{
   fprintf(stderr,"\nDistMat V2.1 (c) 2009-2019, Dr. Andrew C.R. Martin, \
UCL\n");

   fprintf(stderr,"\nUsage: distmat [-p][-c chains][-a | -s] [input [output]]\n");
   fprintf(stderr,"       -p Input is a single PDB file instead \
of a file of files\n");
   fprintf(stderr,"       -c chains Only look at specified chaind\n");
   fprintf(stderr,"       -a Look at all atoms rather than CAs\n");
   fprintf(stderr,"       -s Look at sidechain atoms rather than CAs\n");
   fprintf(stderr,"\nI/O Through stdin/stdout if not specified\n");

   fprintf(stderr,"\nDistMat analyses inter-CA distances in one or \
more PDB files\n");

   fprintf(stderr,"\nThe default input file simply contains a list of \
the PDB files to be processed.\n");
   fprintf(stderr,"If -p is specified the input is a single PDB \
file.\n");

   fprintf(stderr,"\nIf -c is specified it is followed by a comma-separated \
list if chain names to analyze.\n\n");
}

