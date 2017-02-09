/************************************************************************/
/**

   \file       pdbsymm.c
   
   \version    V1.0
   \date       06.02.17
   \brief      Program to apply non-crystollographic symmetry operations
               from the REMARK 350 data in a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin / UCL 2017
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
   Note 1: this code currently ends up with PDB files that are not fully
   valid:
   - the headers don't reflect the additional chains
   - there is no MASTER or CONECT record
   - HETATMs are not all moved to the end of the file

   Note 2: this code only works with single character chain names. It 
   needs updating to deal with multi-character names!

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  06.02.17 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/matrix.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/fsscanf.h"
#include "bioplib/general.h"

#define MAXCHAINS 62

/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
void WriteSymmetryCopies(FILE *out, WHOLEPDB *wpdb);
int ReadSymmetryData(WHOLEPDB *wpdb, REAL matrix[3][3], REAL trans[3], 
                     char chains[MAXCHAINS][8]);
BOOL IsIdentityMatrix(REAL matrix[3][3], REAL trans[3]);
char GetNextChainLabel(char chainLabel);
void ApplyMatrixAndWriteCopy(FILE *out, PDB *pdb, char oldChain[8], 
                             char newChain[8],
                             REAL matrix[3][3], REAL trans[3]);
char FindLastChainLabel(PDB *pdb);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for applying non-crystallographic symmertry operators

-  09.02.17 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE     *in      = stdin,
            *out     = stdout;
   WHOLEPDB *wpdb;

   argc--;
   argv++;

   /* Handle all switches                                               */
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            Usage();
            return(0);
         default:
            Usage();
            return(1);
         }
      }
      else
      {
         break;
      }
      
      argc--;
      argv++;
   }

   /* Open input file if specified                                      */
   if(argc)
   {
      if((in=fopen(argv[0],"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input file: %s\n",argv[0]);
         return(1);
      }
      argc--;
      argv++;
   }

   /* Open output file if specified                                     */
   if(argc)
   {
      if((out=fopen(argv[0],"w"))==NULL)
      {
         fprintf(stderr,"pdbrotate: Unable to open output file: %s\n",
                 argv[0]);
         return(1);
      }
      argc--;
      argv++;
   }
   
   /* Read in the PDB file                                              */
   if((wpdb = blReadWholePDB(in))==NULL)
   {
      fprintf(stderr,"pdbrotate: Unable to read from PDB file\n");
      return(1);
   }

   /* Write the initial copy (unmodified)                               */
   blWriteWholePDBHeader(out, wpdb);
   blWritePDB(out, wpdb->pdb);
   
   WriteSymmetryCopies(out, wpdb);

   /* Need to build everything into one PDB linked list and rebuild conect
      data to do this properly
   */
   /* blWriteWholePDBTrailer(out, wpdb);                                */ 

   return(0);
}


/************************************************************************/
/*>char FindLastChainLabel(PDB *pdb)
   ---------------------------------
*//**
   \param[in]    *pdb   PDB linked list
   \return              The last chain label used

   Identifies the 'alphabetically' last chain used given that the order
   used by the PDB is capital letters, digits 1...0 and lower case
   letters.

-  09.02.17  Original   By: ACRM
*/
char FindLastChainLabel(PDB *pdb)
{
   char *permittedChains =
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890abcdefghijklmnopqrstuvwxyz";
   int  lastIdx          = -1;
   PDB  *p;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      int  idx;
      idx = blChindex(permittedChains, p->chain[0]);
      if(idx > lastIdx)
         lastIdx = idx;
   }

   return(permittedChains[lastIdx]);
}


/************************************************************************/
/*>void WriteSymmetryCopies(FILE *out, WHOLEPDB *wpdb)
   ---------------------------------------------------
*//**
   \param[in]    *out   Output file pointer
   \param[in]    *wpdb  Pointer to WHOLEPDB structure

   Does the work of writing non-crystallographic symmetry related copies
   of the structure.

-  09.02.17  Original   By: ACRM
*/
void WriteSymmetryCopies(FILE *out, WHOLEPDB *wpdb)
{
   REAL  matrix[3][3];
   char  chains[MAXCHAINS][8];
   int   nchains;
   REAL  trans[3];
   char  lastChainLabel;
   char  chainLabel;
   PDB   *pdb;

   pdb=wpdb->pdb;

   lastChainLabel = FindLastChainLabel(pdb);
   chainLabel     = GetNextChainLabel(lastChainLabel);

   while((nchains=ReadSymmetryData(wpdb, matrix, trans, chains))!=0)
   {
      if(!IsIdentityMatrix(matrix, trans))
      {
         int chainNum;
         for(chainNum=0; chainNum<nchains; chainNum++)
         {
            char chainLabelString[8];
            chainLabelString[0] = chainLabel;
            chainLabelString[1] = '\0';

            ApplyMatrixAndWriteCopy(out, pdb, 
                                    chains[chainNum], /* Old label      */
                                    chainLabelString, /* New label      */
                                    matrix, trans);
            chainLabel = GetNextChainLabel(chainLabel);
         }
      }
   }
}


/************************************************************************/
/*>int ReadSymmetryData(WHOLEPDB *wpdb, REAL matrix[3][3], REAL trans[3], 
                        char chains[MAXCHAINS][8])
   -----------------------------------------------------------------------
*//**
   \param[in]    *wpdb      Pointer to WHOLEPDB structure
   \param[out]   matrix[][] 3x3 rotation matrix
   \param[out]   trans[]    translation vector
   \param[out]   chains[][] Chain labels to which the symmetry operators
                            are applied  
   \return                  The number of chains to which the operators
                            are applied (0 = no more NC symmetry records)

   Reads the symmetry data from
   REMARK 350 APPLY THE FOLLOWING TO CHAINS:
   records which give a comma separated list of chains to which the
   operations must be applied and from
   REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000            
   12345612341231234511234123456789012345678901234567890123451234567890
   %6s   %4d %3x%5s%1d%4d %10.6f    %10.6f    %10.6f    %5x  %10.6f
   records which specify the rotatation matrix and translation vector.

   A static variable is used to keep track of where we are in the header
   data. Successive calls will give the next set of operations. The 
   routine returns the number of chains affected or zero when there are
   no more NC symmetry records.
 
-  09.02.17  Original   By: ACRM
*/
int ReadSymmetryData(WHOLEPDB *wpdb, REAL matrix[3][3], REAL trans[3], 
                     char chains[MAXCHAINS][8])
{
   static STRINGLIST *sSymOp = NULL;
   static BOOL       called  = FALSE;
   static char       sChains[MAXCHAINS][8];
   static int        sNChains = 0;
   int               i;

   /* First call - step to the start of the symmertry operators         */
   if((sSymOp == NULL) && !called)
   {
      called = TRUE;

      sSymOp = wpdb->header;
   }
   
   /* Not first call and we have run out of records                     */
   if(sSymOp == NULL)
      return(0);

   /* Step on until the next 'APPLY' or the next 'BIOMT1'               */
   while((sSymOp != NULL) &&
         (sSymOp->string != NULL) &&
         strncmp(sSymOp->string, 
                 "REMARK 350 APPLY THE FOLLOWING TO CHAINS:", 41) &&
         strncmp(sSymOp->string, "REMARK 350   BIOMT1", 19))
   {
      NEXT(sSymOp);
   }
   
   /* The list has ended                                                */
   if((sSymOp == NULL) || (sSymOp->string == NULL))
      return(0);

   /* We are now pointing to the first symmetry operator                */
   
   if(!strncmp(sSymOp->string, 
               "REMARK 350 APPLY THE FOLLOWING TO CHAINS:", 41))
   {
      char *chp;
      
      /* Populate the static list of chains                             */
      chp = sSymOp->string + 41;
      TERMINATE(chp);
      KILLTRAILSPACES(chp);
      KILLLEADSPACES(chp, chp);
      for(sNChains=0; sNChains<MAXCHAINS; sNChains++)
      { 
         if(chp == NULL)
            break;
         
         chp = blGetWord(chp, sChains[sNChains], 8);
      }
      NEXT(sSymOp);
   }

   for(i=0; i<sNChains; i++)
   {
      strcpy(chains[i], sChains[i]);
   }
   
   if((sSymOp != NULL) && (sSymOp->string != NULL))
   {
      /* This should not happen!                                        */
      if(strncmp(sSymOp->string, "REMARK 350   BIOMT1", 19))
         return(0);

      /* Read the three REMARK 350 BIOMTx lines                         */
      for(i=0; i<3; i++)
      {
         if((sSymOp != NULL) && (sSymOp->string != NULL))
         {
            char recordType[8],
                 subType[8];
            int  recordNum,
                 subTypeNum,
                 instanceNum;
            REAL x,y,z,t;

            fsscanf(sSymOp->string, "%6s%4d%3x%5s%1d%4d%10lf%10lf%10lf%5x%10lf",
                    recordType, &recordNum, subType, &subTypeNum,
                    &instanceNum, &x, &y, &z, &t);
            if(!strncmp(recordType, "REMARK", 6) &&
               (recordNum == 350) &&
               !strncmp(subType, "BIOMT", 5))
            {
               matrix[i][0] = x;
               matrix[i][1] = y;
               matrix[i][2] = z;
               trans[i]     = t;
            }
            else
            {
               fprintf(stderr,"Warning (pdbsymm): Unexpected record - %s\n",
                       sSymOp->string);
            }
         }
         NEXT(sSymOp);
      }
   }

   return(sNChains);
}


/************************************************************************/
/*>BOOL IsIdentityMatrix(REAL matrix[3][3], REAL trans[3])
   -------------------------------------------------------
*//**
   \param[in]    matrix[][]   A 3x3 rotation matrix
   \param[in]    trans[3]     A translation vector
   \return                    True if it's an identity matrix with a
                              0,0,0 translation vector; False otherwise

   Determines whether this is an identity matrix with no translation -
   i.e. the matrix and vector would have no effect on the structure.

-  09.02.17  Original   By: ACRM
*/
BOOL IsIdentityMatrix(REAL matrix[3][3], REAL trans[3])
{
   int i, j;
   
   /* Test the matrix                                                   */
   for(i=0; i<3; i++)
   {
      for(j=0; j<3; j++)
      {
         if((i==j) && (matrix[i][j] != (REAL)1.0))  
         {  /* Diagonal not 1.0                                         */
            return(FALSE);
         }
         else if((i!=j) && (matrix[i][j] != (REAL)0.0))
         {  /* Off-diagonal not 0.0                                     */
            return(FALSE);
         }
      }
   }

   /* Test the translation vector                                       */
   for(i=0; i<3; i++)
   {
      if(trans[i] != (REAL)0.0)
      {
         return(FALSE);
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>char GetNextChainLabel(char chainLabel)
   ---------------------------------------
*//**
   \param[in]    chainlabel   A chain label (single character!)
   \return                    A new chain label (single character!)

   Bumps the chain label using the standard PDB order of capitals,
   digits 1...0 and lower case letters.

-  09.02.17  Original   By: ACRM
*/
char GetNextChainLabel(char chainLabel)
{
   if((chainLabel >= 'A') && (chainLabel < 'Z'))
   {
      chainLabel++;
   }
   else if(chainLabel == 'Z')
   {
      chainLabel = '1';
   }
   else if((chainLabel >= '1') && (chainLabel < '9'))
   {
      chainLabel++;
   }
   else if(chainLabel == '9')
   {
      chainLabel = '0';
   }
   else if(chainLabel == '0')
   {
      chainLabel = 'a';
   }
   else if((chainLabel >= 'a') && (chainLabel <'z'))
   {
      chainLabel++;
   }
   else
   {
      chainLabel = 'A';
      fprintf(stderr, "Warning (pdbsymm): More than 61 chains so reusing chain names!\n");
   }
   
   return(chainLabel);
}



/************************************************************************/
/*>void ApplyMatrixAndWriteCopy(FILE *out, PDB *pdb, char oldChain[8], 
                                char newChain[8], REAL matrix[3][3], 
                                REAL trans[3])
   --------------------------------------------------------------------
*//**
   \param[in]    *out       Output file pointer
   \param[in]    *pdb       PDB linked list
   \param[in]    oldChain   The name of the chain that we are moving and 
                            writing
   \param[in]    newChain   The name that we will use for the new chain
   \param[in]    matrix[][] The rotation matrix
   \param[in]    trans[]    The translation vector

   Applies the rotation matrix followed by the translation vector to the
   specified chain writing it out with the new chain label.

-  09.02.17  Original   By: ACRM
*/
void ApplyMatrixAndWriteCopy(FILE *out, PDB *pdb, char oldChain[8], 
                             char newChain[8], REAL matrix[3][3], 
                             REAL trans[3])
{
   PDB *pdb2;
   
   if((pdb2 = blGetPDBChainAsCopy(pdb, oldChain))!=NULL)
   {
      VEC3F transVec;
      PDB   *p;
      
      transVec.x = trans[0];
      transVec.y = trans[1];
      transVec.z = trans[2];
      
      /* Apply the rotations                                            */
      blApplyMatrixPDB(pdb2, matrix);
      /* And the translation                                            */
      blTranslatePDB(pdb2, transVec);
      
      /* Reset the chain label                                          */
      for(p=pdb2; p!=NULL; NEXT(p))
      {
         strcpy(p->chain, newChain);
      }

      /* Write the new PDB file                                         */
      blWritePDB(out,pdb2);

      FREELIST(pdb2, PDB);
   }
}

   
/************************************************************************/
/*>void Usage(void)
   ----------------
*//**
   Prints a usage message

-  09.02.17 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\npdbsymm V1.0 (c) 2017 Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"Usage: pdbsymm [in.pdb [out.pdb]]\n");

   fprintf(stderr,"\nI/O is to stdin/stdout if not specified\n\n");
   fprintf(stderr,"Applies non crystollographic symmetry to a PDB file \
given REMARK 350\n");
   fprintf(stderr,"(BIOMT) records in the PDB file.\n");
   fprintf(stderr,"\n");
   fprintf(stderr,"Note 1: this code currently ends up with PDB files \
that are not fully\n");
   fprintf(stderr,"valid:\n");
   fprintf(stderr,"- the headers don't reflect the additional chains\n");
   fprintf(stderr,"- there is no MASTER or CONECT record\n");
   fprintf(stderr,"- HETATMs are not all moved to the end of the \
file\n");
   fprintf(stderr,"\n");
   fprintf(stderr,"Note 2: this code only works with single character \
chain names. It \n");
   fprintf(stderr,"needs updating to deal with multi-character names!\n");
   fprintf(stderr,"\n");
}


