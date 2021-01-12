/*************************************************************************

   Program:    rangecontacts
   File:       rangecontacts.c
   
   Version:    V1.1
   Date:       12.01.21
   Function:   Finds residues contacting a specified range of residues
   
   Copyright:  (c) Prof. Andrew C. R. Martin 2020-21
   Author:     Prof. Andrew C. R. Martin
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  26.03.20 Original
   V1.1  12.01.21 Added -i for internal contacts, -c to show
                  counts of contacts and -m for mainchain contacts

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define DEF_RAD 3.0
#define BONDED_NOT 0
#define BONDED_NTER 1
#define BONDED_CTER 2

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *radsq, char *startres, char *stopres,
                  BOOL *doMainChain, BOOL *doInternal, BOOL *showCounts);
void DoAnalysis(FILE *out, PDB *pdb, REAL RadSq, char *startres,
                char *stopres, BOOL doMainChain, BOOL doInternal,
                BOOL showCounts);
int MakesContact(PDB *res, PDB *nextRes, PDB *pStart, PDB *pStop,
                 BOOL doMainChain, BOOL showCounts, REAL RadSq,
                 BOOL isInternal);
void PrintContact(FILE *out, PDB *p, BOOL showCounts, int nContacts);
BOOL IsSidechain(PDB *p);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for performing contact analysis

   26.03.20 Original    By: ACRM
   12.01.21 Added -m -c -i
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        startres[MAXBUFF],
        stopres[MAXBUFF];
   PDB  *pdb;
   int  natom;
   FILE *in         = stdin,
        *out        = stdout;
   REAL radsq       = DEF_RAD * DEF_RAD;
   BOOL doMainChain = FALSE,
        doInternal  = FALSE,
        showCounts  = FALSE;

   
   if(ParseCmdLine(argc, argv, infile, outfile, &radsq,
                   startres, stopres, &doMainChain, &doInternal,
                   &showCounts))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDBAtoms(in, &natom))!=NULL)
         {
            DoAnalysis(out, pdb, radsq, startres, stopres,
                       doMainChain, doInternal, showCounts);
         }
         else
         {
            fprintf(stderr,"Warning: No atoms read from PDB file\n");
         }
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
                     REAL *radsq, char *startres, char *stopres,
                     BOOL *doMainChain, BOOL *doInternal, 
                     BOOL *showCounts)
   -----------------------------------------------------------------------
   Input:   int      argc         Argument count
            char     **argv       Argument array
   Output:  char     *infile      Input filename (or blank string)
            char     *outfile     Output filename (or blank string)
            REAL     *radsq       Contact radius squared
            char     *startres    Starting residue of range
            char     *stopres     Starting residue of range
            BOOL     *doMainChain Include mainchain atoms
            BOOL     *doInternal  Include internal atoms in the range
            BOOL     *showCounts  Show counts of contacts

   Returns: BOOL                  Success

   Parse the command line

   26.03.20 Original    By: ACRM
   12.01.21 Added -m -i -c
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *radsq, char *startres, char *stopres,
                  BOOL *doMainChain, BOOL *doInternal, BOOL *showCounts)
{
   argc--;
   argv++;
   
   infile[0]   = outfile[0] = '\0';
   startres[0] = stopres[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'i':
            *doInternal = TRUE;
            break;
         case 'm':
            *doMainChain = TRUE;
            break;
         case 'c':
            *showCounts = TRUE;
            break;
         case 'r':
            argv++;
            argc--;
            sscanf(argv[0],"%lf",radsq);
            (*radsq) *= (*radsq);
            break;
         default:
            return(FALSE);
            break;
         }
         argc--;
         argv++;
      }
      else
      {
         /* Check that there are only 2-4 arguments left                */
         if((argc < 2) || (argc > 4))
            return(FALSE);

         strncpy(startres, argv[0], MAXBUFF);
         argc--; argv++;
         strncpy(stopres,  argv[0], MAXBUFF);
         argc--; argv++;

         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--; argv++;
            if(argc)
            {
               strcpy(outfile, argv[0]);
               argc--; argv++;
            }
         }
      }
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void DoAnalysis(FILE *out, PDB *pdb, REAL RadSq, char *startres,
                   char *stopres, BOOL doMainChain, BOOL doInternal,
                   BOOL showCounts)
   ----------------------------------------------------------------
   Input:   FILE    *out        Output file pointer
            PDB     *pdb        PDB linked list
            REAL    RadSq       Squared radius for contact
            char    *startres   Start residue for range
            char    *stopres    Stop residue for range
            BOOL    doMainChain Include mainchain atoms
            BOOL    doInternal  Include residues within the range
            BOOL    showCounts  Display the contact

   Main routine to do the contacts analysis between the range and the 
   rest of the protein

   26.03.20 Original   By: ACRM
   12.01.21 Added doMainChain, doInternal and showCounts
*/   
void DoAnalysis(FILE *out, PDB *pdb, REAL RadSq, char *startres,
                char *stopres, BOOL doMainChain, BOOL doInternal,
                BOOL showCounts)
{
   PDB  *p,
        *pStart, *pLast, *pStop,
        *nextRes;

   /* Find the first residue of the range                               */
   if((pStart = blFindResidueSpec(pdb, startres))==NULL)
   {
      fprintf(stderr,"Error (rangecontacts) - Residue not found: %s\n",
              startres);
      return;
   }
   
   /* and the one after the last residue                                */
   if((pLast = blFindResidueSpec(pdb, stopres))==NULL)
   {
      fprintf(stderr,"Error (rangecontacts) - Residue not found: %s\n",
              stopres);
      return;
   }
   
   pStop  = blFindNextResidue(pLast);
   
   /* Run through the linked list a residue at a time                   */
   for(p=pdb; p!=NULL; p=nextRes)
   {
      nextRes = blFindNextResidue(p);

      if(doInternal || !blInPDBZoneSpec(p, startres, stopres))
      {
         int  nContacts;
         BOOL isInternal = FALSE;

         if(doInternal)
            isInternal = blInPDBZoneSpec(p, startres, stopres);

         if((nContacts = MakesContact(p, nextRes, pStart, pStop,
                                      doMainChain, showCounts,
                                      RadSq, isInternal)) > 0)
         {
            PrintContact(out, p, showCounts, nContacts);
         }
      }
   }
}


/************************************************************************/
/*>void PrintContact(FILE *out, PDB *p, BOOL showCounts, int nContacts)
   ----------------------------------------------------------------------
   Input:   FILE    *out       Output file pointer
            PDB     *p         Residue which makes contact
            BOOL    showCounts Print number of contacts made
            int     nContacts  Number of contacts made

   26.03.20 Original    By: ACRM
   12.01.21 Added printing of the number of contacts
*/
void PrintContact(FILE *out, PDB *p, BOOL showCounts, int nContacts)
{
   if(showCounts)
   {
      fprintf(out,"%s%d%s %d\n",
              p->chain, p->resnum, p->insert, nContacts);
   }
   else
   {
      fprintf(out,"%s%d%s\n", p->chain, p->resnum, p->insert);
   }
}


/************************************************************************/
/*>BOOL IsSidechain(PDB *p)
   ------------------------
   Input:   PDB    *p    PDB pointer
   Returns: BOOL         Is this a sidechain atom?

   Tests if an atom is part of a sidehchain.

   26.03.20 Original    By: ACRM
*/
BOOL IsSidechain(PDB *p)
{
   if(!strncmp(p->atnam,"N   ",4) ||
      !strncmp(p->atnam,"CA  ",4) ||
      !strncmp(p->atnam,"C   ",4) ||
      !strncmp(p->atnam,"O   ",4) ||
      !strncmp(p->atnam,"OXT ",4) ||
      !strncmp(p->atnam,"O1  ",4) ||
      !strncmp(p->atnam,"O2  ",4))
   {
      return(FALSE);
   }

   return(TRUE);
}


/************************************************************************/
/*>int MakesContact(PDB *testResStart, PDB *testResStop, PDB *rangeStart,
                    PDB *rangeStop, BOOL doMainChain, BOOL showCounts,
                    REAL RadSq, BOOL isInternal)
   -----------------------------------------------------------------

   testResStart/testResStop  is the residue we are looking at for contacts
   rangeStart/rangeStop      is the residue range with which we are 
                             looking for contacts

   26.03.20 Original    By: ACRM
   12.01.20 Now works a residue at a time to support mainchain and 
            internal contacts as well as printing.
*/
int MakesContact(PDB *testResStart, PDB *testResStop,
                 PDB *rangeStart, PDB *rangeStop,
                 BOOL doMainChain, BOOL showCounts,
                 REAL RadSq, BOOL isInternal)
{
   PDB  *p, *q,
        *rangeRes     = NULL,
        *nextRangeRes = NULL;
   int  NContacts     = 0,
        bonded        = BONDED_NOT;
   REAL distSq;
   
   
   /* Step through atoms in the residue of interest                     */
   for(p=testResStart; p!=testResStop; NEXT(p))
   {
      if(doMainChain || IsSidechain(p))
      {
         /* Step through residues in the range we are looking at        */
         for(rangeRes=rangeStart;
             rangeRes!=rangeStop;
             rangeRes=nextRangeRes)
         {
            nextRangeRes = blFindNextResidue(rangeRes);
            
            /* If we are looking at internal residues and the current
               test residue is the same as the current range residue
               then continue to the next range residue - i.e. don't
               look for contacts within a residue
            */
            if(rangeRes == testResStart)
               continue;
            
            if(doMainChain)
            {
               bonded = BONDED_NOT;
               /* If this range residue is the one after the current test 
                  residue and they are in the same chain then it's the
                  residue bonded to the N-terminus of the range.
                  This deals with handling of internal residues too.
               */
               if((rangeRes == testResStop) &&
                  PDBCHAINMATCH(rangeRes, testResStop))
               {
                  bonded = BONDED_NTER;
               }
               /* If the current test residue is the one after this range
                  residue and they are in the same chain then it's the
                  residue bonded to the C-terminus of the range.
                  This deals with handling of internal residues too.
               */
               else if((testResStart == nextRangeRes) &&
                       PDBCHAINMATCH(testResStart, nextRangeRes))
               {
                  bonded = BONDED_CTER;
               }
            }
            
            /* Step through atoms in this residue from the range        */
            for(q=rangeRes; q!=nextRangeRes; NEXT(q))
            {
               /* If the current test residue is bonded to the 
                  N-terminus and it's the C of this test residue and
                  the N of the range residue then ignore the 
                  interaction
               */
               if((bonded == BONDED_NTER) &&
                  !strncmp(p->atnam, "C   ", 4) &&
                  !strncmp(q->atnam, "N   ", 4)) 
               {
#ifdef DEBUG
                  fprintf(stderr,"> NTer bond %s%d%s.%s : %s%d%s.%s\n",
                          p->chain, p->resnum, p->insert, p->atnam,
                          q->chain, q->resnum, q->insert, q->atnam);
#endif
                  
                  continue;
               }
               
               
               /* If the current test residue is bonded to the 
                  C-terminus and it's the N of this test residue and
                  the C of the range residue then ignore the 
                  interaction
               */
               if((bonded == BONDED_CTER) &&
                  !strncmp(p->atnam, "N   ", 4) &&
                  !strncmp(q->atnam, "C   ", 4))
               {
#ifdef DEBUG
                  fprintf(stderr,"> CTer bond %s%d%s.%s : %s%d%s.%s\n",
                          p->chain, p->resnum, p->insert, p->atnam,
                          q->chain, q->resnum, q->insert, q->atnam);
#endif
                  
                  continue;
               }
               
               
               /* Otherwise check the distance for a contact            */
               distSq = DISTSQ(p, q);
               if(distSq <= RadSq)
               {
#if (DEBUG > 2)
                  fprintf(stderr,"> %s%d%s.%s : %s%d%s.%s : %f\n",
                          p->chain, p->resnum, p->insert, p->atnam,
                          q->chain, q->resnum, q->insert, q->atnam,
                          sqrt(distSq));
#endif
                  if(!showCounts)
                     return(1);
                  NContacts++;
               }
            }
         }
      }
   }

   return(NContacts);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Print usage message.

   26.03.20 Original    By: ACRM
   12.01.21 V1.1
*/
void Usage(void)
{
   fprintf(stderr,"\nRangeContacts V1.1 (c) 2020-21, Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"Usage: rangecontacts [-r radius][-i][-c] startres \
stopres [in.pdb [out.dat]]\n");
   fprintf(stderr,"       -r Specify contact radius (Default: %.3f)\n",
           DEF_RAD);
   fprintf(stderr,"       -i Do internal contacts within the range as \
well\n");
   fprintf(stderr,"       -c Count and display the number of contacts \
made by each residue\n");
   fprintf(stderr,"       -m Include mainchain as well as sidechain \
atoms\n");
   
   fprintf(stderr,"\nI/O is through stdin/stdout if files are not \
specified.\n\n");
   fprintf(stderr,"Performs a contact analysis at the residue level to \
find residues whose\n");
   fprintf(stderr,"sidechains contact any atom of the residues specified \
in the given range.\n");
   fprintf(stderr,"When used with -m, C-N bonds between adjacent residues \
are ignored.\n");
   fprintf(stderr,"When used with -i, contacts within a residue are \
ignored.\n\n");
}
      

