/*************************************************************************

   Program:    rangecontacts
   File:       rangecontacts.c
   
   Version:    V1.3
   Date:       28.10.15
   Function:   Finds residues contacting a specified range of residues
   
   Copyright:  (c) Prof. Andrew C. R. Martin 2020
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

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *radsq, char *startres, char *stopres);
void DoAnalysis(FILE *out, PDB *pdb, REAL RadSq, char *startres,
                char *stopres);
BOOL MakesContact(PDB *res, PDB *nextRes, PDB *pStart, PDB *pStop,
                  REAL RadSq);
void PrintContact(FILE *out, PDB *p);
BOOL IsSidechain(PDB *p);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for performing contact analysis

   26.03.20 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        startres[MAXBUFF],
        stopres[MAXBUFF];
   PDB  *pdb;
   int  natom;
   FILE *in        = stdin,
        *out       = stdout;
   REAL radsq      = DEF_RAD * DEF_RAD;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &radsq,
                   startres, stopres))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDBAtoms(in, &natom))!=NULL)
         {
            DoAnalysis(out, pdb, radsq, startres, stopres);
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
/*>void Usage(void)
   ----------------
   Print usage message.

   17.10.95 Original    By: ACRM
   04.03.15 V1.2
   28.10.15 V1.3
*/
void Usage(void)
{
   fprintf(stderr,"\nRangeContacts V1.0 (c) 2020, Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"Usage: rangecontacts [-r radius] startres stopres \
[in.pdb [out.dat]]\n");
   fprintf(stderr,"       -r Specify contact radius (Default: %.3f)\n\n",
           DEF_RAD);

   fprintf(stderr,"I/O is through stdin/stdout if files are not \
specified.\n\n");
   fprintf(stderr,"Performs a contact analysis at the residue level to \
find residues whose\n");
   fprintf(stderr,"sidechains contact any atom of the residues specified \
in the given range.\n\n");
}
      

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     REAL *radsq, char *startres, char *stopres)
   ---------------------------------------------------------------------
   Input:   int      argc        Argument count
            char     **argv      Argument array
   Output:  char     *infile     Input filename (or blank string)
            char     *outfile    Output filename (or blank string)
            REAL     *radsq      Contact radius squared
            char     *startres   Starting residue of range
            char     *stopres    Starting residue of range
   Returns: BOOL                 Success

   Parse the command line

   26.03.20 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *radsq, char *startres, char *stopres)
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
                   char *stopres)
   ----------------------------------------------------------------
   Input:   FILE    *out       Output file pointer
            PDB     *pdb       PDB linked list
            REAL    RadSq      Squared radius for contact
            char    *startres  Start residue for range
            char    *stopres   Stop residue for range

   Main routine to do the contacts analysis between the range and the 
   rest of the protein

   26.03.20 Original (based on code from rangecontacts.c)   By: ACRM
*/   
void DoAnalysis(FILE *out, PDB *pdb, REAL RadSq, char *startres,
                char *stopres)
{
   PDB  *p,
      *pStart, *pStop,
      *nextRes;
   






   /* Find the first residue of the range                               */
   pStart = blFindResidueSpec(pdb, startres);
   /* and the one after the last residue                                */
   pStop  = blFindResidueSpec(pdb, stopres);
   pStop  = blFindNextResidue(pStop);
   
   /* Run through the linked list a residue at a time                   */
   for(p=pdb; p!=NULL; p = nextRes)
   {
      nextRes = blFindNextResidue(p);

      if(!blInPDBZoneSpec(p, startres, stopres))
      {
         if(MakesContact(p, nextRes, pStart, pStop, RadSq))
         {
            PrintContact(out, p);
         }
      }
   }
}

/************************************************************************/
BOOL MakesContact(PDB *res, PDB *nextRes, PDB *pStart, PDB *pStop,
                  REAL RadSq)
{
   PDB *p, *q;
   
   for(p=res; p!=nextRes; NEXT(p))
   {
      if(IsSidechain(p))
      {
         for(q=pStart; q!=pStop; NEXT(q))
         {
            if(DISTSQ(p, q) <= RadSq)
            {
               return(TRUE);
            }
         }
      }
      
   }
   return(FALSE);
   
}

/************************************************************************/
/*>void PrintContact(FILE *out, PDB *p)
   ---------------------------------------------------------------
   Input:   FILE    *out       Output file pointer
            PDB     *p         Residue which maes contact
*/
void PrintContact(FILE *out, PDB *p)
{
   fprintf(out,"%s%d%s\n", p->chain, p->resnum, p->insert);
}

/************************************************************************/
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
