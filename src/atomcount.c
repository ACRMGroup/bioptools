/************************************************************************/
/**

   \file       atomcount.c
   
   \version    V1.4
   \date       19.08.14
   \brief      Count atoms neighbouring each atom in a PDB file
               Results output in B-val column
   
   \copyright  (c) Dr. Andrew C. R. Martin 1994-2014
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0  05.07.94 Original    By: ACRM
-  V1.1  24.08.94 Changed to call OpenStdFiles()
-  V1.2  30.04.08 Added -c (contact) and -n (normcontact). Now strips
                  water by default (-w to keep them)
-  V1.3  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.4  19.08.14 Fixed call to renamed function: blStripWatersPDBAsCopy()
                  By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define MAXCONTACT 160
#define DEFRAD  ((REAL)5.0)
#define DEFCRAD ((REAL)3.5)
#define TYP_ALL         0  /* Atom counting schemes                     */
#define TYP_DIFFRES     1
#define TYP_NONBOND     2
#define TYP_CONTACT     3
#define TYP_NORMCONTACT 4

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *radius, int *CountType, BOOL *StripWater);
void CountNeighbours(PDB *pdb, REAL radius, int CountType);
void Usage(void);
void doResidueContacts(PDB *pdb, REAL RadSq, int CountType);
BOOL ResSep(PDB *pdb, PDB *pr, PDB *qr);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for neighbour counting

-  05.07.94 Original    By: ACRM
-  24.08.94 Changed to call OpenStdFiles()
-  30.04.08 Now strips water by default - added option to keep them
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  19.08.14 Fixed call to renamed function blStripWatersPDBAsCopy() 
            By: CTP
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   REAL radius = DEFRAD;
   PDB  *pdb;
   int  natoms,
        CountType;
   BOOL StripWater = TRUE;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &radius, &CountType,
                   &StripWater))
   {
      /* Square the radius to save on distance sqrt()s                  */
      radius *= radius;
      
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDB(in,&natoms)) != NULL)
         {
            if(StripWater)
            {
               PDB *pdb2;
               int natoms;
               
               pdb2 = blStripWatersPDBAsCopy(pdb, &natoms);
               FREELIST(pdb, PDB);
               pdb = pdb2;
            }
            CountNeighbours(pdb, radius, CountType);
            blWritePDB(out, pdb);
         }
         else
         {
            fprintf(stderr,"No atoms read from PDB file\n");
         }
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
                     REAL *radius, int *CountType, BOOL *StripWater)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *radius      Neighbour radius
   \param[out]     *CountType   Counting scheme
   \param[out]     *StripWater  Strip waters?
   \return                     Success?

   Parse the command line
   
-  05.07.94 Original    By: ACRM
-  29.04.08 Added -c and -n handling
-  30.04.08 Added -w handling
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *radius, int *CountType, BOOL *StripWater)
{
   BOOL GotRad;

   argc--;
   argv++;

   GotRad = FALSE;
   *StripWater = TRUE;
   *CountType = TYP_ALL;
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'r':
            argc--;
            argv++;
            sscanf(argv[0],"%lf",radius);
            GotRad = TRUE;
            break;
         case 'd':
            *CountType = TYP_DIFFRES;
            break;
         case 'b':
            *CountType = TYP_NONBOND;
            break;
         case 'c':
            *CountType = TYP_CONTACT;
            if(!GotRad)
            {
               *radius = DEFCRAD;
            }
            break;
         case 'n':
            *CountType = TYP_NORMCONTACT;
            if(!GotRad)
            {
               *radius = DEFCRAD;
            }
            break;
         case 'w':
            *StripWater = FALSE;
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
/*>void CountNeighbours(PDB *pdb, REAL RadSq, int CountType)
   ---------------------------------------------------------
*//**

   \param[in,out]  *pdb       PDB linked list
   \param[in]      RadSq      Radius squared for neighbour search
   \param[in]      CountType  Counting scheme

   Does the actual work of counting the neighbours. 5 schemes are allowed:
   TYP_ALL         All atoms counted
   TYP_DIFFRES     Only atoms in different residues are counted
   TYP_NONBOND     Only atoms > 2.0A away are counted
   TYP_CONTACT     Counts number of residues which contact each residue
   TYP_NORMCONTACT Counts number of residues which contact each residue
                   and normalize by number of atoms in this residue

-  05.07.94 Original    By: ACRM
-  29.04.08 Added TYP_CONTACT / TYP_NORMCONTACT
*/
void CountNeighbours(PDB *pdb, REAL RadSq, int CountType)
{
   PDB *p,
       *q;
   int count;

   if((CountType == TYP_CONTACT) || (CountType == TYP_NORMCONTACT))
   {
      doResidueContacts(pdb, RadSq, CountType);
   }
   else
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         count = 0;
         for(q=pdb; q!=NULL; NEXT(q))
         {
            /* Skip this comparison if the appropriate conditions apply */
            switch(CountType)
            {
            case TYP_ALL:
               if(p==q) continue;
               break;
            case TYP_DIFFRES:
               if(p->resnum    == q->resnum    &&
                  p->insert[0] == q->insert[0] &&
                  p->chain[0]  == q->chain[0])
                  continue;
               break;
            case TYP_NONBOND:
               /* 29.04.08 Corrected to <4.0 rather than >4.0 !!!       */
               if((p==q) || (DISTSQ(p,q) < (REAL)4.0))
                  continue;
               break;
            }
            
            if(DISTSQ(p,q) < RadSq)
               count++;
         }
         p->bval = (REAL)count;
      }
   }
}

/************************************************************************/
/*>void doResidueContacts(PDB *pdb, REAL RadSq, int CountType)
   -----------------------------------------------------------
*//**

   \param[in]      *pdb        PDB linked list
   \param[in]      RadSq       Squared cutoff distance
   \param[in]      CountType   Counting scheme

   Does residue-by-residue contacts rather than atom-atom contacts
   Allowed counting schemes are
   TYP_CONTACT     Counts number of residues which contact each residue
   TYP_NORMCONTACT Counts number of residues which contact each residue
                   and normalize by number of atoms in this residue
   (though no check is made for invalid types which are treated as
   TYP_CONTACT)

-  29.04.08  Original   By: ACRM
*/

void doResidueContacts(PDB *pdb, REAL RadSq, int CountType)
{
   PDB *p, *q,
       *res_p, *res_q,
       *next_res_p = NULL,
       *next_res_q = NULL;
   
   int atomcount = 0,
       contacts = 0;

   /* Step through each residue                                         */
   for(res_p=pdb; res_p!=NULL; res_p = next_res_p)
   {
      next_res_p=blFindNextResidue(res_p);
      /* Using the occupancy as a flag, clear the contact list          */
      for(p=pdb; p!=NULL; NEXT(p))
      {
         p->occ = 0.0;
      }
      /* Count of atoms in this residue                                 */
      atomcount = 0;
      contacts = 0;
         
      
      /* Step through atoms in this residue                             */
      for(p=res_p; p!=next_res_p; NEXT(p))
      {
         atomcount++;   /* Increment atoms in this residue              */

         /* Step through the other residues                             */
         for(res_q=pdb; res_q!=NULL; res_q = next_res_q)
         {
            next_res_q=blFindNextResidue(res_q);

            /* Check it's a different and not-bonded residue            */
            if(ResSep(pdb, res_p, res_q))
            {
               /* Run through the atoms in this residue                 */
               for(q=res_q; q!=next_res_q; NEXT(q))
               {
                  REAL distsq;
                  distsq = DISTSQ(p,q);
                  if((distsq < RadSq) && (distsq > (REAL)4.0))
                  {
                     q->occ = 1.0;
                  }
               }
            }
         }
      }  /* Atoms in res_p                                              */

      /* Step through the other residues again and count how many
         residues make contact
      */
      for(res_q=pdb; res_q!=NULL; res_q = next_res_q)
      {
         next_res_q=blFindNextResidue(res_q);

         /* Run through the atoms in this residue                       */
         for(q=res_q; q!=next_res_q; NEXT(q))
         {
            /* As soon as we find one atom flagged in this residue,
               bump the counter and jump out of the loop
            */
            if(q->occ > (REAL)0.5)
            {
               contacts++;
               break;
            }
         }  /* Atoms in res_q                                           */
      }  /* res_q                                                       */

      /* Step through atoms in this residue and update the b-value      */
      for(p=res_p; p!=next_res_p; NEXT(p))
      {
         if(CountType == TYP_NORMCONTACT)
         {
            p->bval = (REAL)contacts / (REAL)atomcount;
         }
         else
         {
            p->bval = (REAL)contacts;
         }
      }  /* Atoms in res_p                                              */
   }  /* res_p                                                          */

   /* Reset the occupancy                                               */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->occ = 1.0;
   }
}


/************************************************************************/
/*>BOOL ResSep(PDB *pdb, PDB *pr, PDB *qr)
   ---------------------------------------
*//**

   \param[in]      *pdb       PDB linked list
   \param[in]      *pr        Pointer to first residue of interest
   \param[in]      *qr        Pointer to second residue of interest
   \return                     Are the residues separated?

-  29.04.08  Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
BOOL ResSep(PDB *pdb, PDB *pr, PDB *qr)
{
   int sep;
   PDB *p, *q;
   BOOL GotPR = FALSE,
        GotQR = FALSE;
   
   /* If they are more than 1 resnum apart immediately return TRUE      */
   if(pr->resnum < (qr->resnum - 1))
      return(TRUE);
   if(pr->resnum > (qr->resnum + 1))
      return(TRUE);

   /* Otherwise count through the PDB linked list and see how far 
      apart they are
   */
   sep = 0;
   for(p=pdb; p!=NULL; p = q)
   {
      q = blFindNextResidue(p);
      if(p == pr)
      {
         if(!GotQR)
         {
            sep = 0;
         }
         GotPR = TRUE;
         if(sep > 1)
         {
            return(TRUE);
         }
      }
      if(p == qr)
      {
         if(!GotPR)
         {
            sep = 0;
         }
         GotQR = TRUE;
         if(sep > 1)
         {
            return(TRUE);
         }
      }
      sep++;
   }
   return(FALSE);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  05.07.94 Original    By: ACRM
-  24.08.94 V1.1
-  29.04.08 V1.2
-  22.07.14 V1.3 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\nAtomCount V1.3 (c) 1994-2014, Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Usage: atomcount [-r <rad>] [-d|-b|-c|-n] [-w] \
[<in.pdb>] [<out.pdb>]\n");
   fprintf(stderr,"                 -r Specify radius (Default: %.2f or \
%.2f with -c/-n)\n", DEFRAD, DEFCRAD);
   fprintf(stderr,"                 -d Ignore atoms in current \
residue\n");
   fprintf(stderr,"                 -b Ignore bonded atoms (<2.0A)\n");
   fprintf(stderr,"                 -c Count residue contacts\n");
   fprintf(stderr,"                 -n Normalized residue contacts\n");
   fprintf(stderr,"                 -w Keep waters\n\n");

   fprintf(stderr,"Counts the number of atoms within the specified \
radius of each atom in\n");
   fprintf(stderr,"a PDB structure. The results are placed in the \
B-value column.\n\n");
   fprintf(stderr,"With residue contacts, the number of residues which \
make contact with\n");
   fprintf(stderr,"the current residue is calculated. The residues \
either side of the\n");
   fprintf(stderr,"current residue are not included in the count. When \
normalized, the\n");
   fprintf(stderr,"residue contact counts are divided by the number of \
atoms in the \n");
   fprintf(stderr,"current residue.\n\n");
}

