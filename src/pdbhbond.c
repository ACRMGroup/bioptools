/*** TODO
Finish commenting
Handle PGP file in data directory
 ***/
/************************************************************************/
/**

   \file       pdbhbond.c
   
   \version    V2.1
   \date       08.09.17
   \brief      List hydrogen bonds
   
   \copyright  (c) UCL, Dr. Andrew C.R. Martin, 2014-2017
   \author     Dr. Andrew C.R. Martin
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
   Displays a list of hydrogen bonds based on calculated distances and
   angles according to the Baker and Hubbard criteria


   Rules for finding Hydrogen Bonds
   ================================
   Within the protein we apply the rules exactly as defined by Baker and
   Hubbard.

   If the location of H is known:
      D-H-A and H-A-P are between 90 and 180 degrees
      H-A < 2.5A

   If the location of H is not known:
      D-A-P is between 90 and 180 degrees
      D-A < 3.5A

   Note that the code now uses 3.35A instead of 3.5A - this can be set 
   on the command line.

   For protein ligand interactions, we apply the following rules:

   In FindProtLigandHBonds()
   -------------------------
   1. Look for atom pairs which are
      (i)  Non-protein/nucleotide -to- protein/nucleotide
      (ii) Nucleotide             -to- protein
   In TestForHBond()
   -----------------
   2. See if one is a donor and the other an acceptor
   3. Check the distance is OK
   4. Find the antecedent for the acceptor
   5. Count the antecedents for the donor.
   6. If the atom has the potential to accept HBonds (i.e. number of
      antecedents is less than the max) then test the geometry of the 
      HBond


  Rules for finding non-bonded contacts
  =====================================
  1. For atom pairs within a distance cutoff: 2.7 - 3.35A (centre-centre)
  2. a) Not covalently bonded
     b) Not H-bonded
     c) Not in the same residue

   The following rules (John Mitchell, personal communication) may be
   useful in future:
   * H-bond donors:
   *    Any O, S or N potentially donates the number of Hs bound to it.
   *
   * H-bond acceptors:
   *    O/S: sp3 - can accept up to 2
   *         sp2 - can accept up to 2
   *
   *    N:   planar   - no acceptance
   *         sp3      - can accept 1
   *         sp2      - can accept 1
   *         aromatic - can accept (3 - no. of bonds)
   *         amide    - no acceptance

**************************************************************************

   Notes and issues:
   =================

   In some structures extra HBonds are being formed
   Asp95.OD1 - Sam400.O2*
   Asp95.OD2 - Sam400.O3*
   Asp95.OD1 - Sam400.O3*

   If you look at the structure the Sam400.O3* could actually
   be making an HBond with either OD1 or OD2 of Asp95. Clearly
   OD2 is a better bet (especially since OD1 is binding to O2*)
   but the code is not designed to look ahead for things like
   this. Seems that the extra one is using ASP as a Donor

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-   V1.0  07.06.99 Original   By: ACRM
-   V1.1  04.11.99 Improved checking of potential oxygen donors to see
                   if the antecedent carbon was trigonal planar
                   rather than a tetrahedral. Now checks two 
                   pre-antecedent angles rather than one to avoid cases 
                   where one angle is just above our 115degree threshold.
-   V1.2  12.05.00 plhbonds are now output with a 'relaxed' field - this
                   is 0 for H-bonds which pass all the tests we use and
                   1 if the H-bond fails on the basis of geometry of a
                   donor oxygen (i.e. this oxygen may not be a donor)
-   V2.0  20.07.15 Now a standalone program that uses PDB files
                   rather than based on XMAS. Now uses internal PDB 
                   CONECT information rather than keeping its own version
                   of the CONECT data
-   V2.1  08.09.17 Changed comment in output and spacing of fields

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/hbond.h"
#include "bioplib/hash.h"
#include "bioplib/angle.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/

/*** OPTIONS                                                          ***/
#define NO_WATERS_IN_MOL_LIST  /* Don't treat waters as separate mols   */
#undef  USE_METAL_IN_POLYHET   /* Metals can be treated as part of a
                                  polymer of HETATMs                    */
/*** END-OPTIONS                                                      ***/


#define MAXBUFF          160
#define MAX_TYPE_STRING   24
#define MAX_NAME_STRING  180
#define MAX_CHAIN_STRING   8
#define MAX_START_STRING   8
#define MAX_PEPTIDE_LENGTH 30

#define BOND_TOL             DEFCONECTTOL   /* Tolerance for a bond     */

#define MAXHBONDDISTSQ       11.2225        /* 3.35A max HBond distance */
#define MAXBONDSQ             3.0           /* 1.732A max bond distance */
#define MINNBDISTSQ           8.41          /* 2.9A min NB distance     */
#define MAXNBDISTSQ          15.21          /* 3.9A max NB distance     */
#define MAXTETRAHEDRALANGLE 115.0           /* Max allowed angle for a
                                               tetrahedral carbon. Assumed
                                               to be trigonal planar if
                                               larger than this         */

/* Used to store HBonding information for each element type             */
typedef struct
{
   char *element;
   int  donor,
        acceptor;
}  HBONDING;

/* PDB.extras structure used for original atom numbers (before adding
   hydrogens) and molecule IDs                                          */
typedef struct _pdbextras
{
   int origAtnum;
   int molid;
}  PDBEXTRAS;



/************************************************************************/
/* Globals
*/
static HBONDING sHBonding[] = 
{
   /* Elements capable of true hydrogen bonds                           */
   {"N",  3, 1}, {"O",  2, 2}, {"F",  1, 3}, 
   {NULL, 0, 0} 
} ;

static HBONDING sPseudoHBonding[] =
{  
   /* Elements not able to make anything like a hydrogen bond           */
   {"C",  0, 0},  {"H",  0, 0},  {"HE",  0, 0}, {"NE",  0, 0}, 
   {"AR",  0, 0}, {"KR",  0, 0}, {"XE",  0, 0}, {"RN",  0, 0},

   /* Pseudo-donors Group IA metals                                     */
   {"LI",  6, 0}, {"NA",  6, 0}, {"K",  6, 0},  {"RB",  6, 0}, 
   {"CS",  6, 0}, {"FR",  6, 0}, 

   /* Pseudo-donors Group IIA metals                                    */
   {"BE",  6, 0}, {"MG",  6, 0}, {"CA",  6, 0}, {"SR",  6, 0}, 
   {"BA",  6, 0}, {"RA",  6, 0}, 

   /* Pseudo-acceptors Group VIB                                        */
   {"S",  0, 2},  {"SE",  0, 0}, {"TE",  0, 0}, {"PO",  0, 0}, 

   /* Pseudo-acceptors Group VIIB                                       */
   {"CL",  0, 0}, {"BR",  0, 0}, {"I",  0, 0},  {"AT",  0, 0}, 

   /* Group IIIB                                                        */
   {"B",  0, 0},  {"AL",  6, 0}, {"GA", -1, 0}, {"IN", -1, 0}, 
   {"TL", -1, 0},  

   /* Group IVB                                                         */
   {"SI",  0, 0}, {"GE", -1, 0}, {"PB", -1, 0}, {"SN", -1, 0},

   /* Group VB                                                          */
   {"P",  0, 0},  {"AS", -1, 0}, {"SB", -1, 0}, {"BI", -1, 0}, 

   /* Transition elements                                               */
   {"SC", -1, 0}, {"TI", -1, 0}, {"V", -1, 0},  {"CR", -1, 0}, 
   {"MN",  6, 0}, {"FE",  6, 0}, {"CO", -1, 0}, {"NI", -1, 0}, 
   {"CU",  6, 0}, {"ZN",  6, 0}, {"Y", -1, 0},  {"ZR", -1, 0}, 
   {"NB", -1, 0}, {"MO", -1, 0}, {"TC", -1, 0}, {"RU", -1, 0}, 
   {"RH", -1, 0}, {"PD", -1, 0}, {"AG", -1, 0}, {"CD",  6, 0}, 
   {"LA", -1, 0}, {"HF", -1, 0}, {"TA", -1, 0}, {"W",  -1, 0}, 
   {"RE", -1, 0}, {"OS", -1, 0}, {"IR", -1, 0}, {"PT", -1, 0}, 
   {"AU", -1, 0}, {"HG",  6, 0}, {"AC", -1, 0}, {"KU", -1, 0}, 

   /* Lanthanides                                                       */
   {"CE", -1, 0}, {"PR", -1, 0}, {"ND", -1, 0}, {"PM", -1, 0}, 
   {"SM", -1, 0}, {"EU", -1, 0}, {"GD", -1, 0}, {"TB", -1, 0}, 
   {"DY", -1, 0}, {"HO", -1, 0}, {"ER", -1, 0}, {"TM", -1, 0}, 
   {"YB", -1, 0}, {"LU", -1, 0}, 

   /* Actinides                                                         */
   {"TH", -1, 0}, {"PA", -1, 0}, {"U",  -1, 0}, {"NP", -1, 0}, 
   {"PU", -1, 0}, {"AM", -1, 0}, {"CM", -1, 0}, {"BK", -1, 0}, 
   {"CF", -1, 0}, {"ES", -1, 0}, {"FM", -1, 0}, {"MD", -1, 0}, 
   {"NO", -1, 0}, {"LR", -1, 0}, 

   {NULL, 0, 0} 
};


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *pgpfile, REAL *minNBDistSq, REAL *maxNBDistSq, 
                  REAL *maxHBDistSq);
HBLIST *FindProtProtHBonds(PDB *pdb);
HBLIST *FindProtLigandHBonds(PDB *pdb, PDB **pdbarray,
                             BOOL pseudo, REAL maxHBDistSq);
HBLIST *FindLigandLigandHBonds(PDB *pdb, 
                               PDB **pdbarray, BOOL pseudo,
                               REAL maxHBDistSq);
void PrintHBList(FILE *out, HBLIST *hblist, char *type, BOOL relaxed);
HBLIST *TestForHBond(PDB *pdb, PDB *p, PDB *q, PDB **pdbarray,
                     BOOL pseudo, REAL maxHBDistSq);
int IsDonor(PDB *p, BOOL allowPseudo, BOOL *pseudo);
int IsAcceptor(PDB *p, BOOL allowPseudo, BOOL *pseudo);
PDB *FindAntecedent(PDB *atom, PDB **pdbarray, int *count, int nth);
PDB *FindBondedHydrogen(PDB *pdb, PDB *donor, PDB *acceptor);
HBLIST *doTestForHBond(PDB *pdb, PDB *donor, PDB *acceptor, 
                       PDB **pdbarray, int donMax, REAL maxHBDistSq);
HBLIST *FindNonBonds(PDB *pdb, PDB **pdbarray,
                     HBLIST *hbonds, REAL minNBDistSq, REAL maxNBDistSq);
BOOL IsListedAsHBonded(PDB *p, PDB *q, HBLIST *hbonds);
BOOL isAPeptide(PDB *pdb, PDB *atm);
void SetAtomNumExtras(PDB *pdb);
BOOL UpdatePDBExtras(PDB *pdb);
BOOL SetMolecules(PDB *pdb);
void MarkLinkedResidues(PDB *chainStart, PDB *resStart, 
                        PDB *nextChain, int id);
void DeleteMetalConects(PDB *pdb);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program for calculating HBonds and non-bonds from PDB files

-  07.06.99 Original   By: ACRM
-  08.06.99 Added molecule reading
            Added ligand-ligand HBonds
-  09.06.99 Fixed joining of HBond lists where first list is NULL
-  16.06.99 Added min and max NB/HB distances as variables
-  22.07.15 Modified to use PDB files and standard BiopLib structures
            and functions
*/
int main(int argc, char **argv)
{
   FILE       *in = stdin, 
              *out = stdout,
              *pgp;
   PDB        *pdb,
              **pdbarray;
   int        nhyd,
              indexSize;
   char       infile[MAXBUFF],
              outfile[MAXBUFF],
              pgpfile[MAXBUFF];
   HBLIST     *ppHBonds = NULL,
              *plHBonds = NULL,
              *llHBonds = NULL,
              *pplHBonds = NULL,
              *nbContacts = NULL,
              *hb;
   WHOLEPDB   *wpdb = NULL;
   REAL       minNBDistSq = MINNBDISTSQ,
              maxNBDistSq = MAXNBDISTSQ,
              maxHBDistSq = MAXHBONDDISTSQ;
   STRINGLIST *warnings = NULL;
   
   if(ParseCmdLine(argc, argv, infile, outfile, pgpfile, &minNBDistSq, 
                   &maxNBDistSq, &maxHBDistSq))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         /* Open the PGP file                                           */
         if((pgp = blOpenPGPFile(pgpfile, FALSE))==NULL)
         {
            fprintf(stderr,"pdbhbond: (error) Unable to open PGP file\n");
            return(1);
         }
            
         if((wpdb = blReadWholePDB(in))==NULL)
         {
            fprintf(stderr,"pdbhbond: (error) Unable to PDB file\n");
            return(1);
         }

         if((pdb = wpdb->pdb)==NULL)
         {
            fprintf(stderr,"pdbhbond: (error) No atoms read from PDB \
file\n");
            return(1);
         }

         /* Store the original atom numbers in the extras field         */
         if(!UpdatePDBExtras(pdb))
         {
            fprintf(stderr,"pdbhbond: (error) No memory for extra PDB \
data\n");
            return(1);
         }
         
         SetAtomNumExtras(pdb);

         /* Add hydrogens to the protein                                */
         if((nhyd = blHAddPDB(pgp, pdb))==0)
         {
            fprintf(stderr,"pdbhbond: (warning) No hydrogens added to \
PDB file\n");
         }

         /* Create extras fields for the extra hydrogen atoms           */
         if(!UpdatePDBExtras(pdb))
         {
            fprintf(stderr,"pdbhbond: (error) No memory for extra PDB \
data\n");
            return(1);
         }

         if((warnings = blSetPDBAtomTypes(pdb))!=NULL)
         {
            STRINGLIST *s;
            for(s=warnings; s!=NULL; NEXT(s))
            {
               fprintf(stderr,"%s\n", s->string);
            }
         }
         
         SetMolecules(pdb);

         if((pdbarray=blIndexAtomNumbersPDB(pdb, &indexSize))==NULL)
         {
            fprintf(stderr,"pdbhbond: (error) Failed to index PDB \
data\n");
            return(1);
         }

         DeleteMetalConects(pdb);
            
         /* Find protein-protein HBonds                                 */
         blSetMaxProteinHBondDADistance((REAL)sqrt(maxHBDistSq));
         ppHBonds = FindProtProtHBonds(pdb);
         PrintHBList(out, ppHBonds, "pphbonds", FALSE);
         FREELIST(ppHBonds, HBLIST);

         /* Find protein-ligand HBonds                                  */
         plHBonds = FindProtLigandHBonds(pdb, pdbarray, FALSE,
                                         maxHBDistSq);
         PrintHBList(out, plHBonds, "plhbonds", TRUE);

         /* Find protein-ligand pseudo-HBonds                           */
         pplHBonds = FindProtLigandHBonds(pdb, pdbarray, TRUE,
                                          maxHBDistSq);
         PrintHBList(out, pplHBonds, "pseudohbonds", FALSE);

         /* Join the pseudo-HBonds list onto the protein-ligand list    */
         hb = plHBonds;
         if(hb != NULL)
         {
            LAST(hb);
            hb->next = pplHBonds;
         }
         else
         {
            plHBonds = hb = pplHBonds;
         }

         /* Find ligand-ligand HBonds                                   */
         llHBonds = FindLigandLigandHBonds(pdb, pdbarray, FALSE,
                                           maxHBDistSq);
         PrintHBList(out, llHBonds, "llhbonds", TRUE);

         /* Join the ligand-ligand HBonds list onto the previous lists  */
         if(hb!=NULL)
         {
            LAST(hb);
            hb->next = llHBonds;
         }
         else
         {
            plHBonds = hb = llHBonds;
         }

         /* Find non-bonded contacts                                    */
         nbContacts = FindNonBonds(pdb, pdbarray, plHBonds,
                                   minNBDistSq, maxNBDistSq);
         PrintHBList(out, nbContacts, "nonbonds", FALSE);

         FREEPDBEXTRAS(pdb);
         FREELIST(pdb, PDB);
      }
   }
   else
   {
      Usage();
   }
   
   
   return(0);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
*//**
   Prints a usage message

-  07.06.99 Original   By: ACRM
-  09.06.99 Added -q
-  16.06.99 Added -n, -x, -b
-  22.07.15 V2.0. Added -p

*/
void Usage(void)
{
   fprintf(stderr,"\npdbhbond V2.0 (c) 2015, Dr. Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Usage: pdbhbond [-n dist][-x dist][-b dist]\
[-p pgpfile] [infile [outfile]]\n");
   fprintf(stderr,"       -n  Minimum NBond distance (Default: %.2f)\n",
           sqrt(MINNBDISTSQ));
   fprintf(stderr,"       -x  Maximum NBond distance (Default: %.2f)\n",
           sqrt(MAXNBDISTSQ));
   fprintf(stderr,"       -b  Maximum HBond distance (Default: %.2f)\n",
           sqrt(MAXHBONDDISTSQ));
   fprintf(stderr,"       -p  Specify PGP file containing data for \
adding hydrogens\n");
   fprintf(stderr,"\nIdentifies hydrogen bonds using simple Baker and \
Hubbard rules for\n");
   fprintf(stderr,"the definition of a hydrogen bond.\n");
   fprintf(stderr,"I/O is to standard input/output if filenames are not \
specified.\n\n");
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     char *pgpfile, REAL *minNBDistSq, REAL *maxNBDistSq,
                     REAL *maxHBDistSq)
   ---------------------------------------------------------------------
*//**
   \param[in]    argc          Argument count
   \param[in]    **argv        Argument array
   \param[out]   *infile       Input filename (or blank string)
   \param[out]   *outfile      Output filename (or blank string)
   \param[out]   *pgpfile     PGP filename
   \param[out]   *minNBDistSq  Min non-bond distance
   \param[out]   *maxNBDistSq  Max non-bond distance
   \param[out]   *maxHBDistSq  Max HBond distance
   \return                     Success

   Parse the command line

-  07.06.99 Original    By: ACRM
-  09.06.99 Added -q
-  16.06.99 Added -n, -x, -b and associated parameters
-  21.07.15 Removed -q
-  22.07.15 Added -p and pgpfile
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *pgpfile, REAL *minNBDistSq, REAL *maxNBDistSq,
                  REAL *maxHBDistSq)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = pgpfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         case 'n':
            if(!(--argc))
               return(FALSE);
            argv++;
            if((sscanf(argv[0], "%lf", minNBDistSq))==0)
               return(FALSE);
            *minNBDistSq = (*minNBDistSq) * (*minNBDistSq);
            break;
         case 'x':
            if(!(--argc))
               return(FALSE);
            argv++;
            if((sscanf(argv[0], "%lf", maxNBDistSq))==0)
               return(FALSE);
            *maxNBDistSq = (*maxNBDistSq) * (*maxNBDistSq);
            break;
         case 'b':
            if(!(--argc))
               return(FALSE);
            argv++;
            if((sscanf(argv[0], "%lf", maxHBDistSq))==0)
               return(FALSE);
            *maxHBDistSq = (*maxHBDistSq) * (*maxHBDistSq);
            break;
         case 'p':
            if(!(--argc))
               return(FALSE);
            argv++;
            strncpy(pgpfile, argv[0], MAXBUFF);
            pgpfile[MAXBUFF-1] = '\0';
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 1-2 arguments left                     */
         if(argc > 2)
            return(FALSE);
         
         if(argc)
         {
            strcpy(infile, argv[0]);
         
            /* If there's another, copy it to outfile                   */
            argc--;
            argv++;
            
            if(argc)
               strcpy(outfile, argv[0]);
         }
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>HBLIST *FindProtProtHBonds(PDB *pdb)
   ------------------------------------
*//**
   \param[in]   *pdb    PDB linked list
   \return              Linked list of protein-protein HBonds

   Create a list of HBonds within the protein

-  07.06.99 Original   By: ACRM
-  22.07.15 Modified to use PDB files and standard BiopLib structures
            and functions
*/
HBLIST *FindProtProtHBonds(PDB *pdb)
{
   PDB           *p, *q;
   static HBLIST *hblist = NULL,
                 *hbl,
                 *hb;
   
   
   /* Loop through each residue                                         */
   for(p=pdb; p!=NULL; p=blFindNextResidue(p))
   {
      /* If it's a protein/nucleotide                                   */
      if(!(p->atomtype & ATOMTYPE_NONRESIDUE) &&
         (p->atomtype != ATOMTYPE_UNDEF))
      {
         
         /* Loop through each following residue                         */
         for(q=blFindNextResidue(p); q!=NULL; q=blFindNextResidue(q))
         {
            if(p==q)
               continue;
            
            /* If it's a protein/nucleotide                             */
            if(!(q->atomtype & ATOMTYPE_NONRESIDUE) &&
               (q->atomtype != ATOMTYPE_UNDEF))
            {
               /* If there is an HBond, add it to the list              */
               if((hb=blListAllHBonds(p, q))!=NULL)
               {
                  if(hblist==NULL)
                  {
                     hblist = hbl = hb;
                  }
                  else
                  {
                     hbl->next = hb;
                     hbl       = hb;
                  }
               }
            }
         }
      }
   }

   return(hblist);
}


/************************************************************************/
/*>void PrintHBList(FILE *out, HBLIST *hblist, char *type, BOOL relaxed)
   ---------------------------------------------------------------------
*//**
   \param[in]  *out    Output file pointer
   \param[in]  *hblist Linked list of atom pairs to be printed
   \param[in]  *type   Type for XMAS output
   \param[in]  relaxed Include the 'relaxed' field from the HBLIST
                       structure

   Prints a list of HBonds (or non-bonds)

-  07.06.99 Original   By: ACRM
-  24.08.99 Calls WrapPrint() on output of atom names to escape any
            double inverted commas
-  12.05.99 Added 'relaxed' handling
-  22.07.15 Modified to use PDB files and standard BiopLib structures
            and functions. Changed output format.
-  08.09.17 Changed comment in output and spacing of fields
*/
void PrintHBList(FILE *out, HBLIST *hblist, char *type, BOOL relaxed)
{
   HBLIST *h;
   

   if(hblist != NULL)
   {
      fprintf(out, "TYPE: %s\n", type);
      if(!strcmp(type, "pphbonds")||
         !strcmp(type, "pseudohbonds"))
      {
         fprintf(out, "#  datom   aatom dresnam dresid  datnam aresnam aresid  aatnam\n");
      }
      else if(!strcmp(type, "plhbonds")||
            !strcmp(type, "llhbonds"))
      {
         fprintf(out, "#  datom   aatom dresnam dresid  datnam aresnam aresid  aatnam relaxed\n");
      }
      else if(!strcmp(type, "nonbonds"))
      {
         fprintf(out, "#  atom1   atom2 resnam1 resid1  atnam1 resnam2 resid2  atnam2\n");
      }

      for(h=hblist; h!=NULL; NEXT(h))
      {
         char donResID[24],
              accResID[24];
         MAKERESID(donResID, h->donor);
         MAKERESID(accResID, h->acceptor);
         
         fprintf(out," %7d %7d %-5s   %-7s %4s   %-5s   %-7s %4s",
                 PDBEXTRASPTR(h->donor, PDBEXTRAS)->origAtnum,
                 PDBEXTRASPTR(h->acceptor, PDBEXTRAS)->origAtnum,
                 h->donor->resnam,
                 donResID,
                 h->donor->atnam,
                 h->acceptor->resnam,
                 accResID,
                 h->acceptor->atnam);

         if(relaxed)
         {
            /* 12.05.00 ACRM - print the relaxed field if required      */
            fprintf(out,"   %s", (h->relaxed?"RELAXED":""));
         }

         fprintf(out,"\n");
      }
   }
}


/************************************************************************/
/*>int IsDonor(PDB *p, BOOL allowPseudo, BOOL *pseudo)
   ---------------------------------------------------
*//**
   \param[in]   *p            PDB pointer
   \param[in]   allowPseudo   Allow pseudo HBonds?
   \param[out]  *pseudo       Is it a pseudo donor?
   \return                    Number of donor positions

   Tests whether an atom is a HBond (pseudo-)donor

-  07.06.99 Original   By: ACRM
-  22.07.15 Modified to use PDB files and standard BiopLib structures
            and functions
*/
int IsDonor(PDB *p, BOOL allowPseudo, BOOL *pseudo)
{
   int i;

   *pseudo = FALSE;

   for(i=0; sHBonding[i].element != NULL; i++)
   {
      if(!strcmp(p->element, sHBonding[i].element))
      {
         return(sHBonding[i].donor);
      }
   }

   if(allowPseudo)
   {
      for(i=0; sPseudoHBonding[i].element != NULL; i++)
      {
         if(!strcmp(p->element, sPseudoHBonding[i].element))
         {
            *pseudo = TRUE;
            return(sPseudoHBonding[i].donor);
         }
      }
   }
   else
   {
      return(0);
   }
   

   return(-1);
}

/************************************************************************/
/*>int IsAcceptor(PDB *p, BOOL allowPseudo, BOOL *pseudo)
   ------------------------------------------------------
*//**
   \param[in]   *p            PDB pointer
   \param[in]   allowPseudo   Allow pseudo HBonds?
   \param[out]  *pseudo       Is it a pseudo acceptor?
   \return                    Number of acceptor positions

   Tests whether an atom is a HBond (pseudo-)acceptor

-  07.06.99 Original   By: ACRM
-  22.07.15 Modified to use PDB files and standard BiopLib structures
            and functions
*/
int IsAcceptor(PDB *p, BOOL allowPseudo, BOOL *pseudo)
{
   int i;

   *pseudo = FALSE;

   for(i=0; sHBonding[i].element != NULL; i++)
   {
      if(!strcmp(p->element, sHBonding[i].element))
      {
         return(sHBonding[i].acceptor);
      }
   }

   if(allowPseudo)
   {
      for(i=0; sPseudoHBonding[i].element != NULL; i++)
      {
         if(!strcmp(p->element, sPseudoHBonding[i].element))
         {
            *pseudo = TRUE;
            return(sPseudoHBonding[i].acceptor);
         }
      }
   }
   else
   {
      return(0);
   }

   return(-1);
}


/************************************************************************/
/*>PDB *FindAntecedent(PDB *atom, PDB **pdbarray, int *count, int nth)
   -------------------------------------------------------------------
*//**
   \param[in]   *atom       The atom of interest
   \param[in]   **pdbarray  Array of PDB pointers indexed by atom number
   \param[in]   nth         The Nth antecedent is requested
   \param[out]  *count      The Nth antecedent found
   \return                  Pointer to Nth antecedent

   Finds the nth antecedent atom (nth==0 is the same as nth==1).
   In other words finds the nth atom that is bound to the atom. 
   count contains the actual antecedent number as there may not be 
   enough atoms bound.

   A complete rewrite of the old XMAS version

-  21.07.15  Original   By: ACRM
*/
PDB *FindAntecedent(PDB *atom, PDB **pdbarray, int *count, int nth)
{
   int i,
       atomnum;
   PDB *antecedent = NULL, 
       *best       = NULL,
       *p; 

   if(nth==0) nth=1;

   *count = 0;
   if(atom != NULL)
   {
      /* If it's a HETATM with CONECT information                       */
      if((atom->atomtype & ATOMTYPE_NONRESIDUE) || 
         (atom->atomtype == ATOMTYPE_MODPROT)   ||
         (atom->atomtype == ATOMTYPE_MODNUC))
      {
         if(atom->nConect == 0)
            return(NULL);

         if(nth > atom->nConect)
            nth = atom->nConect;

         antecedent = atom->conect[nth-1];
         *count = nth;
         return(antecedent);
      }
      else /* It's a normal residue                                     */
      {
         /* Work backwards a maximum of 30 atoms                        */
         atomnum = atom->atnum - 1;
         for(i=atomnum; i>=0 && i>atomnum-30; i--)
         {
            if(pdbarray[i] != NULL)
            {
               if(blIsBonded(atom, pdbarray[i], BOND_TOL))
               {
                  best = pdbarray[i];
                  if(++(*count) >= nth)
                     return(pdbarray[i]);
               }
            }
         }

         /* Work forwards a maximum of 30 atoms                         */
         for(p=atom->next, i=0; p!=NULL && i<30; NEXT(p))
         {
            if(blIsBonded(atom, p, BOND_TOL))
            {
               best = p;
               if(++(*count) >= nth)
                  return(p);
            }
            i++;
         }
      }
   }

   return(best);
}


/************************************************************************/
/*>PDB *FindBondedHydrogen(PDB *pdb, PDB *donor, PDB *acceptor)
   ------------------------------------------------------------
*//**
   \param[in]     *pdb      Start of PDB linked list
   \param[in]     *donor    Donor atom
   \param[in]     *acceptor Acceptor atom
   \return                  Hydrogen bonded to the donor

   Finds a hydrogen bonded to the donor. Selects the hydrogen which is
   closest to acceptor.

-  07.06.99 Original   By: ACRM
-  16.06.99 Initialise bestDistSq to 0
-  22.07.15 Modified to use PDB files and standard BiopLib structures
            and functions - added pdb parameter
*/
PDB *FindBondedHydrogen(PDB *pdb, PDB *donor, PDB *acceptor)
{
   PDB  *hydrogen = NULL,
        *p;
   PDB  *start, 
        *stop;
   REAL bestDistSq = (REAL)0.0,
        distSq;
   
   if(donor != NULL)
   {
      /* If the donor is an oxygen, or a sidechain Nitrogen on a lysine
         it is rotatable, so we don't know where the H is
      */
      if(!strcmp(donor->element, "O") ||
         (!strncmp(donor->atnam, "NZ  ", 4) &&
          !strncmp(donor->resnam, "LYS", 3)))
      {
         return(NULL);
      }
      
      /* Search forward for a hydrogen                                  */
      for(p=donor->next;
          p!=NULL && RESIDMATCH(p, donor);
          NEXT(p))
      {
         if(!strcmp(p->element, "H") &&
            (DISTSQ(p, donor) <= MAXBONDSQ))
         {
            if(hydrogen==NULL)
            {
               hydrogen = p;
               bestDistSq = DISTSQ(p, acceptor);
            }
            else
            {
               distSq = DISTSQ(p, acceptor);
               if(distSq < bestDistSq)
               {
                  bestDistSq = distSq;
                  hydrogen = p;
               }
            }
         }
      }
      
      /* Search back for a hydrogen                                     */
      if(!strncmp(donor->record_type, "HETATM", 6))
      {
         start = blFindHetatmResidue(pdb, donor->chain, donor->resnum,
                                     donor->insert);
      }
      else
      {
         start = blFindResidue(pdb, donor->chain, donor->resnum,
                               donor->insert);
      }

      if(start)
      {
         stop  = blFindNextResidue(start);
         for(p=start; p!=stop; NEXT(p))
         {
            if(!strcmp(p->element, "H") && 
               blIsBonded(p, acceptor, BOND_TOL))
            {
               if(hydrogen == NULL)
               {
                  hydrogen = p;
                  bestDistSq = DISTSQ(p, acceptor);
               }
               else
               {
                  distSq = DISTSQ(p, acceptor);
                  if(distSq < bestDistSq)
                  {
                     bestDistSq = distSq;
                     hydrogen = p;
                  }
               }
            }
         }
      }

      return(hydrogen);
   }
   
   return(NULL);
}


/************************************************************************/
/*>HBLIST *TestForHBond(PDB *pdb, PDB *p, PDB *q, PDB **pdbarray,
                        BOOL pseudo, REAL maxHBDistSq)
   -----------------------------------------------------------------------
*//**
   \param[in]       *pdb          Start of PDB linked list
   \param[in]       *p            Ligand atom
   \param[in]       *q            Protein atom
   \param[in]       **pdbarray    Array of PDB pointers indexed by atom
                                  number
   \param[in]       pseudo        Is this to be a pseudo HBond rather 
                                  than a real one?
   \param[in]       maxHBDistSq   Max allowed D-A HBond distance
   \return                        Allocated hbond structure or NULL

   Tests whether the 2 atoms are linked by a HBond. If found, returns
   a malloc'd HBLIST structure (maybe a list).

-  07.06.99 Original   By: ACRM
-  16.06.99 Added maxHBDistSq parameter
-  22.07.15 Modified to use PDB files and standard BiopLib structures
            and functions - added pdb parameter
*/
HBLIST *TestForHBond(PDB *pdb, PDB *p, PDB *q, PDB **pdbarray,
                     BOOL pseudo, REAL maxHBDistSq)
{
   PDB     *acceptor = NULL,
           *donor    = NULL;
   int     donMax,
           accMax;
   BOOL    pseudoA,
           pseudoD;
   HBLIST  *hb1      = NULL,
           *hb2      = NULL;

   /* Is one a donor and the other an acceptor?
      We only allow pseudo Hbonds where the donor is the ligand and the
      acceptor is the protein.
   */
   if(((donMax=IsDonor(p, pseudo, &pseudoD))!=0) && 
      ((accMax=IsAcceptor(q, pseudo, &pseudoA))!=0))
   {
      donor    = p;
      acceptor = q;

      /* If we are doing pseudo HBonds then one *must* be a pseudo.     */
      if((!pseudo) ||
         (pseudo && (pseudoD || pseudoA)))
      {
         hb1 = doTestForHBond(pdb, donor, acceptor, pdbarray, donMax,
                              maxHBDistSq);
      }
   }

   if(!pseudo)
   {
      if(((donMax=IsDonor(q, FALSE, &pseudoD))!=0) && 
         ((accMax=IsAcceptor(p, FALSE, &pseudoA))!=0))
      {
         donor    = q;
         acceptor = p;
         hb2 = doTestForHBond(pdb, donor, acceptor, pdbarray, donMax,
                              maxHBDistSq);
      }
   }

   /* If one of the HBond is NULL, just return the other (which may also
      be NULL)
   */
   if(hb1==NULL)
   {
      return(hb2);
   }
   else if(hb2==NULL)
   {
      return(hb1);
   }

   /* If both hb1 and hb2 are the same HBond, or the same atoms making
      a bond in the other direction, then delete one of them and return
      the other
   */
   if(((hb1->donor == hb2->donor) && (hb1->acceptor == hb2->acceptor)) ||
      ((hb1->donor == hb2->acceptor) && (hb1->acceptor == hb2->donor)))
   {
      free(hb2);
      return(hb1);
   }

   /* Otherwise, link them into a single list and return it             */
   hb1->next = hb2;
   return(hb1);
}


/************************************************************************/
/*>HBLIST *doTestForHBond(PDB *pdb, PDB *donor, PDB *acceptor, 
                          PDB **pdbarray, int donMax, REAL maxHBDistSq)
   --------------------------------------------------------------------
*//**
   \param[in]     *pdb        Start of PDB linked list
   \param[in]     *donor      Donor PDB pointer
   \param[in]     *acceptor   Acceptor PDB pointer
   \param[in]     **pdbarray  Array of PDB pointers indexed by atom
                              number
   \param[in]     donMax      Max number of donor connections
   \param[in]     maxHBDistSq Max allowed D-A distance for HBond
   \return                    Malloc'd HBond structure

   Does the actual work of testing for an HBond between a donor and
   acceptor and allocating a structure for it.

-  07.06.99 Original   By: ACRM
-  16.06.99 Added maxHBDistSq parameter
-  04.11.99 Modified the test for planarity to test two antecedent
            atoms rather than just one
-  12.05.00 Now if the planarity test fails, we just set relaxed to TRUE 
            in the HBLIST structure
-  21.07.15 Modified to use PDB files and standard BiopLib structures
            and functions - added pdb parameter. 
-  22.07.15 Added check that antecendent isn't the donor
*/
HBLIST *doTestForHBond(PDB *pdb, PDB *donor, PDB *acceptor, 
                       PDB **pdbarray, int donMax, REAL maxHBDistSq)
{
   HBLIST  *hb         = NULL;
   int     donCount    = 0,
           accCount    = 0;
   PDB     *antecedent = NULL,
           *hydrogen   = NULL;
   BOOL    relaxed     = FALSE;
   

   /* If we are using oxygen as a donor, then check if the antecedent is
      carbon and is trigonal planar. If this is the case then it will be
      a double bond to the oxygen which will not be able to donate a
      hydrogen. (It could also be a carboxyl group, but we assume these
      are charged - i.e. the hydrogen isn't present).

      12.05.00 Rather than returning NULL if the antedecent is planar,
      we now set relaxed to true and calculate anyway.
   */
   if((donor!=NULL) && !strcmp(donor->element, "O"))
   {
      PDB  *donAnt,
           *donAnt2;
      REAL ang;
      
      /* Find the antecedent for the donor                              */
      donAnt = FindAntecedent(donor,pdbarray,&accCount,0);
      if((donAnt!=NULL) && (accCount==1) && 
         !strcmp(donAnt->element, "C"))
      {
         /* Find first pre-antecendent                                  */
         donAnt2 = FindAntecedent(donAnt,pdbarray,&accCount,1);
         /* Added 22.07.15                                              */
         if(donAnt2 == donor)
         {
            donAnt2 = FindAntecedent(donAnt,pdbarray,&accCount,2);
         }
            
         if(donAnt2!=NULL)
         {
            REAL sang1;
            
            ang = blAngle(donor->x,   donor->y,   donor->z, 
                          donAnt->x,  donAnt->y,  donAnt->z, 
                          donAnt2->x, donAnt2->y, donAnt2->z);
            /* Test for planarity                                       */
            sang1 = blSimpleangle(ang);
            if(sang1 > MAXTETRAHEDRALANGLE*PI/180.0)
            {
               /* OK seems to be planar, test a second pre-antecedent   */
               donAnt2 = FindAntecedent(donAnt,pdbarray,
                                        &accCount,2);
               /* If NULL, then just believe that the first angle was 
                  really indicating trigonal planar
               */
               if(donAnt2 == NULL)
               {
                  relaxed = TRUE;   /* 12.05.00 ACRM - was return(NULL) */
               }
               else
               {
                  /* Check whether we just found the donor again in which
                     case we find the one before that
                  */
                  if(donAnt2 == donor)
                  {
                     donAnt2 = FindAntecedent(donAnt,pdbarray,
                                              &accCount,3);
                     
                     /* If we still found the donor then there weren't any
                        more atoms so believe the first angle was really
                        indicating trigonal planar
                     */
                     if((donAnt2 == donor) || (donAnt2 == NULL))
                     {
                        /* 12.05.00 ACRM - was return(NULL)             */
                        relaxed = TRUE;
                     }
                  }

                  /* 12.05.00 ACRM - New check on relaxed               */
                  if(!relaxed)
                  {
                     /* Find the angle to the second pre-antecedent     */
                     ang = blAngle(donor->x,   donor->y,   donor->z, 
                                   donAnt->x,  donAnt->y,  donAnt->z, 
                                   donAnt2->x, donAnt2->y, donAnt2->z);
                     
                     /* Test for planarity - the sum of the angles must be
                        greater than 2*MAXTETRAHEDRALANGLE (130deg)
                     */
                     if((sang1+blSimpleangle(ang)) > 
                        2*MAXTETRAHEDRALANGLE*PI/180.0)
                     {
                        /* 12.05.00 ACRM - was return(NULL)             */
                        relaxed = TRUE;
                     }
                  }
               }
            }
         }
      }
   }
   
   /* If we found an acceptor (and a donor)                             */
   if(acceptor != NULL)
   {
      /* If they are in distance range                                  */
      if(DISTSQ(donor, acceptor) <= maxHBDistSq)
      {
         /* Find the antecedent for the acceptor                        */
         antecedent = FindAntecedent(acceptor,pdbarray,
                                     &accCount,0);

         /* If the donor is in the protein we can see if we actually have
            a hydrogen here
         */
         if(donor->atomtype == ATOMTYPE_ATOM)
            hydrogen = FindBondedHydrogen(pdb, donor, acceptor);
         else
            hydrogen = NULL;

         donCount = donor->nConect;

         /* Check whether this atom is already bonded so not available
            for coordination
         */
         if((donMax!=(-1)) && (donCount >= donMax))
            return(NULL);

         if(blValidHBond(hydrogen, donor, acceptor, antecedent))
         {
            INIT(hb, HBLIST);
            if(hb == NULL)
               return(NULL);
            
            hb->next     = NULL;
            hb->donor    = donor;
            hb->acceptor = acceptor;
            hb->relaxed  = relaxed;                    /* 12.05.00 ACRM */
            return(hb);
         }
         else
         {
            return(NULL);
         }
      }
   }
   return(NULL);
}


/************************************************************************/
/*>HBLIST *FindNonBonds(PDB *pdb, PDB **pdbarray, HBLIST *hbonds,
                        REAL minNBDistSq, REAL maxNBDistSq)
   --------------------------------------------------------------
*//**
   \param[in]    *pdb         The PDB linked list
   \param[in]    **pdbarray   Array of PDB structure indexed by atom
                              number
   \param[in]    *hbonds      Linked list of HBonds
   \param[in]    minNBDistSq  Minimum distance for non-bond contact
   \param[in]    maxNBDistSq  Maximum distance for non-bond contact
   \return                    Linked list of non-bonds

   Finds non-bonded contacts between ligand and protein/nucleotide or
   between nucleotide and protein.

-  1. For atom pairs within a distance cutoff: 2.7 - 3.35A (centre-centre)
-  2. a) Not covalently bonded
      b) Not H-bonded
      c) Not in the same residue

-  07.06.99 Original   By: ACRM
-  16.06.99 Does NBond interactions for peptides
            Initialise nb to NULL
            Skips interactions of peptide with itself
            Skips interactions involving hydrogen
            min and max distances now variables (and parameters)
-  21.07.15 Modified to use PDB files and standard BiopLib structures
            and functions - added pdb parameter. 
*/
HBLIST *FindNonBonds(PDB *pdb, PDB **pdbarray, HBLIST *hbonds, 
                     REAL minNBDistSq, REAL maxNBDistSq)
{
   PDB    *p, 
          *q;
   REAL   distSq;
   HBLIST *nblist = NULL,
          *nb     = NULL;
   BOOL   isPeptide;
   

   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* Skip hydrogens                                                 */
      if(!strcmp(p->element, "H"))
         continue;
      
      isPeptide = isAPeptide(pdb, p);
      
      /* If it's a HET/METAL/BOUNDHET or a peptide                      */
      if(((p->atomtype & ATOMTYPE_NONRESIDUE) && 
          (p->atomtype != ATOMTYPE_WATER)) ||
         isPeptide)
      {
         for(q=pdb; q!=NULL; NEXT(q))
         {
            /* Skip hydrogens                                           */
            if(!strcmp(q->element, "H"))
               continue;

            /* If our first molecule is a peptide, then skip if this is
               the same molecule
            */
            if(isPeptide && 
               (PDBEXTRASPTR(p, PDBEXTRAS)->molid == 
                PDBEXTRASPTR(q, PDBEXTRAS)->molid))
               continue;

            if(p==q)
               continue;

            /* If it's a protein/nucleotide                             */
            if(!(q->atomtype & ATOMTYPE_NONRESIDUE) &&
                (q->atomtype != ATOMTYPE_UNDEF))
            {
               distSq = DISTSQ(p,q);
               if(distSq >= minNBDistSq && distSq <= maxNBDistSq)
               {
                  if(!RESIDMATCH(p, q)  &&
                     !blIsConected(p, q) &&
                     !IsListedAsHBonded(p, q, hbonds))
                  {
                     if(nblist==NULL)
                     {
                        INIT(nblist, HBLIST);
                        nb = nblist;
                     }
                     else
                     {
                        ALLOCNEXT(nb, HBLIST);
                     }
                     if(nb==NULL)
                     {
                        FREELIST(nblist, HBLIST);
                        fprintf(stderr,"pdbhbond: (error) No memory for \
Non-bond list\n");
                        return(NULL);
                     }
                     nb->donor    = p;
                     nb->acceptor = q;
                  }
               }
            }
         }
      }
      else if((p->atomtype == ATOMTYPE_NUC) || 
              (p->atomtype == ATOMTYPE_MODNUC))
      {
         /* If it's a nucleotide
            Look for interactions with protein
         */
         for(q=pdb; q!=NULL; NEXT(q))
         {
            if(p==q)
               continue;

            if((q->atomtype == ATOMTYPE_ATOM) ||
               (q->atomtype == ATOMTYPE_MODPROT) ||
               (q->atomtype == ATOMTYPE_NONSTDAA))
            {
               distSq = DISTSQ(p,q);
               if(distSq >= minNBDistSq && distSq <= maxNBDistSq)
               {
                  if(!RESIDMATCH(p, q) &&
                     !blIsConected(p, q) &&
                     !IsListedAsHBonded(p, q, hbonds))
                  {
                     if(nblist==NULL)
                     {
                        INIT(nblist, HBLIST);
                        nb = nblist;
                     }
                     else
                     {
                        ALLOCNEXT(nb, HBLIST);
                     }
                     if(nb==NULL)
                     {
                        FREELIST(nblist, HBLIST);
                        fprintf(stderr,"pdbhbond: (error) No memory for \
Non-bond list\n");
                        return(NULL);
                     }
                     nb->donor    = p;
                     nb->acceptor = q;
                  }
               }
            }
         }
      }
   }
   return(nblist);
}


/************************************************************************/
/*>BOOL IsListedAsHBonded(PDB *p, PDB *q, HBLIST *hbonds)
   ------------------------------------------------------
*//**
   \param[in]     *p      PDB pointer
   \param[in]     *q      PDB pointer
   \param[in]     *hbonds List of HBonds iedntified thus far
   \return                Listed?

   Tests whether the two specified atoms are already listed as being
   hydrogen bonded

-  07.06.99 Original   By: ACRM
-  21.07.15 Modified to use PDB files and standard BiopLib structures
            and functions - added pdb parameter. 
*/
BOOL IsListedAsHBonded(PDB *p, PDB *q, HBLIST *hbonds)
{
   HBLIST *h;

   for(h=hbonds; h!=NULL; NEXT(h))
   {
      if(((h->donor == p) && (h->acceptor == q)) ||
         ((h->donor == q) && (h->acceptor == p)))
      {
         return(TRUE);
      }
   }

   return(FALSE);
}


/************************************************************************/
/*>HBLIST *FindLigandLigandHBonds(PDB *pdb, PDB **pdbarray, BOOL pseudo,
                                  REAL maxHBDistSq)
   ---------------------------------------------------------------------
*//**
   \param[in]      *pdb        PDB linked list
   \param[in]      **pdbarray  Array of PDB pointers indexed by atom
                               number
   \param[in]      pseudo      Pseudo hbonds? (or true HBonds)
   \param[in]      maxHBDistSq Max D-A Hbond distance
   \return                     Linked list of hbonds

   Finds HBonds between ligands. If pseudo is true then it finds 
   pseudo-HBonds rather than real ones. Basically a copy of the first
   part of FindProtLigandHBonds()

-  08.06.99 Original   By: ACRM
-  16.06.99 Added maxHBDistSq parameter
-  21.07.15 Modified to use PDB files and standard BiopLib structures
            and functions - added pdb parameter. 
*/
HBLIST *FindLigandLigandHBonds(PDB *pdb, PDB **pdbarray, BOOL pseudo,
                               REAL maxHBDistSq)
{
   PDB           *p, *q;
   static HBLIST *hblist = NULL,
                 *hbl,
                 *hb;
   static int    prevPseudo = (-5);

   
   if(pseudo != prevPseudo)
   {
      prevPseudo = pseudo;
      hblist = NULL;
      hbl    = NULL;
      hb     = NULL;
   }
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If it's a HET/METAL/BOUNDHET                                   */
      if((p->atomtype & ATOMTYPE_NONRESIDUE) &&
         (p->atomtype != ATOMTYPE_WATER))
      {
         /* Look for interactions with protein or nucleotide            */
         for(q=pdb; q!=NULL; NEXT(q))
         {
            /* Inter-molecule only...                                   */
            if((p == q) || PDBCHAINMATCH(p, q))
               continue;
            
            /* If it's a HET/METAL/BOUNDHET                             */
            if((q->atomtype & ATOMTYPE_NONRESIDUE) && 
               (q->atomtype != ATOMTYPE_WATER))
            {
               /* If they are not already covalently bonded or in the
                  current HBond list
               */
               if(!blIsConected(p, q) &&
                  !IsListedAsHBonded(p, q, hblist))
               {
                  if((hb=TestForHBond(pdb, p, q,pdbarray,pseudo,
                                      maxHBDistSq))!=NULL)
                  {
                     /* Store the Hbond                                 */
                     if(hblist==NULL)
                     {
                        hblist = hbl = hb;
                     }
                     else
                     {
                        hbl->next = hb;
                     }
                     if(hbl!=NULL)
                        LAST(hbl);
                  }
               }
            }
         }
      }
   }
   
   return(hblist);
}


/************************************************************************/
/*>HBLIST *FindProtLigandHBonds(PDB *pdb, PDB **pdbarray, BOOL pseudo,
                                REAL maxHBDistSq)
   -------------------------------------------------------------------
*//**
   \param[in]     *pdb        PDB linked list
   \param[in]     **pdbarray  Array of PDB pointers indexed by atom
                              number
   \param[in]     pseudo      Pseudo hbonds? (or true HBonds)
   \param[in]     maxHBDistSq Max D-A HBond distance
   \return                    Linked list of hbonds

   Finds HBonds between protein and ligand. If pseudo is true then it
   finds pseudo-HBonds rather than real ones.

-  07.06.99 Original   By: ACRM
-  08.06.99 Modified to treat peptides as ligands and not to repeat
            HBonds already in the current list. Note that peptide/
            protein HBonds will also appear in the main protein/protein
            HBonds list
-  16.06.99 Added maxHBDistSq parameter
-  03.11.99 Added check that the lignad and protein are different
            molecules!
-  21.07.15 Modified to use PDB files and standard BiopLib structures
            and functions - added pdb parameter. 
*/
HBLIST *FindProtLigandHBonds(PDB *pdb, PDB **pdbarray, BOOL pseudo,
                             REAL maxHBDistSq)
{
   PDB           *p, *q;
   static HBLIST *hblist = NULL,
                 *hbl,
                 *hb;
   static int    prevPseudo = (-5);

   
   if(pseudo != prevPseudo)
   {
      prevPseudo = pseudo;
      hblist = NULL;
      hbl    = NULL;
      hb     = NULL;
   }
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If it's a HET/METAL/BOUNDHET                                   */
      if((p->atomtype & ATOMTYPE_NONRESIDUE) && 
         (p->atomtype != ATOMTYPE_WATER))
      {
         /* Look for interactions with protein or nucleotide            */
         for(q=pdb; q!=NULL; NEXT(q))
         {
            /* 03.11.99 Check that it's a different molecule as well as 
               a different atom
            */
            if((p == q) ||
               (PDBEXTRASPTR(p, PDBEXTRAS)->molid ==
                PDBEXTRASPTR(q, PDBEXTRAS)->molid))
               continue;
            
            /* If it's a protein/nucleotide                             */
            if(!(q->atomtype & ATOMTYPE_NONRESIDUE) &&
               (q->atomtype != ATOMTYPE_UNDEF))
            {
               /* If they are not already covalently bonded or in the
                  current HBond list
               */
               if(!blIsConected(p, q) &&
                  !IsListedAsHBonded(p, q, hblist))
               {
                  if((hb=TestForHBond(pdb, p,q,pdbarray,pseudo,
                                      maxHBDistSq))!=NULL)
                  {
                     /* Store the Hbond                                 */
                     if(hblist==NULL)
                     {
                        hblist = hbl = hb;
                     }
                     else
                     {
                        hbl->next = hb;
                     }
                     if(hbl!=NULL)
                        LAST(hbl);
                  }
               }
            }
         }
      }
      else if(!pseudo)
      {
         /* If it's a nucleotide or a peptide                           */
         if((p->atomtype == ATOMTYPE_NUC) || 
            (p->atomtype == ATOMTYPE_MODNUC) ||
            isAPeptide(pdb, p))
         {
            /* Look for interactions with protein                       */
            for(q=pdb; q!=NULL; NEXT(q))
            {
               /* 03.11.99 Check that it's a different molecule as well
                  as a different atom
               */
               if((p == q) || 
                  (PDBEXTRASPTR(p, PDBEXTRAS)->molid == 
                   PDBEXTRASPTR(q, PDBEXTRAS)->molid))
                  continue;

               if((q->atomtype == ATOMTYPE_ATOM) ||
                  (q->atomtype == ATOMTYPE_MODPROT) ||
                  (q->atomtype == ATOMTYPE_NONSTDAA))
               {
                  /* If they are not already covalently bonded or in the
                     current HBond list
                  */
                  if(!blIsConected(p, q) &&
                     !IsListedAsHBonded(p, q, hblist))
                  {
                     if((hb=TestForHBond(pdb,p,q,pdbarray,pseudo,
                                         maxHBDistSq))
                        != NULL)
                     {
                        /* Store the Hbond                              */
                        if(hblist==NULL)
                        {
                           hblist = hbl = hb;
                        }
                        else
                        {
                           hbl->next = hb;
                        }
                        if(hbl!=NULL)
                           LAST(hbl);
                     }
                  }
               }
            }
         }
      }
   }
   
   return(hblist);
}


/************************************************************************/
/*>BOOL isAPeptide(PDB *pdb, PDB *atm)
   -----------------------------------
*//**
   \param[in]    pdb    Start of PDB linked list
   \param[in]    atm    An atom in the PDB linked list
   \return              Is it a peptide

   Tests whether the molecule containing the given atom is a peptide
   (i.e. less than MAX_PEPTIDE_LENGTH residues long)

-  21.07.15  Original   By: ACRM
*/
BOOL isAPeptide(PDB *pdb, PDB *atm)
{
   PDB *start, *stop, *res;
   int nRes;
   
   
   /* Find the start of this chain                                      */
   for(start=pdb; start!=atm; NEXT(start))
   {
      if(PDBCHAINMATCH(start, atm))
         break;
   }
   /* Find the next chain                                               */
   stop = blFindNextChain(start);
   
   /* Count the residues in the chain                                   */
   for(res=start, nRes=0; res!=stop; res=blFindNextResidue(res))
   {
      nRes++;
   }

   if(nRes > MAX_PEPTIDE_LENGTH)
      return(FALSE);
   
   return(TRUE);
}


/************************************************************************/
/*>void SetAtomNumExtras(PDB *pdb)
   -------------------------------
*//**
   \param[in,out]   PDB  *pdb   PDB Linked list

   Walks the PDB linked list storing the original atom numbers in a 
   field of the PDB.extras structure.

-  21.07.15  Original   By: ACRM
*/
void SetAtomNumExtras(PDB *pdb)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      PDBEXTRASPTR(p, PDBEXTRAS)->origAtnum = p->atnum;
   }
}

/************************************************************************/
/*>BOOL UpdatePDBExtras(PDB *pdb)
   ------------------------------
*//**
   \param[in, out]   *pdb    PDB linked list
   \return                   Success in allocations

   Walks the PDB linked list creating an 'extra' structure for each
   PDB entry that doesn't already have one. 
   It initializes the PDB.extras.origAtnum to -1
   and the PDB.extras.molid to 0

-  21.07.15  Original   By: ACRM
*/
BOOL UpdatePDBExtras(PDB *pdb)
{
   PDB *p;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->extras == NULL)
      {
         if((p->extras = (APTR)malloc(sizeof(PDBEXTRAS)))==NULL)
            return(FALSE);
         PDBEXTRASPTR(p, PDBEXTRAS)->origAtnum = (-1);
         PDBEXTRASPTR(p, PDBEXTRAS)->molid     = 0;
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>BOOL SetMolecules(PDB *pdb)
   ---------------------------
*//**
   \param[in]   *pdb         PDB linked list
   \return                   Success?

   Identifies all individual molecules within the structure.

-  23.03.99 Original   By: ACRM
-  01.04.99 Only checks for peptide if it's already a protein and
            now also checks for CA only
-  21.05.99 Makes sure that the molecule ID gets incremented for each
            protein chain as well as for HETATMs. Also labels all
            atoms within a protein chain with the molecule ID.
-  21.07.15 Heavily modified to deal with being passed PDB linked list
            No longer stores a list of molecules - just annotates them
            in the PDB.extras.molid field
*/
BOOL SetMolecules(PDB *pdb)
{
   PDB  *chainStart,
        *nextChain,
        *resStart,
        *nextRes,
        *firstAtom,
        *p;
   BOOL GotAtoms;
   int  id=0;
   

   /* Clear all the molid flags                                         */
   for(p=pdb; p!=NULL; NEXT(p))
      PDBEXTRASPTR(p, PDBEXTRAS)->molid = 0;
   
   /* For each chain......                                              */
   for(chainStart=pdb; chainStart!=NULL; chainStart=nextChain)
   {
      nextChain = blFindNextChain(chainStart);

      /* See if we have ATOM records in this chain                      */
      GotAtoms = FALSE;
      for(firstAtom=chainStart; firstAtom!=nextChain; NEXT(firstAtom))
      {
         /* If the high bit isn't set (i.e. it's not a ligand)          */
         if(!(firstAtom->atomtype & ATOMTYPE_NONRESIDUE))
         {
            GotAtoms = TRUE;
            break;
         }
      }

      /* If we do have atom records, then write a molecule entry for the
         chain, checking if it's a peptide rather than a protein chain
      */
      if(GotAtoms)
      {
         /* Set the molecule ID for all atoms in this chain             */
         id++;
         
         for(p=chainStart; p!=nextChain; NEXT(p))
         {
            /* If the high bit isn't set (i.e. it's not a ligand)       */
            if(!(p->atomtype & ATOMTYPE_NONRESIDUE))
            {
               PDBEXTRASPTR(p, PDBEXTRAS)->molid = id;
            }
         }

      }

      /* Work through the chain again, a residue at a time, looking for 
         HETATMs

         For each residue in this chain....
      */
      for(resStart=chainStart; resStart!=nextChain; resStart=nextRes)
      {
         /* 18.06.99 Changed from FindNextIResidue()                    */
         nextRes = blFindNextResidue(resStart);

         /* If it is a het residue                                      */
         if((resStart->atomtype==ATOMTYPE_HETATM)
            || (resStart->atomtype==ATOMTYPE_METAL)
            || (resStart->atomtype==ATOMTYPE_BOUNDHET)
#ifndef NO_WATERS_IN_MOL_LIST
            || (resStart->atomtype==ATOMTYPE_WATER)
#endif
            )
         {
            /* If we haven't set the flag to say this residue has been
               used, then we create a new molecule entry
            */
            if(PDBEXTRASPTR(resStart, PDBEXTRAS)->molid == 0)
            {
               /* Mark all het residues which are linked to this one    */
               id++;
               MarkLinkedResidues(chainStart, resStart, nextChain, 
                                  id);
            }
         }
      }
   }

   return(TRUE);
}


/************************************************************************/
/*>void MarkLinkedResidues(PDB *chainStart, PDB *resStart,
                           PDB *nextChain, int id)
   -------------------------------------------------------
*//**
   \param[in]       *head       Head of PDB structure
   \param[in]       *chainStart Start of this chain
   \param[in]       *resStart   First residue of maybe-linked residues
   \param[in]       *nextChain  Start of next chain
   \param[in]       id          Indentifier for this group
   \param[out]      *m          Molecule structure

   Marks all HET residues linked by CONNECTs to resStart. When a link is
   found, is called recursively to mark HET residues connected to that
   one. Sets m->polymer if any such links were found.

-  23.03.99 Original   By: ACRM
-  21.07.15 Modified to use PDB files and standard BiopLib structures
            and functions - added pdb parameter. 
*/
void MarkLinkedResidues(PDB *chainStart, PDB *resStart, 
                        PDB *nextChain, int id)
{ 
   PDB *p, *q,
         *nextRes,
         *resStart2,
         *nextRes2;

   nextRes = blFindNextResidue(resStart);
   
   for(p=resStart; p!=nextRes; NEXT(p))
   {
      PDBEXTRASPTR(p, PDBEXTRAS)->molid = id;
      
      /* For each following residue....                                 */
      for(resStart2=nextRes; 
          resStart2!=nextChain; 
          resStart2=nextRes2)
      {
         int res2Molid;
         
         nextRes2 = blFindNextResidue(resStart2);

         /* If it's not already marked as linked and is a HET residue   */
         res2Molid = PDBEXTRASPTR(resStart2, PDBEXTRAS)->molid;

         if((res2Molid == 0) &&
            ((resStart2->atomtype==ATOMTYPE_HETATM)
#ifdef USE_METAL_IN_POLYHET
             || (resStart2->atomtype==ATOMTYPE_METAL)
#endif
             || (resStart2->atomtype==ATOMTYPE_BOUNDHET)
#ifndef NO_WATERS_IN_MOL_LIST
             || (resStart2->atomtype==ATOMTYPE_WATER)
#endif
           ))
         {
            /* For each atom in this other residue....                  */
            for(q=resStart2; q!=nextRes2; NEXT(q))
            {
               /* If they are connected, set this molecule as
                  a polymer and mark this residue as used
               */
               if(blIsConected(p, q))
               {
                  PDBEXTRASPTR(resStart2, PDBEXTRAS)->molid = id;
                  MarkLinkedResidues(chainStart, resStart2, 
                                     nextChain, id);
                  break;
               }
            }
         }
      }

      /* For each preceeding residue....                                */
      for(resStart2=chainStart; 
          resStart2!=resStart; 
          resStart2=nextRes2)
      {
         int res2Molid;
         
         nextRes2 = blFindNextResidue(resStart2);
                     
         /* If it's not already marked as linked and is a HET residue   */
         res2Molid = PDBEXTRASPTR(resStart2, PDBEXTRAS)->molid;
         if((res2Molid == 0) &&
            ((resStart2->atomtype==ATOMTYPE_HETATM)
#ifdef USE_METAL_IN_POLYHET
             || (resStart2->atomtype==ATOMTYPE_METAL)
#endif
             || (resStart2->atomtype==ATOMTYPE_BOUNDHET)
#ifndef NO_WATERS_IN_MOL_LIST
             || (resStart2->atomtype==ATOMTYPE_WATER)
#endif
           ))
         {
            /* For each atom in this other residue....                  */
            for(q=resStart2; q!=nextRes2; NEXT(q))
            {
               /* If they are connected, set this molecule as
                  a polymer and mark this residue as used
               */
               if(blIsConected(p, q))
               {
                  PDBEXTRASPTR(resStart2, PDBEXTRAS)->molid = id;
                  MarkLinkedResidues(chainStart, resStart2, 
                                     nextChain, id);
                  break;
               }
            }
         }
      }
   }
}


/************************************************************************/
/*>void DeleteMetalConects(PDB *pdb)
   ---------------------------------
*//**
   \param[in,out]   *pdb    PDB linked list

   Removes CONECT data related to metals. This is needed so that they
   can be processed as pseudo-H-bonds.

-  21.07.15  Original   By: ACRM
*/
void DeleteMetalConects(PDB *pdb)
{
   PDB *p;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->atomtype == ATOMTYPE_METAL)
      {
         blDeleteAtomConects(p);
      }
   }
}

