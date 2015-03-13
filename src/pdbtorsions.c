/************************************************************************/
/**

   \file       pdbtorsions.c
   
   \version    V2.1
   \date       04.03.15
   \brief      Calculate torsion angles for a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin 1996-2015
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
   Calculate torsion angles from a PDB file

**************************************************************************

   Usage:
   ======
   NOTE! If the executable is called 'torsions' rather than 'pdbtorsions'
   this has the effect of setting the '-o' (old style) flag by default.
   Consequently a symbolic link to the program can be used to obtain
   old-style output for backwards compatibility

**************************************************************************

   Revision History:
   =================

-  V1.0  10.06.94 Original
-  V1.1  16.08.94 Added -m option for Martin's output format
-  V1.2  12.06.95 Added -c option for pseudo-CA torsions     
-  V1.3  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.4  19.08.14 Added AsCopy suffix to call to blSelectAtomsPDB() 
                  By: CTP
-  V1.5  06.11.14 Renamed from torsions
-  V1.6  07.11.14 Initialized a variable
-  V2.0  27.11.14 Major rewrite to deal properly with multiple chains
                  and to associate the omega angle with the correct
                  residue. Still makes the old format available
-  V2.1  04.03.15 Improved checking for old name
                  Now done by blCheckProgName()

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/general.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/angle.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF     512
#define ERROR_VALUE 9999.0

/* Macro to set pointers 0..2 in an array to NULL                       */
#define CLEARVALUES(c)  \
do {                    \
   int i;               \
   for(i=0; i<3; i++)   \
   (c)[i] = NULL;       \
} while(0)

/* Macro to shift values in a 0..2 array back one position              */
#define UPDATEVALUES(c) \
do {                    \
   (c)[0] = (c)[1];     \
   (c)[1] = (c)[2];     \
   (c)[2] = NULL;       \
}  while(0)

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *CATorsions, BOOL *terse, BOOL *Radians, 
                  BOOL *oldStyle);
BOOL CalculateAndDisplayTorsions(FILE *out, PDB *fullpdb, 
                                 BOOL CATorsions, BOOL terse, 
                                 BOOL Radians, BOOL oldStyle);
void doCATorsions(FILE *out, PDB *pdb, BOOL terse, BOOL Radians, 
                  BOOL oldStyle);
void doFullTorsions(FILE *out, PDB *pdb, BOOL terse, BOOL Radians, 
                    BOOL oldStyle);
void PrintCARecord(FILE *out, PDB *p, REAL tor, BOOL terse, 
                   BOOL showLabel, BOOL dummy);
BOOL SetOldStyle(char *progname);
REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians);
void PrintFullRecord(FILE *out, PDB *p, REAL phi, REAL psi, REAL omega, 
                     BOOL terse, BOOL showLabel);
void BuildLabel(char *label, PDB *p, int width, BOOL LeftJustify);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for converting a PDB file to torsions.

-  10.06.94 Original   By: ACRM
-  16.08.94 Added -m option
-  12.06.95 Added -c option
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  19.08.14 Added AsCopy suffix to call to blSelectAtomsPDB() By: CTP
-  07.11.14 Initialized TorNum
*/
int main(int argc, char **argv)
{
   FILE    *in  = stdin,
           *out = stdout;
   char    inFile[MAXBUFF],
           outFile[MAXBUFF];
   int     natoms;
   PDB     *pdb;
   BOOL    CATorsions = FALSE;
   BOOL    terse      = FALSE;
   BOOL    Radians    = FALSE;
   BOOL    oldStyle   = FALSE;

   /* Set the default output style based on whether the program is called
      pdbtorsions or torsions
   */
   oldStyle = blCheckProgName(argv[0], "torsions");

   if(ParseCmdLine(argc, argv, inFile, outFile, 
                   &CATorsions, &terse, &Radians, &oldStyle))
   {
      if(blOpenStdFiles(inFile, outFile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natoms))!=NULL)
         {
            if(!CalculateAndDisplayTorsions(out, pdb, CATorsions, terse, 
                                            Radians, oldStyle))
               return(1);
         }
         else
         {
            fprintf(stderr,"pdbtorsions: Error - no atoms read from PDB \
file\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr,"pdbtorsions: Error - unable to open input or \
output file\n");
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
/*>BOOL CalculateAndDisplayTorsions(FILE *out, PDB *fullpdb, 
                                    BOOL CATorsions, BOOL terse, 
                                    BOOL Radians, BOOL oldStyle)
   -----------------------------------------------------------------------
*//**

   \param[in]    *out          Output file pointer
   \param[in]    *fullpdb      Full PDB linked list
   \param[in]    CATorsions    Do CA pseudo-torsions
   \param[in]    terse         Terse (single letter code AAs) output
   \param[in]    Radians       Use radians instead of degrees
   \param[in]    oldStyle      Old style output

   Calculate and display the torsion angles as required

- 27.11.14 Original   By: ACRM
*/
BOOL CalculateAndDisplayTorsions(FILE *out, PDB *fullpdb, BOOL CATorsions,
                                 BOOL terse, BOOL Radians, BOOL oldStyle)
{
   char *sel[4];
   int  natoms;
   PDB  *pdb;

   /* Set up the atom selection and select them                         */
   SELECT(sel[0],"CA  ");
   if(!CATorsions)
   {
      SELECT(sel[1],"N   ");
      SELECT(sel[2],"C   ");
   }

   if((pdb = blSelectAtomsPDBAsCopy(fullpdb,(CATorsions?1:3),sel,&natoms))
      == NULL)
   {
      fprintf(stderr,"pdbtorsions: Error - Unable to select backbone \
atoms from PDB file (no memory?)\n");
      return(FALSE);
   }

   if(CATorsions)
   {
      doCATorsions(out, pdb, terse, Radians, oldStyle);
   }
   else
   {
      doFullTorsions(out, pdb, terse, Radians, oldStyle);
   }

   return(TRUE);
}


/************************************************************************/
/*>void doCATorsions(FILE *out, PDB *pdb, BOOL terse, BOOL Radians, 
                     BOOL oldStyle)
   -----------------------------------------------------------------
*//**
   \param[in]    *out      Output file pointer
   \param[in]    *pdb      PDB linked list
   \param[in]    terse     Terse output
   \param[in]    Radians   Radians in output instead of degrees
   \param[in]    oldStyle  Old style output

   Main routine for doing CA pseudo-torsions

- 27.11.14 Original   By: ACRM
*/
void doCATorsions(FILE *out, PDB *pdb, BOOL terse, BOOL Radians, 
                  BOOL oldStyle)
{
   PDB  *p, 
        *start,
        *stop,
        *p1 = NULL,
        *p2 = NULL,
        *p3 = NULL,
        *p4 = NULL;
   REAL tor;
   
   if(oldStyle)
   {
      terse = TRUE;
   }

   /* Print title                                                       */
   if(oldStyle)
   {
      fprintf(out,"Res_N    CA_N--CA_(N+1)\n");
      fprintf(out,"--------------------------------------\n");
   }
   else
   {
      fprintf(out,"#ResnumN   ResnamN Torsion((N-1)--N)\n");
      fprintf(out,"#-----------------------------------\n");
   }
   
   /* Step through the chains                                           */
   for(start=pdb; start!=NULL; start=stop)
   {
      p1 = NULL;
      p2 = NULL;
      p3 = NULL;
      p4 = NULL;

      stop = blFindNextChain(start);
      
      for(p=start; p!=stop; NEXT(p))
      {
         p1=p2;
         p2=p3;
         p3=p4;
         p4=p;
         
         if(p1 && p2 && p3 && p4)   /* Got all 4 atoms                  */
         {
            PDB *key = p3;
            if(oldStyle)
               key = p2;
            
            tor = CalcTorsion(p1,p2,p3,p4,Radians);
            
            PrintCARecord(out, key, tor, terse, !oldStyle, FALSE);
         }
         else if(p1==NULL && p2==NULL && (!oldStyle || (p3==NULL)))
         {  /* Got 1 atom in old style or 2 atoms new style             */
            PrintCARecord(out, p4, tor, terse, !oldStyle, TRUE);
         }
      }
      if(oldStyle)
         PrintCARecord(out, p3, tor, terse, !oldStyle, TRUE);
      PrintCARecord(out, p4, tor, terse, !oldStyle, TRUE);
   }
}


/************************************************************************/
/*>void PrintCARecord(FILE *out, PDB *p, REAL tor, BOOL terse, 
                      BOOL showLabel, BOOL dummy)
   ------------------------------------------------------------
*//**

   \param[in]    *out        Output file pointer
   \param[in]    *p          Pointer to amino acid of interest
   \param[in]    tor         CA pseudo-torsion angle to print
   \param[in]    terse       Terse output style
   \param[in]    showLabel   Whether to show the residue label
   \param[in]    dummy       This is a dummy record (no torsion)

   Does the work of printing a CA pseudo-torsion in the required format

- 27.11.14 Original   By: ACRM
*/
void PrintCARecord(FILE *out, PDB *p, REAL tor, BOOL terse, 
                   BOOL showLabel, BOOL dummy)
{
   char resnam[16];
   char label[16];

   if(terse)
   {
      resnam[0] = blThrone(p->resnam);
      resnam[1] = '\0';
   }
   else
   {
      strcpy(resnam, p->resnam);
   }

   if(showLabel)
      BuildLabel(label, p, 8, TRUE);
   else
      label[0] = '\0';

   if(dummy)
   {
      fprintf(out,"%s   %s        -\n",label, resnam);
   }
   else
   {
      fprintf(out,"%s   %s    %8.3f\n",label, resnam,tor);
   }
}


/************************************************************************/
/*>void doFullTorsions(FILE *out, PDB *pdb, BOOL terse, BOOL Radians, 
                       BOOL oldStyle)
   -------------------------------------------------------------------
*//**
   \param[in]    *out       Output file pointer
   \param[in]    *pdb       PDB linked list
   \param[in]    terse      Terse output
   \param[in]    Radians    Radians in output instead of degrees
   \param[in]    oldStyle   Old style output

   Main routine for doing normal full torsion angles

- 27.11.14 Original   By: ACRM
*/
void doFullTorsions(FILE *out, PDB *pdb, BOOL terse, BOOL Radians, 
                    BOOL oldStyle)
{
   PDB  *startChain,
        *stopChain,
        *startRes,
        *stopRes,
        *N[3],
        *CA[3],
        *C[3];
   REAL phi, psi, omega, omega2;


   /* Print title                                                       */
   if(oldStyle)
   {
      fprintf(out,"               PHI      PSI     OMEGA\n");
      fprintf(out,"--------------------------------------\n");
   }
   else
   {
      fprintf(out,"#Resnum  Resnam     PHI      PSI     OMEGA\n");
      fprintf(out,"#------------------------------------------\n");
   }

   /* Step through the chains                                           */
   for(startChain=pdb; startChain!=NULL; startChain=stopChain)
   {
      /* Set all atom pointers to NULL                                  */
      CLEARVALUES(N);
      CLEARVALUES(CA);
      CLEARVALUES(C);

      stopChain = blFindNextChain(startChain);

      /* Step through the residues                                      */
      for(startRes=startChain; startRes!=stopChain; startRes=stopRes)
      {
         stopRes = blFindNextResidue(startRes);

         /* Shift atom pointers back one place                          */
         UPDATEVALUES(N);
         UPDATEVALUES(CA);
         UPDATEVALUES(C);
         
         /* Find the atoms of interest                                  */
         N[2]  = blFindAtomInRes(startRes, "N   ");
         CA[2] = blFindAtomInRes(startRes, "CA  ");
         C[2]  = blFindAtomInRes(startRes, "C   ");
         
         /* Calculate the torsions                                      */
         omega  = CalcTorsion(CA[0], C[0],  N[1],  CA[1], Radians);
         phi    = CalcTorsion(C[0],  N[1],  CA[1], C[1],  Radians);
         psi    = CalcTorsion(N[1],  CA[1], C[1],  N[2],  Radians);
         omega2 = CalcTorsion(CA[1], C[1],  N[2],  CA[2], Radians);

         if(oldStyle)
            PrintFullRecord(out, N[1], phi, psi, omega2, terse, oldStyle);
         else
            PrintFullRecord(out, N[1], phi, psi, omega, terse, oldStyle);
      }

      /* Deal with the last amino acid                                  */
      omega  = CalcTorsion(CA[1], C[1],  N[2],  CA[2], Radians);
      phi    = CalcTorsion(C[1],  N[2],  CA[2], C[2],  Radians);
      psi    = ERROR_VALUE;
      omega2 = ERROR_VALUE;

      if(oldStyle)
         PrintFullRecord(out, N[2], phi, psi, omega2, terse, oldStyle);
      else
         PrintFullRecord(out, N[2], phi, psi, omega, terse, oldStyle);
   }
}


/************************************************************************/
/*>void PrintFullRecord(FILE *out, PDB *p, REAL phi, REAL psi, 
                        REAL omega, BOOL terse, BOOL oldStyle)
   -----------------------------------------------------------
*//**
   \param[in]    *out      Output file pointer
   \param[in]    *p        Pointer to PDB record
   \param[in]    phi       Phi angle to print
   \param[in]    psi       Psi angle to print
   \param[in]    omega     Omega angle to print
   \param[in]    terse     Terse output
   \param[in]    oldStyle  Old style output

   Does the work of printing a record for a normal full torsion angle

- 27.11.14 Original   By: ACRM
*/
void PrintFullRecord(FILE *out, PDB *p, REAL phi, REAL psi, REAL omega, 
                     BOOL terse, BOOL oldStyle)
{
   char label[16];
   char resnam[16];

   if(p!=NULL)
   {
      if(terse)
      {
         resnam[0] = blThrone(p->resnam);
         resnam[1] = '\0';
      }
      else
      {
         strcpy(resnam, p->resnam);
      }

      if(oldStyle)
      {
         fprintf(out, "%5d%c %-4s %8.3f %8.3f %8.3f\n", 
                 p->resnum, p->insert[0], resnam, phi, psi, omega);
      }
      else
      {
         BuildLabel(label, p, 6, FALSE);
         fprintf(out, "%-8s %-4s    %8.3f %8.3f %8.3f\n", 
                 label, resnam, phi, psi, omega);
      }
   }
}


/************************************************************************/
/*>void BuildLabel(char *label, PDB *p, int width, BOOL LeftJustify)
   -----------------------------------------------------------------
*//**
   \param[out]   *label      Residue label
   \param[in]    *p          Pointer to PDB record
   \param[in]    width       Width for left-justified labels
   \param[in]    LeftJustify Should we left-justify

   Builds a label by concatenating chain label, residue number and insert
   code

- 27.11.14 Original   By: ACRM
*/
void BuildLabel(char *label, PDB *p, int width, BOOL LeftJustify)
{
   char format[16];
   
   if(p==NULL)
   {
      label[0] = '\0';
   }
   else
   {
      if((strlen(p->chain) > 1) || isdigit(p->chain[0]))
      {
         sprintf(label,"%s.%d%s", p->chain, p->resnum, p->insert);
      }
      else
      {
         sprintf(label,"%s%d%s", p->chain, p->resnum, p->insert);
      }
   }
   
   /* If left justifying then the format string specifies the width     */
   if(LeftJustify)
      sprintf(format, "%%-%ds", width);
   else
      strcpy(format, "%s");

   sprintf(label, format, label);
}


/************************************************************************/
/*>REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians)
   ------------------------------------------------------------------
*//**
   \param[in]    *p1      Pointer to PDB record
   \param[in]    *p2      Pointer to PDB record
   \param[in]    *p3      Pointer to PDB record
   \param[in]    *p4      Pointer to PDB record
   \param[in]    Radians  Output radians rather than degrees?
   \return                Torsion angle

   Wraps around the blPhi() function to take (and check) PDB pointers
   rather than coordinates, handle conversion from radians, etc.

- 27.11.14 Original   By: ACRM
*/
REAL CalcTorsion(PDB *p1, PDB *p2, PDB *p3, PDB *p4, BOOL Radians)
{
   REAL tor;
   
   if((p1==NULL)||(p2==NULL)||(p3==NULL)||(p4==NULL))
   {
      return(ERROR_VALUE);
   }

   tor = blPhi(p1->x, p1->y, p1->z,
               p2->x, p2->y, p2->z,
               p3->x, p3->y, p3->z,
               p4->x, p4->y, p4->z);

   if(!Radians)
      tor *= 180/PI;

   return(tor);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     BOOL *CATorsions, BOOL *terse, BOOL *Radians, 
                     BOOL *oldStyle)
   ---------------------------------------------------------------------
*//**

   \param[in]     argc         Argument count
   \param[in]     **argv       Argument array
   \param[out]    *infile      Input file (or blank string)
   \param[out]    *outfile     Output file (or blank string)
   \param[out]    *CATorsions  Do CA pseudo-torsions
   \param[out]    *terse       Terse (1-letter AA code) output
   \param[out]    *Radians     Output radians rather than degrees
   \param[out]    *oldStyle    Old style output
   \return                     Success?

   Parse the command line
   
-  05.02.96 Original    By: ACRM
-  27.02.14 V2.0
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *CATorsions, BOOL *terse, BOOL *Radians, 
                  BOOL *oldStyle)
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
         case 'c':
            *CATorsions = TRUE;
            break;
         case 't':
            *terse = TRUE;
            break;
         case 'r':
            *Radians = TRUE;
            break;
         case 'o':
            *oldStyle = TRUE;
            break;
         case 'n':
            *oldStyle = FALSE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are <= 2 arguments left                    */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
         }
         
         /* Copy the second to outfile                                  */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
         }
         
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

   Displays a usage message

-  10.06.94 original   By: ACRM
-  16.08.94 Added -m
-  12.06.95 Added -c
-  22.07.14 V1.3 By: CTP
-  06.11.14 V1.5 By: ACRM
-  07.11.14 V1.6
-  27.11.14 V2.0
-  04.03.15 V2.1
*/
void Usage(void)
{
   fprintf(stderr,"\npdbtorsions V2.1 (c) 1994-2015 Andrew Martin, \
UCL.\n");
   fprintf(stderr,"\nUsage: pdbtorsions [-h][-r][-c][-t][-o][-n] \
[in.pdb [out.tor]]\n");
   fprintf(stderr,"       -h   This help message\n");
   fprintf(stderr,"       -r   Give results in radians\n");
   fprintf(stderr,"       -c   Generate CA-CA pseudo-torsions\n");
   fprintf(stderr,"       -t   Terse format - use 1-letter code\n");
   fprintf(stderr,"       -o   Old format (see below)\n");
   fprintf(stderr,"       -n   New format (see below)\n");

   fprintf(stderr,"\nGenerates a set of backbone torsions from a PDB \
file.\n\n");
   fprintf(stderr,"I/O is through stdin/stdout if unspecified.\n");

   fprintf(stderr,"\nV1.x of this program associated the omega torsion angle with the residue\n");
   fprintf(stderr,"before the torsion instead of the standard way of associating it with\n");
   fprintf(stderr,"the residue after. In addition chain labels were not displayed since\n");
   fprintf(stderr,"the code did not handle multiple chains correctly (i.e. it displayed\n");
   fprintf(stderr,"non-existent torsion angles between the residues at the termini of chains\n");
   fprintf(stderr,"since it assumed everything was a single chain).\n");

   fprintf(stderr,"\nV2.x corrects the association of the omega torsion angle and changes the \n");
   fprintf(stderr,"output format to include the chain label in the residue number. It treats\n");
   fprintf(stderr,"multiple chains correctly. The old behaviour of associating the omega\n");
   fprintf(stderr,"angle with the preceding residue and the old output format can be \n");
   fprintf(stderr,"obtained by using the -o (old) flag. However the chain breaks are still\n");
   fprintf(stderr,"handled correctly.\n");

   fprintf(stderr,"\nThe old behaviour is also obtained if the executable is named 'torsions'\n");
   fprintf(stderr,"rather than 'pdbtorsions'. In that case the new behaviour can be obtained\n");
   fprintf(stderr,"by using the -n (new) flag.\n\n");
}


