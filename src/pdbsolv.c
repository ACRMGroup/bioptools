/************************************************************************/
/**

   \file       pdbsolv.c
   
   \version    V1.6
   \date       21.11.17
   \brief      Solvent accessibility using bioplib
   
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-   V1.0   16.07.14 Original   By: ACRM
-   V1.1   19.08.14 Fixed call to renamed function: 
                    blStripWatersPDBAsCopy() By: CTP
-   V1.2   06.11.14 Renamed from solv  By: ACRM
-   V1.3   13.02.15 Modified to use whole PDB
-   V1.4   17.06.15 Residue access file now includes sidechain 
                    accessibility
-   V1.5   08.03.16 Corrected insert code printing so it is left-justified
                    and now touches the residue number
-   V1.7   21.11.17 Added -x flag to add radii in occupancy column

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/access.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define DEF_PROBERADIUS 1.4
#define DEF_RADFILE "radii.dat"
#define DATA_ENV "DATADIR"

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *integrationAccuracy, REAL *rad, char *radfile,
                  BOOL *doAccessibility, char *resfile, BOOL *noAtoms,
                  BOOL *addRadii);
void Usage(void);
void PopulateBValWithAccess(PDB *pdb);
void PopulateOccWithRadii(PDB *pdb);
void PrintResidueAccessibility(FILE *out, PDB *pdb, RESRAD *resrad);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//** 
   Main program for solvent accessibility calculations

-  17.07.14 Original   By: ACRM
-  19.08.14 Fixed call to renamed function: blStripWatersPDBAsCopy()
                  By: CTP
-  13.02.15 Modified to use whole PDB   By: ACRM

*/
int main(int argc, char **argv)
{
   RESRAD   *resrad;
   FILE     *in     = stdin,
            *out    = stdout,
            *resout = stdout,
            *fpRad  = NULL;
   int      natoms;
   WHOLEPDB *wpdb;
   PDB      *pdb;
   BOOL     doAccessibility = FALSE,
            noenv           = FALSE,
            noAtoms         = FALSE,
            doResaccess     = FALSE,
            addRadii        = FALSE;
   REAL     integrationAccuracy,
            probeRadius;
   char     infile[MAXBUFF],
            outfile[MAXBUFF],
            radfile[MAXBUFF],
            resfile[MAXBUFF];
   
   if(!ParseCmdLine(argc, argv, infile, outfile, 
                    &integrationAccuracy, &probeRadius, 
                    radfile, &doAccessibility, resfile, &noAtoms,
                    &addRadii))
   {
      Usage();
      return(0);
   }

   if(resfile[0] != '\0')
   {
      doResaccess = TRUE;
      if((resout = blOpenOrPipe(resfile))==NULL)
      {
         fprintf(stderr, "Error (pdbsolv): Unable to open file or pipe \
for residue accessibility data (%s)\n", resfile);
         return(1);
      }
   }

   if(!blOpenStdFiles(infile, outfile, &in, &out))
   {
      fprintf(stderr, "Error (pdbsolv): Unable to open input or output \
file\n");
      return(1);
   }

   if((wpdb = blReadWholePDB(in))==NULL)
   {
      fprintf(stderr, "Error (pdbsolv): No atoms read from PDB \
file, %s\n", infile);
      return(1);
   }

   /* Strip waters                                                      */
   if((pdb = blStripWatersPDBAsCopy(wpdb->pdb, &natoms))==NULL)
   {
      fprintf(stderr, "Error (pdbsolv): No memory to strip waters from \
PDB file, %s\n",
              infile);
      return(1);
   }

   /* Free the original linked list of atoms and patch in the new one   */
   FREELIST(wpdb->pdb, PDB);
   wpdb->pdb = pdb;

   /* Open the radius file                                              */
   if((fpRad=blOpenFile(radfile, DATA_ENV, "r", &noenv))==NULL)
   {
      fprintf(stderr, "Error (pdbsolv): Unable to open radius file, \
%s\n", radfile);
      if(noenv)
      {
         fprintf(stderr, "              Environment variable %s \
not set\n", DATA_ENV);
      }
      return(1);
   }

   /* Set the atom radii in the linked list                             */
   resrad = blSetAtomRadii(pdb, fpRad);

   /* Do the actual accessibility calculations                          */
   if(!blCalcAccess(pdb, natoms, 
                    integrationAccuracy, probeRadius,
                    doAccessibility))
   {
      fprintf(stderr,"Error: (pdbsolv) No memory for accessibility \
arrays\n");
      return(1);
   }

   /* And populate the B-values with the accessibility and write the
      new PDB file
   */
   if(!noAtoms)
   {
      PopulateBValWithAccess(pdb);
      if(addRadii)
      {
         PopulateOccWithRadii(pdb);
      }
      
      blWriteWholePDB(out, wpdb);
   }

   if(doResaccess)
   {
      PrintResidueAccessibility(resout, pdb, resrad);
      blCloseOrPipe(resout);
   }

   /* Free up the memory for the PDB linked list                        */
   FREELIST(pdb, PDB);
   /* Free up the memory from the residue radii                         */
   FREELIST(resrad, RESRAD);

   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     REAL *p, REAL *rad, char *radfile,
                     BOOL *doAccessibility, char *resfile, BOOL *noAtoms,
                     BOOL *addRadii)
   ----------------------------------------------------------------------
*//**
   \param[in]   int    argc              Argument count
   \param[in]   char   **argv            Argument array
   \param[out]  char   *infile           Input filename (or blank string)
   \param[out]  char   *outfile          Output filename (or blank string)
   \param[out]  REAL   *p                Integration accuracy parameter
   \param[out]  REAL   *rad              Radius of probe
   \param[out]  char   *radfile          File of atom radii
   \param[out]  BOOL   *doAccessibility  Calculate accessibility rather than 
                                         contact area
   \param[out]  char   *resfile          File for storing residue 
                                         accessibilities
   \param[out]  BOOL   *noAtoms          Do not write atom accessibilities
   \param[out]  BOOL   *addRadii         Add radii to occupancy column
   \return      BOOL                     Success

   Parse the command line

   17.07.14 Original    By: ACRM
   21.11.17 Added -x addRadii
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *integrationAccuracy, REAL *rad, char *radfile,
                  BOOL *doAccessibility, char *resfile, BOOL *noAtoms,
                  BOOL *addRadii)
{
   argc--;
   argv++;
   
   *doAccessibility     = TRUE;
   *integrationAccuracy = ACCESS_DEF_INTACC;
   *rad                 = DEF_PROBERADIUS;
   *noAtoms             = FALSE;

   infile[0] = outfile[0] = radfile[0] = resfile[0] = '\0';
   strcpy(radfile, DEF_RADFILE);
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'i':
            if(!(--argc) || !sscanf((++argv)[0],"%lf",
                                    integrationAccuracy))
               return(FALSE);
            break;
         case 'p':
            if(!(--argc) || !sscanf((++argv)[0],"%lf",rad))
               return(FALSE);
            break;
         case 'f':
            if(!(--argc))
               return(FALSE);
            strncpy(radfile,(++argv)[0],MAXBUFF);
            break;
         case 'r':
            if(!(--argc))
               return(FALSE);
            strncpy(resfile,(++argv)[0],MAXBUFF);
            break;
         case 'n':
            *noAtoms = TRUE;
            break;
         case 'c':
            *doAccessibility = FALSE;
            break;
         case 'x':
            *addRadii = TRUE;
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
   Prints a usage message

-   17.07.14 Original   By: ACRM
-   06.11.14 V1.2
-   13.02.15 V1.3
-   17.06.15 V1.4
-   08.03.16 V1.5
-   21.11.17 V1.6
*/
void Usage(void)
{
   fprintf(stderr,"\npdbsolv V1.6 (c) 2014-2017 UCL, Dr. Andrew C.R. \
Martin\n");

   fprintf(stderr,"\nUsage: pdbsolv [-i val] [-p val] [-f radfile] \
[-r resfile] [-n] [-c] [-x] [in.pdb [out.pdb]]\n");
   fprintf(stderr,"            -i val      Specify integration accuracy \
(Default: %.2f)\n",ACCESS_DEF_INTACC);
   fprintf(stderr,"            -p val      Specify probe radius \
(Default: %.2f)\n",DEF_PROBERADIUS);
   fprintf(stderr,"            -f radfile  Specify radius file\n");
   fprintf(stderr,"                        (Default: %s)\n",
           DEF_RADFILE);
   fprintf(stderr,"            -r resfile  Specify file for saving \
residue accessibility \n");
   fprintf(stderr,"                        data. If the file is \
specified as 'stdout' then\n");
   fprintf(stderr,"                        data will be written to \
standard output. If the\n");
   fprintf(stderr,"                        filename starts with a pipe \
symbol (|), data will\n");
   fprintf(stderr,"                        be piped to the specified \
program.\n");
   fprintf(stderr,"            -n          Do not print atom \
accessibility. Used with -r\n");
   fprintf(stderr,"            -c          Do contact area instead of \
accessibility\n");
   fprintf(stderr,"            -x          Add radii in occupancy \
column of PDB file\n");


   fprintf(stderr,"\nPerforms solvent accessibility calculations \
according to the method of\n");
   fprintf(stderr,"Lee and Richards. Reads and writes PDB format files. \
Input/output is\n");
   fprintf(stderr,"to standard input/output if files are not \
specified.\n\n");
}


/************************************************************************/
/*>void PopulateBValWithAccess(PDB *pdb)
   -------------------------------------
*//**
   \param   PDB  *pdb    PDB linked list

   Copies the accessibility inforation into the B-Value column for output

-  17.07.14  Original   By: ACRM   
*/
void PopulateBValWithAccess(PDB *pdb)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->bval = p->access;
   }
}


/************************************************************************/
/*>void PopulateOccWithRadii(PDB *pdb)
   -----------------------------------
*//**
   \param   PDB  *pdb    PDB linked list

   Copies the radius information into the occupancy column for output

-  21.11.17  Original   By: ACRM   
*/
void PopulateOccWithRadii(PDB *pdb)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->occ = p->radius;
   }
}


/************************************************************************/
/*>void PrintResidueAccessibility(FILE *out, PDB *pdb, RESRAD *resrad)
   -------------------------------------------------------------------
*//**
   \param[in]  FILE   *out     Output file pointer
   \param[in]  PDB    *pdb     PDB linked list
   \param[in]  RESRAD *resrad  Radius information and standard 
                               accessibilities

   Calls blCalcResAccess() to calculate residue accessibilities and 
   prints the results

-  17.07.14  Original   By:ACRM
-  17.06.15  Added sidechain accessibility printing
-  08.03.16  Corrected insertion printing so it is left justified
*/
void PrintResidueAccessibility(FILE *out, PDB *pdb, RESRAD *resrad)
{
   RESACCESS *resaccess;
   if((resaccess = blCalcResAccess(pdb, resrad))==NULL)
   {
      fprintf(stderr, "Error: (pdbsolv) Unable to allocate memory for \
residue accessibilities\n");
   }
   else
   {
      RESACCESS *r;
      fprintf(out, "#       RESIDUE  AA   ACCESS  RELACC  SCACC   SCRELACC\n");

      for(r=resaccess; r!=NULL; NEXT(r))
      {
         fprintf(out, "RESACC %2s%5d%-2s %s %7.3f %7.3f %7.3f %7.3f\n",
                 r->chain, r->resnum, r->insert, r->resnam, 
                 r->resAccess, r->relAccess,
                 r->scAccess,  r->scRelAccess);
      }
   }
}

