/************************************************************************/
/**

   \file       pdbline.c
   
   \version    V1.2
   \date       24.10.14
   \brief      Draws a best fit line through a specified set of CA atoms
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 2014
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
-  V1.0   08.10.14 Original   By: ACRM Based on code by Abhi Raghavan
                   and Saba Ferdous
-  V1.1   13.10.14 Added -r and -a options to specify the residue label
                   and atom label for the line residues
-  V1.2   24.10.14 Fixed bug in plotting line when Eigen vector had X 
                   dimension close to zero. Added -v option
                   Moved regression code into BiopLib

*************************************************************************/
/* Includes
*/
#include <stdio.h> 
#include <stdlib.h>
#include <unistd.h>

#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/array.h"
#include "bioplib/macros.h"
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/regression.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 540
#define DEF_RESNAM "LIN"
#define DEF_ATNAM  "X"

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
REAL **BuildCaCoordArray(PDB *pdb, int *numCa);
BOOL DrawPDBRegressionLine(FILE *wfp, double **coordinates,
                           double *eigenVector, int numberOfPoints,
                           char *chainLabel, char *resnam, char *atnam);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *firstRes, char *lastRes, 
                  char *resnam, char *atnam, BOOL *verbose);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program

-  08.10.14  Original   By: ACRM
-  13.10.14  Added resnam and atnam
*/
int main(int argc, char **argv)
{
   int     natoms       = 0,
           retval       = 0;
   PDB     *pdb         = NULL,
           *zone        = NULL;
   REAL    **coordArray = NULL;
   FILE    *in          = stdin,
           *out         = stdout;
   BOOL    verbose      = FALSE;
   char    infile[MAXBUFF], 
           outfile[MAXBUFF],
           firstRes[MAXBUFF], 
           lastRes[MAXBUFF],
           resnam[MAXBUFF],
           atnam[MAXBUFF];
   

   strcpy(resnam, DEF_RESNAM);
   PADMINTERM(resnam, 4);
   strcpy(atnam,  DEF_ATNAM);
   PADMINTERM(resnam, 4);

   if(ParseCmdLine(argc, argv, infile, outfile, firstRes, lastRes, resnam,
                   atnam, &verbose))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         /* Read PDB file                                               */
         if((pdb = blReadPDB(in, &natoms)) == NULL)
         {
            fprintf(stderr, "No atoms read from PDB file.\n");
            retval = 1;
         }
         else
         {
            /* Extract the zone of interest                             */
            if((zone = blExtractZoneSpecPDBAsCopy(pdb, firstRes,
                                                  lastRes))==NULL)
            {
               fprintf(stderr, "Unable to extract specified zone from \
PDB file\n");
               retval = 1;
            }
            else
            {
               int  numCa;
               REAL EigenVector[3],
                    Centroid[3];
               
               coordArray = BuildCaCoordArray(zone, &numCa);
               blCalculateBestFitLine(coordArray, numCa, 3, 
                                      Centroid, EigenVector);
               if(verbose)
               {
                  fprintf(stderr,"Centroid:     %8.3f %8.3f %8.3f\n",
                          Centroid[0], Centroid[1], Centroid[2]);
                  fprintf(stderr,"Eigen vector: %8.3f %8.3f %8.3f\n",
                          EigenVector[0], EigenVector[1], EigenVector[2]);
               }
               
               DrawPDBRegressionLine(out, coordArray, EigenVector, 
                                     numCa, "X", resnam, atnam);
               blWritePDB(out,zone);
               blFreeArray2D((char **)coordArray, numCa, 3);
               FREELIST(zone, PDB);
            }
            FREELIST(pdb, PDB);
         }
      }
   }
   else
   {
      Usage();
   }
   
   return(retval);
}


/************************************************************************/
/*>REAL **BuildCaCoordArray(PDB *pdb, int *numCa)
   ----------------------------------------------
*//**
   \param[in]   pdb    PDB linked list for region of interest
   \param[out]  numCa  Number of CA atoms in the region
   \return             2D array of coordinates for the CA atoms
                       in the input PDB linked list

   Extracts the CA atoms from the PDB linked list and obtains their
   coordinates. These are stored in an allocated 2D array which is
   returned by the routine

-  08.10.14  Original   By: ACRM
*/
REAL **BuildCaCoordArray(PDB *pdb, int *numCa)
{
   char *sel[1];
   PDB  *capdb = NULL;
   int  count = 0;
   REAL **coordArray = NULL;
   PDB  *p;

   SELECT(sel[0],"CA  ");
   if((capdb = blSelectAtomsPDBAsCopy(pdb, 1, sel, numCa))!=NULL)
   {
      if((coordArray = (REAL **)blArray2D(sizeof(REAL), *numCa, 3))==NULL)
      {
         *numCa = 0;
         return(NULL);
      }

      for(p=capdb, count=0; p!=NULL; NEXT(p))
      {
         coordArray[count][0] = p->x;
         coordArray[count][1] = p->y;
         coordArray[count][2] = p->z;
         count++;
      }
      FREELIST(capdb, PDB);
   }

   return(coordArray);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     char *firstRes, char *lastRes,
                     char *resnam, char *atnam, BOOL *verbose)
   ----------------------------------------------------------------------
*//**
   \param[in]     argc       Argument count
   \param[in]     argv       Argument array
   \param[out]    infile     Input filename (or blank string)
   \param[out]    outfile    Output filename (or blank string)
   \param[out]    firstRes   First residue spec
   \param[out]    lastRes    Last residue spec
   \param[out]    resnam     Residue name for the line
   \param[out]    atnam      Atom name for the line
   \param[out]    verbose    Verbose?
   \return                   Success?

   Parse the command line

-  08.10.14  Original    By: ACRM
-  13.10.14  Added -r resnam and -a atnam
-  24.10.14  Added -v verbose
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *firstRes, char *lastRes,
                  char *resnam, char *atnam, BOOL *verbose)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   *verbose  = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'r':
            argc--; argv++;
            if(!argc) return(FALSE);
            strncpy(resnam, argv[0], MAXBUFF);
            UPPER(resnam);
            break;
         case 'a':
            argc--; argv++;
            if(!argc) return(FALSE);
            strncpy(atnam, argv[0], MAXBUFF);
            UPPER(atnam);
            break;
         case 'v':
            *verbose = TRUE;
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
         /* Check that there are 2-4 arguments left                     */
         if((argc < 2) || (argc > 4))
            return(FALSE);
         
         /* Copy the first two to firstRes and lastRes                  */
         strcpy(firstRes, argv[0]);
         argc--; argv++;
         strcpy(lastRes,  argv[0]);
         argc--; argv++;
         
         /* If there's another, copy it to infile                       */
         if(argc)
            strcpy(infile, argv[0]);
         argc--; argv++;

         /* If there's another, copy it to outfile                      */
         if(argc)
            strcpy(outfile, argv[0]);

         return(TRUE);
      }
      argc--; argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>BOOL DrawPDBRegressionLine(FILE *wfp, REAL **coordinates, 
                              REAL *eigenVector,
                              int numberOfPoints, char *chainLabel,
                              char *resnam, char *atnam)
   ----------------------------------------------------------------
*//**
   \param[in]  wfp            File pointer to write into a PDB file.
   \param[in]  coordinates    Coordinates of points used to compute 
                              the regression line
   \param[in]  eigenVector    Eigen vector components of the 
                              regression line
   \param[in]  numberOfPoints Number of points (number of rows in 
                              coordinates array).
   \param[in]  chainLabel     Chain label when writing to PDB file.
   \param[in]  resnam         Residue name for the line
   \param[in]  atnam          Atom name for the line
   \return                    Success?

   Outputs a line in the form of PDB records. 

   From the coordinates supplied, the centroid is calculated. The eigen 
   vector components of the regression line are used to find other 
   points on the line that passes through the centroid.

-  06.10.14  Original   By: ACRM
-  13.10.14  Added resnam and atnam
-  24.10.14  Modified to walk along the direction of the largest Eigen
             vector dimension rather than always walking along X. The
             previous code broke if there X dimension was close to zero
*/
BOOL DrawPDBRegressionLine(FILE *wfp, REAL **coordinates, 
                           REAL *eigenVector,
                           int numberOfPoints, char *chainLabel,
                           char *resnam, char *atnam)
{
   PDB  *p            = NULL,
        *pdb          = NULL;
   int  kmin          = 0,
        kmax          = 0,
        i             = 0,
        largestDim    = 0;
   REAL smallestCoord = 0.0,
        largestCoord  = 0.0,
        centroid[3];

   /* Find the centroid of the atom positions                           */
   blFindCentroid(coordinates, numberOfPoints, 3, centroid);

   /* Determine which dimension of the eigen vector is largest          */
   largestDim = 0;
   for(i=1; i<3; i++)
   {
      if(ABS(eigenVector[i]) > ABS(eigenVector[largestDim]))
         largestDim = i;
   }

   /* Find the smallest and largest coordinate in the largest dimension */
   smallestCoord = largestCoord = coordinates[0][largestDim];
   for(i=1; i<numberOfPoints; i++)
   {
      if(coordinates[i][largestDim] > largestCoord)
         largestCoord = coordinates[i][largestDim];
      if(coordinates[i][largestDim] < smallestCoord)
         smallestCoord = coordinates[i][largestDim];
   }
      
   /* Find values of k such that points on the line of best fit may be 
      plotted 
   */
   kmin = (int)((smallestCoord - centroid[largestDim]) /
                eigenVector[largestDim]);
   kmax = (int)((largestCoord  - centroid[largestDim]) /
                eigenVector[largestDim]);
   if(kmin > kmax)
   {
      kmin = -1 * kmin;
      kmax = -1 * kmax;
   }

   /* Create a PDB linked list for the points along the line            */
   for(i=kmin; i<=kmax; i++)
   {
      if(pdb == NULL)
      {
         INIT(pdb, PDB);
         p=pdb;
      }
      else
      {
         ALLOCNEXT(p, PDB);
      }
      if(p==NULL)
      {
         FREELIST(pdb, PDB);
         return(FALSE);
      }
      CLEAR_PDB(p);
      
      strcpy(p->chain,chainLabel);
      p->x      = centroid[0]+((REAL)i * eigenVector[0]);
      p->y      = centroid[1]+((REAL)i * eigenVector[1]);
      p->z      = centroid[2]+((REAL)i * eigenVector[2]);
      p->occ    = 1.0;
      p->bval   = 1.0;
      p->next   = NULL;
      p->atnum  = i-kmin+1;
      p->resnum = i-kmin+1;

      strcpy(p->record_type,"ATOM");

      /* Atom name                                                      */
      strncpy(p->atnam, atnam, 8);
      p->atnam[4] = '\0';

      /* Raw atom name                                                  */
      p->atnam_raw[0] = ' ';
      strncpy(p->atnam_raw+1,atnam,3);
      p->atnam_raw[4] = '\0';

      /* Residue name                                                   */
      strncpy(p->resnam,resnam, 8);
      p->resnam[4] = '\0';

      strcpy(p->insert," ");
      p->altpos=' ';
   }

   /* Write the points along the line out in PDB format                 */
   for(p=pdb; p!=NULL; NEXT(p))
      blWritePDBRecord(wfp,p);

   /* Free memory for linked list and return                            */
   FREELIST(pdb, PDB);

   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**
   Prints a usage message

-  08.10.14  Original   By: ACRM
-  13.10.14  Added -r and -a options
-  24.10.14  Added -v option
*/
void Usage(void)
{
   printf("\npdbline V1.2 (c) 2014 UCL, Dr. Andrew C.R. Martin\n");
   printf("        With contributions from Abhi Raghavan and Saba \
Ferdous\n");

   printf("\nUsage: pdbline [-r resnam][-a atnam][-v] firstres lastres \
[in.pdb [out.pdb]]\n");
   printf("       firstres - a residue identifier of the form \
[chain]resnum[insert]\n");
   printf("                  representing the first residue of \
interest\n");
   printf("       lastres  - a residue identifier of the form \
[chain]resnum[insert]\n");
   printf("                  representing the last residue of \
interest\n");
   printf("       -r Specify the residue name for the line \
(Default: %s)\n", DEF_RESNAM);
   printf("       -a Specify the atom name for the line \
(Default: %s)\n", DEF_ATNAM);
   printf("       -v Verbose: Prints the Eigen vector and centroid\n");

   printf("\nGenerates a set of atom positions along a best fit line \
through a\n");
   printf("specified set of C-alpha atoms. Input and output are through \
standard\n");
   printf("input/output if files are not specified\n\n");
}


