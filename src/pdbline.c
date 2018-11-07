/************************************************************************/
/**

   \file       pdbline.c
   
   \version    V1.0
   \date       08.10.14
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

*************************************************************************/
/* Includes
*/
/* Includes
*/
#include <stdio.h> 
#include <stdlib.h>
#include <unistd.h>

#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/array.h"
#include "regression.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 540

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
REAL **BuildCaCoordArray(PDB *pdb, int *numCa);
PDB *ExtractZoneSpecPDB(PDB *pdb, char *firstRes, char *lastRes);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *firstRes, char *lastRes);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program

-  08.10.14  Original   By: ACRM
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
   char    infile[MAXBUFF], outfile[MAXBUFF],
           firstRes[MAXBUFF], lastRes[MAXBUFF];

  if(ParseCmdLine(argc, argv, infile, outfile, firstRes, lastRes))
  {
     if(OpenStdFiles(infile, outfile, &in, &out))
     {
        /* Read PDB file                                                */
        if((pdb = ReadPDB(in, &natoms)) == NULL)
        {
           fprintf(stderr, "No atoms read from PDB file.\n");
           retval = 1;
        }
        else
        {
           /* Extract the zone of interest                              */
           if((zone = ExtractZoneSpecPDB(pdb, firstRes, lastRes))==NULL)
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
              ComputeBestFitLine(coordArray, numCa, 3, 
                                 Centroid, EigenVector);
              DrawRegressionLine(out, coordArray, EigenVector, 
                                 numCa, "X");
              WritePDB(out,zone);
              FreeArray2D((char **)coordArray, numCa, 3);
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
   \param[in]   PDB  *pdb    PDB linked list for region of interest
   \param[out]  int  *numCa  Number of CA atoms in the region
   \return      REAL **      2D array of coordinates for the CA atoms
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
   if((capdb = SelectAtomsPDB(pdb, 1, sel, numCa))!=NULL)
   {
      if((coordArray = (REAL **)Array2D(sizeof(REAL), *numCa, 3))==NULL)
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
/*>PDB *ExtractZoneSpecPDB(PDB *pdb, char *firstRes, char *lastRes)
   ----------------------------------------------------------------
*//**
   \param[in]   PDB  *pdb       PDB linked list
   \param[in]   char *firstRes  Residue spec ([chain]resnum[insert])
   \param[in]   char *lastRes   Residue spec ([chain]resnum[insert])

   Extracts a zone from a PDB linked list, making a copy of the original
   list.

-  08.10.14  Original   By: ACRM
*/
PDB *ExtractZoneSpecPDB(PDB *pdb, char *firstRes, char *lastRes)
{
   char chain1[8],  chain2[8],
        insert1[8], insert2[8];
   int  resnum1,    resnum2;
   PDB  *zone = NULL;

   if(ParseResSpec(firstRes, chain1, &resnum1, insert1) &&
      ParseResSpec(lastRes,  chain2, &resnum2, insert2))
   {
      zone = ExtractZonePDB(pdb, 
                            chain1, resnum1, insert1,
                            chain2, resnum2, insert2);
   }
   return(zone);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     char *firstRes, char *lastRes)
   ----------------------------------------------------------------------
*//**
   \param[in]  int    argc        Argument count
   \param[in]  char   **argv      Argument array
   \param[out] char   *infile     Input filename (or blank string)
   \param[out] char   *outfile    Output filename (or blank string)
   \param[out] char   *firstRes   First residue spec
   \param[out] char   *lastRes    Last residue spec
   \return     BOOL               Success

   Parse the command line

-  08.10.14 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *firstRes, char *lastRes)
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
/*>void Usage(void)
   ----------------
*//**
   Prints a usage message

-  08.10.14  Original   By: ACRM
*/
void Usage(void)
{
   printf("\npdbline V1.0 (c) 2014 UCL, Dr. Andrew C.R. Martin\n");
   printf("        With contributions from Abhi Raghavan and Saba \
Ferdous\n");

   printf("\nUsage: pdbline firstres lastres [in.pdb [out.pdb]]\n");
   printf("       firstres - a residue identifier of the form \
[chain]resnum[insert]\n");
   printf("                  representing the first residue of \
interest\n");
   printf("       lastres  - a residue identifier of the form \
[chain]resnum[insert]\n");
   printf("                  representing the last residue of \
interest\n");

   printf("\nGenerates a set of atom positions along a best fit line \
through a\n");
   printf("specified set of C-alpha atoms. Input and output are through \
standard\n");
   printf("input/output if files are not specified\n\n");
}
