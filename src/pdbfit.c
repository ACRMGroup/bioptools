/*************************************************************************

   Program:    pdbfit
   File:       pdbfit.c
   
   Version:    V2.0
   Date:       03.11.17
   Function:   Simple program to fit 2 sets of coordinates for an 
               identical protein
   
   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 2001-2017
   Author:     Dr. Andrew C. R. Martin
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
   V1.0   12.12.01  Original pdbcafit (12.05.10) and pdbfit (12.12.01)
   V2.0   03.11.17  Combined pdbfit and pdbcafit
 
*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF  512
#define TYPE_ALL 0
#define TYPE_CA  1
#define TYPE_BB  2

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile1, char *infile2, 
                  int *type, BOOL *showCoords);
REAL CalcRMSPDBOverType(PDB *pdbin1, PDB *pdbin2, int type);

/************************************************************************/
int main(int argc, char **argv)
{
   char infile1[MAXBUFF],
        infile2[MAXBUFF];
   int  type       = TYPE_ALL;
   BOOL showCoords = FALSE;
   
   if(ParseCmdLine(argc, argv, infile1, infile2, &type, &showCoords))
   {
      FILE *in1, *in2;
      int  natoms1, natoms2;
      PDB  *pdb1,   *pdb2;
      REAL rm[3][3], rms;
      
      
      if((in1 = fopen(infile1, "r")) == NULL)
      {
         fprintf(stderr, "Error: Unable to read file (%s)\n", infile1);
         exit(1);
      }
      if((in2 = fopen(infile2, "r")) == NULL)
      {
         fprintf(stderr, "Error: Unable to read file (%s)\n", infile2);
         exit(1);
      }
      
            
      if((pdb1=blReadPDB(in1, &natoms1))==NULL)
      {
         fprintf(stderr,"Error: Can't read atoms from %s\n",argv[1]);
         return(1);
      }
      if((pdb2=blReadPDB(in2, &natoms2))==NULL)
      {
         fprintf(stderr,"Error: Can't read atoms from %s\n",argv[2]);
         return(1);
      }

      if(natoms1 != natoms2)
      {
         fprintf(stderr,"Non-identical PDB lists\n");
         return(1);
      }
   
      switch(type)
      {
      case TYPE_CA:
         if(!blFitCaPDB(pdb1, pdb2, rm))
         {
            fprintf(stderr,"Error: Unable to fit structures\n");
            return(1);
         }
         break;
      case TYPE_BB:
         if(!blFitNCaCPDB(pdb1, pdb2, rm))
         {
            fprintf(stderr,"Error: Unable to fit structures\n");
            return(1);
         }
         break;
      default:
         if(!blFitPDB(pdb1, pdb2, rm))
         {
            fprintf(stderr,"Error: Unable to fit structures\n");
            return(1);
         }
         break;
      }

      rms = CalcRMSPDBOverType(pdb1, pdb2, type);
      
      printf("RMSD  %.3f\n", rms);
      
      if(showCoords)
         blWritePDB(stdout, pdb2);
      
      FREELIST(pdb1, PDB);
      FREELIST(pdb2, PDB);
   }
   else
   {
      Usage();
   }
   
   return(0);
}


/************************************************************************/
/*>REAL CalcRMSPDBOverType(PDB *pdbin1, PDB *pdbin2, int type)
   -----------------------------------------------------------
*//**
   \param[in]   *pdbin1   First PDB linked list
   \param[in]   *pdbin2   Second PDB linked list
   \param[in]   type      Type of fit (atoms over which to calculate
   \return                RMSD

   Calculates the RMSD over the approproate set of atoms

-  03.11.17 Original   By: ACRM
*/
REAL CalcRMSPDBOverType(PDB *pdbin1, PDB *pdbin2, int type)
{
   char *sel[8];
   int  natoms1,
        natoms2;
   PDB  *pdb1 = NULL,
        *pdb2 = NULL;
   REAL rms;
   

   switch(type)
   {
   case TYPE_CA:
      SELECT(sel[0], "CA  ");
      pdb1 = blSelectAtomsPDBAsCopy(pdbin1, 1, sel, &natoms1);
      pdb2 = blSelectAtomsPDBAsCopy(pdbin2, 1, sel, &natoms2);
      rms = blCalcRMSPDB(pdb1, pdb2);
      free(sel[0]);
      FREELIST(pdb1, PDB);
      FREELIST(pdb2, PDB);
      break;
   case TYPE_BB:
      SELECT(sel[0], "N   ");
      SELECT(sel[1], "CA  ");
      SELECT(sel[2], "C   ");
      pdb1 = blSelectAtomsPDBAsCopy(pdbin1, 3, sel, &natoms1);
      pdb2 = blSelectAtomsPDBAsCopy(pdbin2, 3, sel, &natoms2);
      free(sel[0]);
      free(sel[1]);
      free(sel[2]);
      rms = blCalcRMSPDB(pdb1, pdb2);
      FREELIST(pdb1, PDB);
      FREELIST(pdb2, PDB);
      break;
   default:
      rms = blCalcRMSPDB(pdbin1, pdbin2);
      break;
   }

   return(rms);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     int *type, BOOL *showCoords)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *type        Type of fit
   \param[out]     *showCoords  Show the fitted coordinates?

   \return                      Success?

   Parse the command line
   
-  03.11.17 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile1, char *infile2, 
                  int *type, BOOL *showCoords)
{
   argc--;
   argv++;

   *type = TYPE_ALL;
   infile1[0] = infile2[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'c':
            if(*type != TYPE_ALL)
            {
               fprintf(stderr, "Warning: -c ignored as -b already specified\n");
            }
            else
            {
               *type = TYPE_CA;
            }
            break;
         case 'b':
            if(*type != TYPE_ALL)
            {
               fprintf(stderr, "Warning: -b ignored as -c already specified\n");
            }
            else
            {
               *type = TYPE_BB;
            }
            break;
         case 'w':
            *showCoords = 1;
            break;
         case 'h':
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2 arguments left                    */
         if(argc != 2)
            return(FALSE);
         
         /* Copy the next to infile                                  */
         strcpy(infile1, argv[0]);
         argc--;
         argv++;

         strcpy(infile2, argv[0]);
         argc--;
         argv++;
         
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

-  03.11.17  Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\npdbfit V2.0 (c) 2001-2017, UCL, Dr. Andrew C.R. \
Martin\n");
   fprintf(stderr,"\nUsage: pdbcafit [-c|-b][-w] file1.pdb file2.pdb\n");
   fprintf(stderr,"       -c Fit only C-alphas\n");
   fprintf(stderr,"       -b Fit only backbone (N,CA,C)\n");
   fprintf(stderr,"       -w Write the results to fitted coordinates to \
the output file\n");
           
   fprintf(stderr,"\nSimple program to fit two PDB files containing \
identical atoms but with\n");
   fprintf(stderr,"different coordinates.\n\n");
}

