/*************************************************************************

   Program:    avbr
   File:       avbr.c
   
   Version:    V1.0
   Date:       07.10.94
   Function:   Calc means and SDs of BValues by residue type
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0372) 275775
   EMail:      INTERNET: amartin@scitec.adsp.sub.org
                         martin@bsm.bioc.ucl.ac.uk
               UUCP:     ...{uunet|rutgers}!cbmehq!cbmuk!scitec!amartin
               JANET:    martin@uk.ac.ucl.bioc.bsm
               
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
   V1.0  07.10.94 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "bioplib/MathUtil.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXRES   24
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/
static char *sTypes[] = 
{
   "ALA ", "CYS ", "ASP ", "GLU ", "PHE ", "GLY ", "HIS ", "ILE ",
   "LYS ", "LEU ", "MET ", "ASN ", "PRO ", "GLN ", "ARG ", "SER ",
   "THR ", "VAL ", "TRP ", "TYR ", "UNK ", "GLX ", "ASX ", "PCA ",
   NULL
} ;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *FindMax, REAL *MaxVal, BOOL *Normalise, 
                  int *NBin);
void DoMeanSD(FILE *out, PDB *pdb);
BOOL DoBarchart(FILE *out, PDB *pdb, BOOL FindMax, REAL MaxVal, 
                BOOL Normalise, int NBin);
BOOL DoBarsForRes(FILE *out, PDB *pdb, char *resnam, BOOL FindMax, 
                  REAL MaxVal, BOOL Normalise, int NBin);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for average B-value by residue.

   07.10.94 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   PDB  *pdb;
   int  natoms,
        NBin;
   char InFile[MAXBUFF],
        OutFile[MAXBUFF];
   REAL MaxVal;
   BOOL FindMax   = TRUE,
        Normalise = FALSE;
   
   if(ParseCmdLine(argc, argv, InFile, OutFile, &FindMax, &MaxVal, 
                   &Normalise, &NBin))
   {
      if(OpenStdFiles(InFile, OutFile, &in, &out))
      {
         if((pdb = ReadPDB(in, &natoms)) != NULL)
         {
            DoMeanSD(out, pdb);
            if(!DoBarchart(out, pdb, FindMax, MaxVal, Normalise, NBin))
            {
               return(1);
            }
         }
         else
         {
            fprintf(stderr,"abvr: No atoms read from PDB file\n");
            return(1);
         }
      }
      else
      {
         return(1);
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
   Prints usage message

   07.10.94 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nabvr V1.0 (c) 1994, Dr. Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: abvr [-n] [-m maxval] [-b nbin] [in.pdb] \
[output.txt]\n");
   fprintf(stderr,"       -n  Normalise output bars (sum will be 1.0)\n");
   fprintf(stderr,"       -m  Specify max value on x-axis\n");
   fprintf(stderr,"       -b  Specify number of bins (default: 10)\n\n");
   fprintf(stderr,"Calculates means and standard deviations for \
B-values per residue.\n");
   fprintf(stderr,"I/O through standard input/output if files not \
specified.\n\n");
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile
                     BOOL *FindMax, REAL *MaxVal, BOOL *Normalise, 
                     int *NBin)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            BOOL   *FindMax     Find max value automatically?
            REAL   *MaxVal      Supplied max value
            BOOL   *Normalise   Normalise bars?
            int    *NBin        Number of bars (bins)
   Returns: BOOL                Success?

   Parse the command line
   
   07.10.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *FindMax, REAL *MaxVal, BOOL *Normalise, int *NBin)
{
   argc--;
   argv++;

   infile[0]  = outfile[0] = '\0';
   *FindMax   = TRUE;
   *MaxVal    = (REAL)0.0;
   *Normalise = FALSE;
   *NBin      = 10;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'n':
            *Normalise = TRUE;
            break;
         case 'm':
            *FindMax = FALSE;
            argc--;
            argv++;
            if(argc < 0)
               return(FALSE);
            if(!sscanf(argv[0],"%lf",MaxVal))
               return(FALSE);
            break;
         case 'b':
            argc--;
            argv++;
            if(argc < 0)
               return(FALSE);
            if(!sscanf(argv[0],"%d",NBin))
               return(FALSE);
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
/*>void DoMeanSD(FILE *out, PDB *pdb)
   ----------------------------------
   Calculates means & SDs on a per-residue basis.

   07.10.94 Original    By: ACRM
*/
void DoMeanSD(FILE *out, PDB *pdb)
{
   PDB *p;
   REAL Sx[MAXRES],
        SxSq[MAXRES],
        mean,
        sd;
   int  resnum,
        NValues[MAXRES];

   /* Print header                                                      */
   fprintf(out,"Means and standard deviations\n");
   fprintf(out,"=============================\n");

   /* Initialise variable arrays                                        */
   for(resnum=0; sTypes[resnum] != NULL; resnum++)
   {
      Sx[resnum]      = (REAL)0.0;
      SxSq[resnum]    = (REAL)0.0;
      NValues[resnum] = 0;
   }
   
   /* Sample the value                                                  */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(resnum=0; sTypes[resnum] != NULL; resnum++)
      {
         if(!strncmp(p->resnam,sTypes[resnum],4))
            break;
      }
      if(sTypes[resnum] == NULL)
      {
         fprintf(stderr,"Unknown residue type: %s\n",p->resnam);
         return;
      }
      
      CalcExtSD(p->bval, 0, &(Sx[resnum]), &(SxSq[resnum]),
                &(NValues[resnum]), &mean, &sd);
   }

   /* Display the values                                                */
   for(resnum=0; sTypes[resnum] != NULL; resnum++)
   {
      CalcExtSD((REAL)0.0, 1, &(Sx[resnum]), &(SxSq[resnum]),
                &(NValues[resnum]), &mean, &sd);
      fprintf(out,"%4s Mean: %f SD: %f\n",
              sTypes[resnum], mean, sd);
   }
}

/************************************************************************/
BOOL DoBarchart(FILE *out, PDB *pdb, BOOL FindMax, REAL MaxVal, 
                BOOL Normalise, int NBin)
{
   int resnum;
   
   /* Print header                                                      */
   fprintf(out,"\n\nBarchart\n");
   fprintf(out,"========\n");

   fprintf(out,"Res  MaxVal Bars...\n");
   fprintf(out,"-------------------\n");
   for(resnum=0; sTypes[resnum] != NULL; resnum++)
      if(!DoBarsForRes(out, pdb, sTypes[resnum], FindMax, MaxVal, 
                       Normalise, NBin))
         return(FALSE);
   return(TRUE);
}

/************************************************************************/
BOOL DoBarsForRes(FILE *out, PDB *pdb, char *resnam, BOOL FindMax, 
                  REAL MaxVal, BOOL Normalise, int NBin)
{
   PDB *p;
   int  nvalues = 0,
        i;
   static int *bins = NULL;

   if(bins==NULL)
   {
      /* Allocate memory for bins                                       */
      if((bins = (int *)malloc(NBin * sizeof(int)))==NULL)
      {
         fprintf(stderr,"abvr: No memory for barchart bins\n");
         return(FALSE);
      }
   }
   
   /* Zero the bins                                                     */
   for(i=0; i<NBin; i++)
      bins[i] = 0;
   
   if(FindMax)
   {
      /* Find the max B-value                                           */
      MaxVal = (REAL)(-1000000.0);
      
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(!strncmp(p->resnam,resnam,4) && (p->bval > MaxVal))
            MaxVal = p->bval;
      }
   }
   
   /* Run through again and place values into bins                      */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->resnam,resnam,4))
      {
         i = (int)(NBin * p->bval / MaxVal);
         if(i==NBin) i--;
         
         (bins[i])++;
         nvalues++;
      }
   }

   /* Print the bins                                                    */
   fprintf(out,"%4s %6.3f ",resnam, MaxVal);
   for(i=0; i<NBin; i++)
   {
      if(Normalise)
         fprintf(out,"%5.4f ",(REAL)bins[i]/(REAL)nvalues);
      else
         fprintf(out,"%5d ",bins[i]);
   }
   fprintf(out,"\n");

   return(TRUE);
}
