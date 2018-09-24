/************************************************************************/
/**

   \file       pdbhphob.c
   
   \version    V1.0
   \date       24.09.18
   \brief      Patches hydrophobicity data into B-value column of PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin 2018
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
   Patches hydrophobicity data into B-value column of PDB file

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

-  V1.0  24.09.18 Original based on older FORTRan code

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
#define DEF_HPHOBFILE "consensus.hpb"
#define MAXBUFF       512
#define MAXAATYPES     40

typedef struct _hphob
{
   char resnam[8];
   REAL value;
}  HPHOB;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *hphobFile);
BOOL ReadHPhobFile(char *filename, HPHOB *hphob, int maxhphob);
void Usage(void);
void PatchHPhob(PDB *pdb, HPHOB *hphob);
REAL FindHPhobValue(HPHOB *hphob, char *resnam);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program for patching hydrophobicity values into the B-value
   column of a PDB file.

-  24.09.18 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE    *in  = stdin,
           *out = stdout;
   char    inFile[MAXBUFF],
           outFile[MAXBUFF],
           hphobFile[MAXBUFF];
   int     natoms  = 0,
           retval  = 0;
   PDB     *pdb    = NULL;
   HPHOB   hphob[MAXAATYPES];

   
   if(ParseCmdLine(argc, argv, inFile, outFile, hphobFile))
   {
      if(ReadHPhobFile(hphobFile, hphob, MAXAATYPES))
      {
         if(blOpenStdFiles(inFile, outFile, &in, &out))
         {
            if((pdb=blReadPDB(in, &natoms))!=NULL)
            {
               PatchHPhob(pdb, hphob);
               blWritePDB(out, pdb);
            }
            else
            {
               fprintf(stderr,"pdbhphob: Error - no atoms read from \
PDB file\n");
               return(1);
            }
         }
         else
         {
            fprintf(stderr,"pdbhphob: Error - unable to open input or \
output file\n");
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
   }
   
   return(retval);
}


/************************************************************************/
/*>void PatchHPhob(PDB *pdb, HPHOB *hphob)
   ---------------------------------------
*//**
   \param[in,out] *pdb      PDB linked list
   \param[in]     *hphob    Array of hydrophobicity data

   Does the actual work of patching hydrophobicity values into the 
   temperature factor data

-  24.09.18 Original   By: ACRM
*/
void PatchHPhob(PDB *pdb, HPHOB *hphob)
{
   PDB *p,
       *nextRes;

   for(p=pdb; p!=NULL; p=nextRes)
   {
      REAL value;
      PDB  *q;
      
      nextRes = blFindNextResidue(p);
      value = FindHPhobValue(hphob, p->resnam);
      if(value > 9998.0)
         value = 0.0;
      for(q=p; q!=nextRes; NEXT(q))
      {
         q->bval = value;
      }
   }
}


/************************************************************************/
/*>REAL FindHPhobValue(HPHOB *hphob, char *resnam)
   -----------------------------------------------
*//**
   \param[in]  *hphob   Hydrophobicity data
   \param[in]  *resnam  Residue name for which we need the hydrophobicity
   \return              The hydrophobicity value or 9999.0 if the residue
                        name is not found

   Finds the hydrophobicity value for a given residue name from the
   data read from a hydrophobcity file

-  24.09.18 Original   By: ACRM
*/
REAL FindHPhobValue(HPHOB *hphob, char *resnam)
{
   int i=0;
   for(i=0; hphob[i].value < 9998.0; i++)
   {
      if(!strncmp(hphob[i].resnam, resnam, 3))
         return(hphob[i].value);
   }

   return(9999.0);
}


/************************************************************************/
/*>BOOL ReadHPhobFile(char *filename, HPHOB *hphob, int maxhphob)
   --------------------------------------------------------------
*//**
   \param[in]  *filename  Hydrophobicity file name
   \param[out] *hphob     The array of hydrophobicity data
   \param[in]  maxhphob   The max size of the hydrophibicity array

   Reads a hydrophobicity file using specified name or looking in
   $DATADIR

-  24.09.18 Original   By: ACRM
*/
BOOL ReadHPhobFile(char *filename, HPHOB *hphob, int maxhphob)
{
   BOOL noenv;
   FILE *fp;
   
   if((fp=blOpenFile(filename, "DATADIR", "r", &noenv))!=NULL)
   {
      char buffer[MAXBUFF];
      int linenum = 0;

      /* Initialize with dummy values                                   */
      for(linenum=0; linenum<maxhphob; linenum++)
      {
         hphob[linenum].resnam[0] = '\0';
         hphob[linenum].value     = 9999.0;
      }
      linenum = 0;
      
      /* Read the file skipping the first line                          */
      while(fgets(buffer, MAXBUFF, fp))
      {
         if(linenum && (linenum<maxhphob))
         {
            sscanf(buffer, "%s %lf",
                   hphob[linenum-1].resnam,
                   &(hphob[linenum-1].value));
         }
         linenum++;
      }
      
      fclose(fp);
   }
   else
   {
      fprintf(stderr,"pdbhphob: Error - unable to read hydrophobicity \
file\n");
      if(noenv)
      {
         fprintf(stderr, "          Environment variable, DATADIR not \
set.\n");
      }
         
      return(FALSE);
   }
   return(TRUE);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     char *hphobfile)
   ---------------------------------------------------------------------
*//**
   \param[in]     argc         Argument count
   \param[in]     **argv       Argument array
   \param[out]    *infile      Input file (or blank string)
   \param[out]    *outfile     Output file (or blank string)
   \param[out]    *hphobfile   Hydrophobicity filename
   \return                     Success?

   Parse the command line
   
-  24.09.18 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char * hphobfile)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';

   strcpy(hphobfile, DEF_HPHOBFILE);
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd':
            argc--;
            argv++;
            if(argc < 1)
               return(FALSE);
            strcpy(hphobfile, argv[0]);
            break;
         case 'h':
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

-  24.09.18 Original   By: ACRM
*/
void Usage(void)
{
   printf("\npdbhphob V1.0 (c) 2018 UCL, Dr. Andrew C.R. Martin\n");

   printf("\nUsage: pdbhphob [-d datafile] [in.pdb [out.txt]]\n");
   printf("       -d Specify hydrophobicity data file [Default: %s]\n",
          DEF_HPHOBFILE);

   printf("\npdbhphob takes a PDB file and patches the residue \
hydrophobicity values\n");
   printf("into the B-value (temperature factor) column of the PDB data, \
writing\n");
   printf("a new PDB file. Colouring by temperature factor in a \
molecular graphics\n");
   printf("program will then colour by residue hydrophobicity.\n\n");

   printf("The hydrophobicity file format consists of a required \
comment line which\n");
   printf("is simply skipped by the code, followed by (normally) 20 \
lines each of\n");
   printf("which contains two fields: the amino acid name (3-letter \
code in\n");
   printf("capitals), and the hydrophobicity value.\n\n");
}

