/************************************************************************/
/**

   \file       pdblistss.c
   
   \version    V1.1
   \date       13.03.19
   \brief      List disulphide bonds
   
   \copyright  (c) UCL, Dr. Andrew C.R. Martin, 2014-2019
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
   Displays a list of disulphides based on calculated distances rather 
   than SSBOND or CONECT record data

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-   V1.0   20.07.15 Original   By: ACRM
-   V1.1   13.03.19 Fixed some buffer sizes

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

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
/* Ideal disulphide S-S length is 2.03A - we'll allow 2.25A             */
#define DISULPHIDE_CUTOFFSQ    5.0625


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);
void ListDisulphides(FILE *out, PDB *pdb);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//** 
   Main program for finding disulphides

-  20.07.15 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   FILE     *in     = stdin,
            *out    = stdout;
   int      natoms;
   PDB      *pdb;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];
   
   if(!ParseCmdLine(argc, argv, infile, outfile))
   {
      Usage();
      return(0);
   }

   if(!blOpenStdFiles(infile, outfile, &in, &out))
   {
      fprintf(stderr, "Error (pdblistss): Unable to open input or output \
file\n");
      return(1);
   }

   if((pdb = blReadPDBAtoms(in, &natoms))==NULL)
   {
      fprintf(stderr, "Error (pdblistss): No atoms read from PDB \
file, %s\n", infile);
      return(1);
   }

   ListDisulphides(out, pdb);

   /* Free up the memory for the PDB linked list                        */
   FREELIST(pdb, PDB);

   return(0);
}


/************************************************************************/
/*>void ListDisulphides(FILE *out, PDB *pdb)
   -----------------------------------------
*//**
   \param[in]   *out   Output file pointer
   \param[in]   *pdb   PDB linked list

   Does the actual work of finding and printing the disulphides

-  20.07.15   Original   By: ACRM
-  13.03.19   Increased resid size from 16 to 32
*/
void ListDisulphides(FILE *out, PDB *pdb)
{
   PDB *p, *q;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->resnam,"CYS",3) &&
         !strncmp(p->atnam, "SG  ",4))
      {
         for(q=p->next; q!=NULL; NEXT(q))
         {
            if(!strncmp(q->resnam,"CYS",3) &&
               !strncmp(q->atnam, "SG  ",4))
            {
               if(DISTSQ(p,q) < DISULPHIDE_CUTOFFSQ)
               {
                  char resid1[32],
                       resid2[32];
                  MAKERESID(resid1, p);
                  MAKERESID(resid2, q);

                  fprintf(out, "%6s Atom %5d : %6s Atom %5d : %.3f\n", 
                          resid1, p->atnum,
                          resid2, q->atnum,
                          DIST(p,q));
               }
            }
         }
      }
   }

}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ----------------------------------------------------------------------
*//**
   \param[in]   int    argc              Argument count
   \param[in]   char   **argv            Argument array
   \param[out]  char   *infile           Input filename (or blank string)
   \param[out]  char   *outfile          Output filename (or blank string)
   \return      BOOL                     Success

   Parse the command line

-  20.07.14 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
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

-  20.07.15 Original   By: ACRM
-  13.03.19 V1.1
*/
void Usage(void)
{
   fprintf(stderr,"\npdblistss V1.1 (c) 2015-2019 UCL, Dr. Andrew C.R. \
Martin\n");

   fprintf(stderr,"\nUsage: pdblistss [in.pdb [out.txt]]\n");

   fprintf(stderr,"\nDisplays a list of disulphides based on calculated \
distances rather\n");
   fprintf(stderr,"than SSBOND or CONECT record data.\n\n");
}


