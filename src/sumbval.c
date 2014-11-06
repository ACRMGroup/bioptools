/************************************************************************/
/**

   \file       sumbval.c
   
   \version    V1.3
   \date       22.07.14
   \brief      Sum B-vals over each residue and replace with the
               summed or average value
   
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
-  V1.0  05.07.94 Original
-  V1.1  06.07.94 Fixed potential /0 bug
-  V1.2  24.08.94 Changed to call OpenStdFiles()
-  V1.3  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/MathUtil.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *average, BOOL *sidechain, BOOL *quiet);
void Usage(void);
void SumBVals(PDB *pdb, BOOL average, BOOL sidechain, BOOL quiet);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for summing B-values

-  05.07.94 Original    By: ACRM
-  24.08.94 Changed to call OpenStdFiles()
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE *in       = stdin,
        *out      = stdout;
   BOOL average   = FALSE,
        sidechain = FALSE,
        quiet     = FALSE;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   int  natoms;
   PDB  *pdb;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &average, &sidechain,
                   &quiet))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDB(in,&natoms)) != NULL)
         {
            SumBVals(pdb, average, sidechain, quiet);
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
                     BOOL *average, BOOL *sidechain, BOOL *quiet)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *average     Average the b-vals
   \param[out]     *sidechain   Separate s/c and m/c (N,CA,C,O)
   \param[out]     *quiet       Do not display values
   \return                     Success?

   Parse the command line
   
-  05.07.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *average, BOOL *sidechain, BOOL *quiet)
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
         case 'a':
            *average = TRUE;
            break;
         case 's':
            *sidechain = TRUE;
            break;
         case 'q':
            *quiet = TRUE;
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

-  05.07.94 Original    By: ACRM
-  06.07.94 V1.1
-  24.08.94 V1.2
-  22.07.14 V1.3 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\nSumBVal V1.3 (c) 1994-2014, Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: sumbval [-a] [-s] [-q] [<in.pdb>] \
[<out.pdb>]\n");
   fprintf(stderr,"                -a Average over the residues\n");
   fprintf(stderr,"                -s Separate s/c and m/c\n");
   fprintf(stderr,"                -q Do not display overall mean and \
standard deviation\n\n");
   fprintf(stderr,"Sums the b-values over each residue and places the \
summed values in the\n");
   fprintf(stderr,"b-value column. Averaging causes averages rather than \
sums to be used\n");
   fprintf(stderr,"while separation causes the mainchain (N,CA,C,O) for \
each residue to be\n");
   fprintf(stderr,"treated separately from the sidechain.\n\n");
}

/************************************************************************/
/*>void SumBVals(PDB *pdb, BOOL average, BOOL sidechain, BOOL quiet)
   -----------------------------------------------------------------
*//**

   \param[in,out]  *pdb       PDB linked list
   \param[in]      average    Do averaging
   \param[in]      sidechain  Treat m/c and s/c separately
   \param[in]      quiet      Do not display overall means and SDs

   Does the actual work of summing and averaging B-values in each residue.

-  05.07.94 Original    By: ACRM
-  06.07.94 Added check before divides
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
void SumBVals(PDB *pdb, BOOL average, BOOL sidechain, BOOL quiet)
{
   PDB  *p, *start, *end;
   REAL SumSC,
        SumMC,
        mean,
        sd,
        Sx[3],
        SxSq[3];
   int  CountSC,
        CountMC,
        NValues[3];

   /* Initialise valiables for calculating overall mean and SD          */
   Sx[0]      = Sx[1]      = Sx[2]      = (REAL)0.0;
   SxSq[0]    = SxSq[1]    = SxSq[2]    = (REAL)0.0;
   NValues[0] = NValues[1] = NValues[2] = 0;
   
   /* Step through the PDB linked list                                  */
   for(start=pdb; start!=NULL; start=end)
   {
      /* Find start of next residue                                     */
      end = blFindEndPDB(start);

      /* Zero the totals and counts                                     */
      SumSC   = SumMC   = (REAL)0.0;
      CountSC = CountMC = 0;

      /* Step through atoms in this residue                             */
      for(p=start; p!=end; NEXT(p))
      {
         if(!strncmp(p->atnam,"N   ",4) ||
            !strncmp(p->atnam,"CA  ",4) ||
            !strncmp(p->atnam,"C   ",4) ||
            !strncmp(p->atnam,"O   ",4))
         {
            CountMC++;
            SumMC += p->bval;
            blCalcExtSD(p->bval, 0, &(Sx[0]), &(SxSq[0]), &(NValues[0]),
                      &mean, &sd);
         }
         else
         {
            CountSC++;
            SumSC += p->bval;
            blCalcExtSD(p->bval, 0, &(Sx[1]), &(SxSq[1]), &(NValues[1]),
                      &mean, &sd); 
         }

         blCalcExtSD(p->bval, 0, &(Sx[2]), &(SxSq[2]), &(NValues[2]),
                   &mean, &sd);
      }
      
      /* Calculate averages if required                                 */
      if(average)
      {
         if(sidechain)
         {
            /* Calc average separately for m/c and s/c                  */
            if(CountMC) SumMC /= CountMC;
            if(CountSC) SumSC /= CountSC;
         }
         else
         {
            /* Calc average for m/c and s/c together and copy into both
               the SumMC and SumSC variables
            */
            if(CountMC + CountSC)
               mean = (SumMC + SumSC)/(CountMC + CountSC);
            SumMC = SumSC = mean;
         }
      }
      else                        /* Not averaging                      */
      {
         if(!sidechain)
         {
            /* Not treating sidechain separately, so add the values
               together and equate them
            */
            SumMC += SumSC;
            SumSC = SumMC;
         }
      }

      /* Step through the residue again, filling in the values          */
      for(p=start; p!=end; NEXT(p))
      {
         if(!strncmp(p->atnam,"N   ",4) ||
            !strncmp(p->atnam,"CA  ",4) ||
            !strncmp(p->atnam,"C   ",4) ||
            !strncmp(p->atnam,"O   ",4))
         {
            p->bval = SumMC;
         }
         else
         {
            p->bval = SumSC;
         }
      }
   }

   /* Print output values to stderr                                     */
   if(!quiet)
   {
      blCalcExtSD((REAL)0.0, 1, &(Sx[0]), &(SxSq[0]), &(NValues[0]), 
                &mean, &sd);
      fprintf(stderr,"Mean B-value over backbone (N,CA,C,O) = %6.3f, \
SD = %6.3f (%d atoms)\n",mean,sd,NValues[0]);
      
      blCalcExtSD((REAL)0.0, 1, &(Sx[1]), &(SxSq[1]), &(NValues[1]), 
                &mean, &sd);
      fprintf(stderr,"Mean B-value over sidechains          = %6.3f, \
SD = %6.3f (%d atoms)\n",mean,sd,NValues[1]);

      blCalcExtSD((REAL)0.0, 1, &(Sx[2]), &(SxSq[2]), &(NValues[2]), 
                &mean, &sd);
      fprintf(stderr,"Mean B-value over all atoms           = %6.3f, \
SD = %6.3f (%d atoms)\n",mean,sd,NValues[2]);
   }
}

