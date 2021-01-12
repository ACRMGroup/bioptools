/************************************************************************/
/**

   \file       pdbcount.c
   
   \version    V1.6.1
   \date       24.11.20
   \brief      Count residues and atoms in a PDB file
   
   \copyright  (c) Prof. Andrew C. R. Martin 1994-2020
   \author     Prof. Andrew C. R. Martin
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
-  V1.0   16.08.94 Original
-  V1.1   24.08.94 Changed to call OpenStdFiles()
-  V1.2   30.05.02 Changed PDB field from 'junk' to 'record_type'
-  V1.3   22.07.14 Renamed deprecated functions with bl prefix.
                   Added doxygen annotation. By: CTP
-  V1.4   06.11.14 Renamed from countpdb  By: ACRM
-  V1.5   12.03.15 Changed to allow multi-character chain names
-  V1.6   21.10.20 Added -c for by-chain calculation
-  V1.6.1 24.11.20 Fixed display of by-chain calculations

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
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
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *byChain);
void Usage(void);
void CountAndPrint(FILE *out, PDB *pdb, BOOL byChain);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for counting residues & atoms

-  16.08.94 Original    By: ACRM
-  24.08.94 Changed to call OpenStdFiles()
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   PDB  *pdb;
   BOOL byChain = FALSE;
   int  natom;
        
   if(ParseCmdLine(argc, argv, infile, outfile, &byChain))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natom))==NULL)
         {
            fprintf(stderr,"No atoms read from input file\n");
         }
         else
         {
            CountAndPrint(out, pdb, byChain);
         }
      }
      else
      {
         Usage();
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *byChain)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \param[out]     *outfile    Output filename (or blank string)
   \param[out]     *byChain    Calculate for each chain separately
   \return                     Success

   Parse the command line

-  16.08.94 Original    By: ACRM
-  21.10.20 Added byChain
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *byChain)
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
            *byChain = TRUE;
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

   Print a usage message

-  16.08.94 Original    By: ACRM
-  22.07.14 V1.3 By: CTP
-  06.11.14 V1.4 By: ACRM
-  12.03.15 V1.5
-  21.10.20 V1.6
-  24.11.20 V1.6.1
*/
void Usage(void)
{
   fprintf(stderr,"\npdbcount V1.6.1 (c) 1994-2020 Prof. Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"\nUsage: pdbcount [-c] [in.pdb [out.txt]]\n");
   fprintf(stderr,"       -c Calculate for each chain separately\n\n");
   
   fprintf(stderr,"If files are not specified, stdin and stdout are \
used.\n");
   fprintf(stderr,"Counts chains, residues & atoms in a PDB file.\n\n");
}


/************************************************************************/
/*>void CountAndPrint(FILE *out, PDB *pdb, BOOL byChain)
   --------------------------------------------------------
*//**

   Does the actual work of counting and printing

-  16.08.94 Original    By: ACRM
-  12.03.15 Changed to allow multi-character chain names
-  21.10.20 Now does the printing too - Added byChain paramater
*/
void CountAndPrint(FILE *out, PDB *pdb, BOOL byChain)
{
   PDB *p;
   char LastChain[9],
        LastIns = '-';
   int  LastRes = -9999,
        nchainTotal = 0,
        nresTotal   = 0,
        natomTotal  = 0,
        nhydTotal   = 0,
        nhetTotal   = 0,
        nres        = 0,
        natom       = 0,
        nhyd        = 0;
   
   strcpy(LastChain, "-");

   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->record_type,"ATOM  ",6))
      {
         if(!CHAINMATCH(p->chain, LastChain))
         {
            /* Chain has changed

               Print information for the previous chain 
            */
            if(nres && byChain)
            {
               fprintf(out,"Chain: %s Residues: %d Atoms: %d \
Hydrogens: %d\n", LastChain, nres, natom, nhyd);
            }

            /* Update total counts                                      */
            (nchainTotal)++;
            nresTotal  += nres;
            natomTotal += natom;
            nhydTotal  += nhyd;

            /* Reset counts for new chain                               */
            nres  = 1;
            natom = 0;
            nhyd  = 0;

            /* Update current chain and residue info                    */
            strncpy(LastChain, p->chain, 8);
            LastRes   = p->resnum;
            LastIns   = p->insert[0];
         }
         else if((p->insert[0] != LastIns) ||
                 (p->resnum    != LastRes))
         {
            /* Residue has changed                                      */
            nres++;
         
            LastRes   = p->resnum;
            LastIns   = p->insert[0];
         }

         /* Update hydrogen and atom counts                             */
         if(p->atnam[0] == 'H') 
            nhyd++;
      
         natom++;
      }
      else if(!strncmp(p->record_type,"HETATM",6))
      {
         nhetTotal++;
      }
   }
   
   /* Print information for the last chain                              */
   if(nres && byChain)
   {
      fprintf(out,"Chain: %s Residues: %d Atoms: %d Hydrogens: %d\n",
              LastChain, nres, natom, nhyd);
   }
   
   /* Add data for last chain                                           */
   nresTotal  += nres;
   natomTotal += natom;
   nhydTotal  += nhyd;

   /* Print total counts                                                */
   if(byChain)
   {
      fprintf(out,"NChains: %d TotResidues: %d TotAtoms: %d \
TotHetAtoms: %d TotHydrogens: %d\n",
              nchainTotal, nresTotal, natomTotal, nhetTotal, nhydTotal);
   }
   else
   {
      fprintf(out,"Chains: %d Residues: %d Atoms: %d Het Atoms: %d \
Total Hydrogens: %d\n",
              nchainTotal, nresTotal, natomTotal, nhetTotal, nhydTotal);
   }
}
