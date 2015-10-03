/************************************************************************/
/**

   \file       pdbgetzone.c
   
   \version    V1.8
   \date       02.10.15
   \brief      Extract a numbered zone from a PDB file
   
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
-  V1.0   22.07.96  Original
-  V1.1   29.09.05  Added a command line option to prevent the residue
                    specifications being uppercased. This is required as 
                    some PDBs (eg 1rxsm) use upper and lower case letters 
                    as chain codes. (Tony Lewis) By: TL
-  V1.2   05.11.07  Improved error checking of command line
   V1.2.1 29.10.10  Fixed Bioplib bug where end of zone matched end of 
                    chain
-  V1.3   28.08.13  Modified for new ParseResSpec()
-  V1.4   22.07.14  Renamed deprecated functions with bl prefix.
                    Added doxygen annotation. By: CTP
-  V1.5   19.08.14  Fixed call to renamed function: 
                    blExtractZonePDBAsCopy() By: CTP
-  V1.6   06.11.14  Renamed from getpdb  By: ACRM
-  V1.7   13.02.15  Removed handling of -l option since there are now
                    too many PDB files with lower case chain labels for
                    it to make sense
-  V1.8   02.10.15  Added -w parameter

*************************************************************************/
/* Includes
*/
#include <stdio.h>
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
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *Zone1, char *Zone2,
                  char *infile, char *outfile, int *width);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   main routine

-  22.07.96 Original    By: ACRM
-  29.09.05 Modified for -l By: TL
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  19.08.14 Added AsCopy suffix for call to blExtractZonePDBAsCopy() 
            By: CTP
-  13.02.15 Removed -l handling - this is now the only option By: ACRM
-  02.10.15 Added -w (width) handling
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char InFile[MAXBUFF],
        OutFile[MAXBUFF],
        Zone1[MAXBUFF],
        Zone2[MAXBUFF],
        chain1[8],  chain2[8],
        insert1[8], insert2[8];
   int  res1,    res2,
        natom,   
        width = 0;
   PDB  *pdb;
   

   if(ParseCmdLine(argc, argv, Zone1, Zone2, InFile, OutFile, &width))
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         BOOL ParseResSpec1Result, ParseResSpec2Result;

         if((pdb=blReadPDB(in, &natom))==NULL)
         {
            fprintf(stderr,"pdbgetzone: No atoms read from PDB file\n");
            return(1);
         }
         
         ParseResSpec1Result = blParseResSpec(Zone1, chain1, 
                                              &res1, insert1);
         ParseResSpec2Result = blParseResSpec(Zone2, chain2, 
                                              &res2, insert2);

         if(!ParseResSpec1Result)
         {
            fprintf(stderr,"pdbgetzone: Illegal residue specification \
(%s)\n", Zone1);
            return(1);
         }
         if(!ParseResSpec2Result)
         {
            fprintf(stderr,"pdbgetzone: Illegal residue specification \
(%s)\n", Zone2);
            return(1);
         }

         if(width)
         {
            UpdateResRange(pdb, width,
                           chain1, &res1, insert1,
                           chain2, &res2, insert2);
         }

         if((pdb = blExtractZonePDBAsCopy(pdb, chain1, res1, insert1,
                                          chain2, res2, insert2))==NULL)
         {
            fprintf(stderr,"pdbgetzone: Zone not found (%s or %s)\n",
                    Zone1, Zone2);
            return(1);
         }
         
         blWritePDB(out,pdb);
      }
   }
   else
   {
      Usage();
   }
   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *Zone1, char *Zone2,
                     char *infile, char *outfile, int *width)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *Zone1       First end of zone
   \param[out]     *Zone2       Second end of zone
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *width       Amount to expand the range
   \return                      Success?

   Parse the command line
   
-  22.07.96 Original    By: ACRM
-  29.09.05 Added uppercaseresspec param and handling of -l  By: TL
-  05.11.07 Added first check that at least one parameter is on the
            command line  By: ACRM
-  13.02.15 -l option is now ignored - the program is always case 
            sensitive
-  02.10.15 Added -w and width parameter
*/
BOOL ParseCmdLine(int argc, char **argv, char *Zone1, char *Zone2,
                  char *infile, char *outfile, int *width)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   Zone1[0]  = Zone2[0]   = '\0';
   width     = 0;

   if(!argc)               /* 05.11.07 Added this                       */
   {
      return(FALSE);
   }
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         case 'l':
            fprintf(stderr, "-l option is now deprecated\n");
            break;
         case 'w':
            argc--; argv++;
            if(!argc) return(FALSE);
            if(!sscanf(argv[0], "%d", width)) return(FALSE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2, 3 or 4 arguments left               */
         if(argc < 2 || argc > 4)
            return(FALSE);
         
         /* Copy the first to Zone1 and second to Zone2                 */
         strcpy(Zone1, argv[0]);
         argc--;
         argv++;
         strcpy(Zone2, argv[0]);
         argc--;
         argv++;
         
         /* Copy the first to infile                                    */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
            argv++;
         }
         
         /* If there's another, copy it to outfile                      */
         if(argc)
         {
            strcpy(outfile, argv[0]);
            argc--;
            argv++;
         }
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}

/************************************************************************/
BOOL UpdateResRange(PDB *pdb, int width,
                    char *chain1, int *res1, char *insert1,
                    char *chain2, int *res2, char *insert2)
{
   PDBSTRUCT  *pdbs;
   PDBCHAIN   *pdbc;
   PDBRESIDUE *pdbr, *endr, *startr;
   char       chain[8], insert[8];
   int        resnum, i;
   
   
   if((pdbs = blAllocPDBStructure(pdb))==NULL)
      return(FALSE);

   for(pdbc = pdbs->chains; pdbc!=NULL; NEXT(pdbc))
   {
      if(CHAINMATCH(pdbc->chain, chain1))
      {
         for(pdbr = pdbc->residues; pdbr!=NULL; NEXT(pdbr))
         {
            if((pdbr->resnum == resnum1) && 
               INSERTMATCH(pdbr->insert, insert1))
            {
               startr = endr = pdbr;
               for(i=0; i<width; i++)
               {
                  if(startr!=NULL) startr = startr->prev;
               }
               if(startr == NULL) return(FALSE);
               
               *resnum1 = startr->resnum;
               strcpy(insert1, startr->insert);
               
               return(0);
            }
         }
         break;
      }
   }
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  22.07.14 V1.4 By: CTP
-  06.11.14 V1.6 By: ACRM
-  13.02.15 V1.7 By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\n");
   fprintf(stderr,"pdbgetzone V1.7 (c) 1996-2015, Dr. Andrew C.R. \
Martin, UCL.\n");
   fprintf(stderr,"                    Modified by Tony Lewis, \
UCL, 2005\n");

   fprintf(stderr,"\nUsage: pdbgetzone [-l] start end [in.pdb \
[out.pdb]]\n");
   fprintf(stderr,"       -l  Redundant - kept for backwards \
compatibility\n");
   fprintf(stderr,"\n");
   fprintf(stderr,"Start and end are resspec residue \
specifications:\n");
   blPrintResSpecHelp(stderr);
   fprintf(stderr,"\npdbgetzone extracts a specified zone from a PDB \
file writing it out in\n");
   fprintf(stderr,"PDB format. I/O is through standard input/output if \
filenames are\n");
   fprintf(stderr,"not specified.\n\n");
   fprintf(stderr,"Note that the residue specification is case \
sensitive. The -l option\n");
   fprintf(stderr,"used to be required for case sensitivity.\n\n");
}
