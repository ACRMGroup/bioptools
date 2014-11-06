/************************************************************************/
/**

   \file       countpdb.c
   
   \version    V1.3
   \date       22.07.14
   \brief      Count residues and atoms in a PDB file
   
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
-  V1.0  16.08.94 Original
-  V1.1  24.08.94 Changed to call OpenStdFiles()
-  V1.2  30.05.02 Changed PDB field from 'junk' to 'record_type'
-  V1.3  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

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
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);
void DoCount(PDB *pdb, int *nchain, int *nres, int *natom, int *nhyd,
             int *nhet);

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
   int  nchain, nres, natom, nhyd, nhet;
        
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natom))==NULL)
         {
            fprintf(stderr,"No atoms read from input file\n");
         }
         else
         {
            DoCount(pdb, &nchain, &nres, &natom, &nhyd, &nhet);
            fprintf(out,"Chains: %d Residues: %d Atoms: %d Het Atoms: %d \
Total Hydrogens: %d\n", nchain, nres, natom, nhet, nhyd);
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \param[out]     *outfile    Output filename (or blank string)
   \return                     Success

   Parse the command line

-  16.08.94 Original    By: ACRM
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

   Print a usage message

-  16.08.94 Original    By: ACRM
-  22.07.14 V1.3 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\nCountPDB V1.3 (c) 1994-2014 Dr. Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Freely distributable if no profit is made\n\n");
   fprintf(stderr,"Usage: countpdb [<infile>] [<outfile>]\n\n");
   fprintf(stderr,"If files are not specified, stdin and stdout are \
used.\n");
   fprintf(stderr,"Counts chains, residues & atoms in a PDB file.\n\n");
}

/************************************************************************/
/*>void DoCount(PDB *pdb, int *nchain, int *nres, int *natom, int *nhyd,
                int *nhet)
   ---------------------------------------------------------------------
*//**

   Does the actual work of counting 

-  16.08.94 Original    By: ACRM
*/
void DoCount(PDB *pdb, int *nchain, int *nres, int *natom, int *nhyd,
             int *nhet)
{
   PDB *p;
   char LastChain,
        LastIns;
   int  LastRes;

   LastChain = '-';
   LastIns   = '-';
   LastRes   = -999;

   *nchain   = 0;
   *nres     = 0;
   *natom    = 0;
   *nhyd     = 0;
   *nhet     = 0;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(!strncmp(p->record_type,"ATOM  ",6))
         (*natom)++;
      else if(!strncmp(p->record_type,"HETATM",6))
         (*nhet)++;
      
      if(p->atnam[0] == 'H') 
         (*nhyd)++;
      
      if(p->chain[0] != LastChain)
      {
         /* Chain has changed & so, by definition has the residue       */
         if(!strncmp(p->record_type,"ATOM  ",6))
         {
            (*nchain)++;
            (*nres)++;
         }

         LastChain = p->chain[0];
         LastRes   = p->resnum;
         LastIns   = p->insert[0];
      }
      else if((p->insert[0] != LastIns) ||
              (p->resnum    != LastRes))
      {
         if(!strncmp(p->record_type,"ATOM  ",6))
            (*nres)++;
         
         LastRes   = p->resnum;
         LastIns   = p->insert[0];
      }
   }
}
