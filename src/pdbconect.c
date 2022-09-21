/************************************************************************/
/**

   \file       pdbconect.c
   
   \version    V1.1
   \date       21.09.22
   \brief      Rebuild CONECT records for a PDB file
   
   \copyright  (c) Prof. Andrew C. R. Martin 2015-2022
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
-  V1.0  26.02.15 Original
-  V1.1  21.09.22 Added -m option to merge chains that are connected

*************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 256
#define DEF_TOL 0.2

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *tol, BOOL *merge);
void MergeConnectedChains(PDB *pdb);
void RelabelChain(PDB *pdb, char *oldLabel, char *newLabel);
void RenumberHetResidues(PDB *pdb);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program

-  26.02.15 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE     *in      = stdin,
            *out     = stdout;
   WHOLEPDB *wpdb;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];
   REAL     tol = DEF_TOL;
   BOOL     merge = FALSE;

   if(ParseCmdLine(argc, argv, infile, outfile, &tol, &merge))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb=blReadWholePDB(in))!=NULL)
         {
            blBuildConectData(wpdb->pdb, tol);
            if(merge)
               MergeConnectedChains(wpdb->pdb);
            blWriteWholePDB(out, wpdb);
         }
         else
         {
            fprintf(stderr,"No atoms read from PDB file\n");
            return(1);
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
/*>void Usage(void)
   ----------------
*//**

-  26.02.15 Original    By: ACRM
-  21.09.22 Added -m
*/
void Usage(void)
{
   fprintf(stderr,"\npdbconect V1.1  (c) 2015-22 UCL, Andrew C.R. \
Martin\n");
   fprintf(stderr,"Usage: pdbconect [-t x][-m] [<input.pdb> \
[<output.pdb>]]\n");
   fprintf(stderr,"       -t specify tolerance [Default: %.1f]\n", 
           DEF_TOL);
   fprintf(stderr,"       -m merge chains connected via CONECTs\n");

   fprintf(stderr,"\nGenerates CONECT records for a PDB file from the \
covalent radii of the\n");
   fprintf(stderr,"elements involved. Existing CONECT records are \
discarded first\n");
   fprintf(stderr,"I/O is to stdin/stdout if not specified\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     REAL *tol, BOOL *merge)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *tol         Covalent bond tolerance
   \return                      Success?

   Parse the command line
   
-  26.02.15 Original    By: ACRM
-  21.09.22 Added -m / merge
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *tol, BOOL *merge)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *tol      = DEF_TOL;
   *merge    = FALSE;
   
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 't':
            argc--;
            argv++;
            if(argc)
            {
               if(!sscanf(argv[0], "%lf", tol))
                  return(FALSE);
            }
            break;
         case 'm':
            *merge = TRUE;
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
/*>void MergeConnectedChains(PDB *pdb)
   -----------------------------------
*//**
   \param[in,out]  *pdb         PDB linked list

   Examine the CONECT data and relabel connected chains so they all 
   match the first one.
   
-  21.09.22 Original    By: ACRM
*/
void MergeConnectedChains(PDB *pdb)
{
   PDB *p;
   BOOL changed = TRUE;
   
   while(changed)
   {
      changed = FALSE;
      
      /* Step through the PDB file                                      */
      for(p=pdb; p!=NULL; NEXT(p))
      {
         int i;
         /* Step through the CONECTs for this atom                      */
         for(i=0; i<p->nConect; i++)
         {
            PDB *q = p->conect[i];
            if(!PDBCHAINMATCH(p,q))
            {
               changed = TRUE;
               RelabelChain(pdb, q->chain, p->chain);
            }
         }
      }
   }
   RenumberHetResidues(pdb);
}

/************************************************************************/
/*>void RelabelChain(PDB *pdb, char *oldLabel, char *newLabel)
   -----------------------------------------------------------
*//**
   \param[in,out]  *pdb         PDB linked list
   \param[in]      *oldLabel    Old label to be replaced
   \param[in]      *newLabel    New label

   Relabels a specified chain
   
-  21.09.22 Original    By: ACRM
*/
void RelabelChain(PDB *pdb, char *oldLabel, char *newLabel)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(CHAINMATCH(p->chain, oldLabel))
      {
         strcpy(p->chain, newLabel);
      }
   }
}

/************************************************************************/
/*>void RenumberHetResidues(PDB *pdb)
   ----------------------------------
*//**
   \param[in,out]  *pdb         PDB linked list

   Renumber HET residues such that they follow on from the proceding
   chain
   
-  21.09.22 Original    By: ACRM
*/
void RenumberHetResidues(PDB *pdb)
{
   int prevRes = pdb->resnum,
       thisRes = pdb->resnum,
       newRes = 0;
   PDB *p;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      /* If the residue has changed, set prevRes and update thisRes     */
      if(p->resnum != thisRes)
      {
         prevRes = thisRes;
         thisRes = p->resnum;
         newRes  = prevRes+1;
      }

      /* If it's a HET residue, renumber it based on the number after 
         the last one in the PDB ATOM chain
      */
      if(!strncmp(p->record_type, "HETATM", 6))
      {
         PDB *q;
         int startRes = p->resnum;
         
         for(q=p; (q!=NULL) && (q->resnum == startRes); NEXT(q))
         {
            q->resnum = newRes;
         }
         thisRes=newRes;
      }
   }
}

