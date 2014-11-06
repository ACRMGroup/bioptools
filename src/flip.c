/************************************************************************/
/**

   \file       flip.c
   
   \version    V1.2
   \date       19.08.14
   \brief      Standardise equivalent atom labelling
   
   \copyright  (c) Dr. Andrew C. R. Martin 1996-2014
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
   
-  V1.1   22.07.14 Renamed deprecated functions with bl prefix.
                   Added doxygen annotation. By: CTP
-  V1.2   19.08.14 Removed unused variable in DoFlipping() By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/angle.h"


/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

typedef struct
{
   char *resnam,
        *atom1,
        *atom2,
        *atom3,
        *atom4,
        *atom4b,
        *connect4,
        *connect4b;
} TORSION;


/************************************************************************/
/* Globals
*/


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *verbose, BOOL *quiet);
void Usage(void);
void DoFlipping(PDB *pdb, BOOL verbose, BOOL quiet);
void DoAFlip(PDB *atom4, PDB *atom4b, PDB *connect4, PDB *connect4b);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program to flip sidechain equivalent atom naming

-  08.11.96 Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   PDB  *pdb;
   int  natoms;
   BOOL verbose, quiet;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &verbose, &quiet))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDB(in,&natoms)) != NULL)
         {
            DoFlipping(pdb,verbose,quiet);
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
                     BOOL *verbose, BOOL *quiet)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *verbose     Report each flip made?
   \param[out]     *quiet       Don't report skipped residues
   \return                      Success?

   Parse the command line
   
-  08.11.96 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *verbose, BOOL *quiet)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *verbose  = FALSE;
   *quiet    = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            *verbose = TRUE;
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

-  08.11.96 Original   By: ACRM
-  22.07.14 V1.1 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\nflip V1.1 (c) 2014 Dr. Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"\nUsage: flip [in.pdb [out.pdb]]\n");

   fprintf(stderr,"\nFlip is a rather crude and simple program for \
correcting the atom\n");
   fprintf(stderr,"naming of equivalent atoms about freely rotable \
bonds. Currently\n");
   fprintf(stderr,"it handles ARG, ASP, GLU, PHE, and TYR and \
assumes that\n");
   fprintf(stderr,"the connectivity is correct (e.g. in PHE, CE1 is \
connected to\n");
   fprintf(stderr,"CD1 and CE2 is connected to CD2). A more \
sophisticated version\n");
   fprintf(stderr,"should also check connectivity of atom names and \
should therefore\n");
   fprintf(stderr,"also handle ILE and TRP. Also LEU and VAL must be \
defined.\n\n");
}


/************************************************************************/
/*>void DoFlipping(PDB *pdb, BOOL verbose, BOOL quiet)
   ---------------------------------------------------
*//**

   Does the main work of isolating residues and seeing if a flip needs
   to be made.

   Doesn't handle Leu or Val as these are not planar. However, the
   naming convention needs to be checked (not yet done).

-  08.11.96 Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  19.08.14 Removed unused variable. By: CTP
*/
void DoFlipping(PDB *pdb, BOOL verbose, BOOL quiet)
{
   PDB  *p, *q,
        *end,
        *atom1, *atom2, *atom3, *atom4, *atom4b,
        *connect4, *connect4b;
   int  res;
   REAL tor1, tor2;
   
   static TORSION torsions[] = 
   {
   {  "ARG ", "CD  ", "NE  ", "CZ  ", "NH1 ", "NH2 ", NULL  , NULL  },
   {  "ASP ", "CA  ", "CB  ", "CG  ", "OD1 ", "OD2 ", NULL  , NULL  },
   {  "GLU ", "CB  ", "CG  ", "CD  ", "OE1 ", "OE2 ", NULL  , NULL  },
/* {  "LEU ", "CA  ", "CB  ", "CG  ", "CD1 ", "CD2 ", NULL  , NULL  }, */
   {  "PHE ", "CA  ", "CB  ", "CG  ", "CD1 ", "CD2 ", "CE1 ", "CE2 "},
   {  "TYR ", "CA  ", "CB  ", "CG  ", "CD1 ", "CD2 ", "CE1 ", "CE2 "},
/* {  "VAL ", "N   ", "CA  ", "CB  ", "CG1 ", "CG2 ", NULL  , NULL  }, */
   {  NULL,   NULL,   NULL,   NULL,   NULL,   NULL,   NULL,   NULL  }
   };
   
   for(p=pdb; p!=NULL; p=end)
   {
      end = blFindNextResidue(p);

      /* Run through the list of residues which need treating           */
      for(res=0; torsions[res].resnam!=NULL; res++)
      {
         atom1 = atom2 = atom3 = atom4 = atom4b = NULL;
         connect4 = connect4b = NULL;
         
         /* If this one is in the list                                  */
         if(!strncmp(torsions[res].resnam, p->resnam, 4))
         {
            /* Find all the atoms of interest                           */
            for(q=p; q!=end; NEXT(q))
            {
               if(!strncmp(q->atnam,torsions[res].atom1,4))
               {
                  atom1 = q;
               }
               else if(!strncmp(q->atnam,torsions[res].atom2,4))
               {
                  atom2 = q;
               }
               else if(!strncmp(q->atnam,torsions[res].atom3,4))
               {
                  atom3 = q;
               }
               else if(!strncmp(q->atnam,torsions[res].atom4,4))
               {
                  atom4 = q;
               }
               else if(!strncmp(q->atnam,torsions[res].atom4b,4))
               {
                  atom4b = q;
               }
               else if(torsions[res].connect4 != NULL)
               {
                  if(!strncmp(q->atnam,torsions[res].connect4,4))
                     connect4 = q;
                  else if(!strncmp(q->atnam,torsions[res].connect4b,4))
                     connect4b = q;
               }
            }
            
            /* Check we found all atoms                                 */
            if((atom1==NULL) || (atom2==NULL) ||
               (atom3==NULL) || (atom4==NULL) ||
               (atom4b==NULL))
            {
               if(!quiet)
                  fprintf(stderr,"Warning: Missing atoms in %s %c%d%c, \
not processed.\n",p->resnam,p->chain[0],p->resnum,p->chain[0]);
               continue;  /* With next residue                          */
            }
            
            /* Calculate the two torsions                               */
            tor1 = blPhi(atom1->x,  atom1->y,  atom1->z,
                         atom2->x,  atom2->y,  atom2->z,
                         atom3->x,  atom3->y,  atom3->z,
                         atom4->x,  atom4->y,  atom4->z);
            tor2 = blPhi(atom1->x,  atom1->y,  atom1->z,
                         atom2->x,  atom2->y,  atom2->z,
                         atom3->x,  atom3->y,  atom3->z,
                         atom4b->x, atom4b->y, atom4b->z);

            /* See if we need to do a flip                              */
            if(ABS(tor2) < ABS(tor1))
            {
               if(verbose)
                  fprintf(stderr,"Flipping %s %c%d%c\n",
                          p->resnam, p->chain[0],p->resnum,p->insert[0]);
               DoAFlip(atom4,atom4b,connect4,connect4b);
            }
         }
      }
   }
}


/************************************************************************/
/*>void DoAFlip(PDB *atom4, PDB *atom4b, PDB *connect4, PDB *connect4b)
   --------------------------------------------------------------------
*//**

   Does the actual flip of atom names

-  08.11.96 Original   By: ACRM
*/
void DoAFlip(PDB *atom4, PDB *atom4b, PDB *connect4, PDB *connect4b)
{
   char temp[8];
   
   strcpy(temp, atom4->atnam);
   strcpy(atom4->atnam, atom4b->atnam);
   strcpy(atom4b->atnam, temp);
   
   if((connect4!=NULL) && (connect4b!=NULL))
   {
      strcpy(temp, connect4->atnam);
      strcpy(connect4->atnam, connect4b->atnam);
      strcpy(connect4b->atnam, temp);
   }
}

