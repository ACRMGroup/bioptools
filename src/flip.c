/*************************************************************************

   Program:    flip
   File:       flip.c
   
   Version:    V1.0
   Date:       08.11.96
   Function:   Standardise equivalent atom labelling
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1996
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0)1372 275775
               (Work) +44 (0)171 419 3890
   EMail:      INTERNET: martin@biochem.ucl.ac.uk
               
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
   Main program to flip sidechain equivalent atom naming

   08.11.96 Original   By: ACRM
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
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = ReadPDB(in,&natoms)) != NULL)
         {
            DoFlipping(pdb,verbose,quiet);
            WritePDB(out, pdb);
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
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            BOOL   *verbose     Report each flip made?
            BOOL   *quiet       Don't report skipped residues
   Returns: BOOL                Success?

   Parse the command line
   
   08.11.96 Original    By: ACRM
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
   Prints a usage message

   08.11.96 Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nflip V1.0 (c) Dr. Andrew C.R. Martin, UCL\n");
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
   Does the main work of isolating residues and seeing if a flip needs
   to be made.

   Doesn't handle Leu or Val as these are not planar. However, the
   naming convention needs to be checked (not yet done).

   08.11.96 Original   By: ACRM
*/
void DoFlipping(PDB *pdb, BOOL verbose, BOOL quiet)
{
   PDB  *p, *q,
        *end,
        *atom1, *atom2, *atom3, *atom4, *atom4b,
        *connect4, *connect4b;
   int  i, res;
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
      end = FindNextResidue(p);

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
            tor1 = phi(atom1->x,  atom1->y,  atom1->z,
                       atom2->x,  atom2->y,  atom2->z,
                       atom3->x,  atom3->y,  atom3->z,
                       atom4->x,  atom4->y,  atom4->z);
            tor2 = phi(atom1->x,  atom1->y,  atom1->z,
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
   Does the actual flip of atom names

   08.11.96 Original   By: ACRM
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

