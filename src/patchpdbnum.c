/*************************************************************************

   Program:    PatchPDBNum
   File:       patchpdbnum.c
   
   Version:    V1.3
   Date:       28.08.13
   Function:   Patch the numbering of a PDB file from a file of numbers
               and sequence (as created by KabatSeq, etc)
   
   Copyright:  (c) Dr. Andrew C. R. Martin / UCL 1995-2013
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
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
   V1.0  09.08.95 Original    By: ACRM
   V1.1  12.02.97 Fixed NULL pointer reference which was shown up under
                  Linux
   V1.2  14.06.06 Now skips records in the patch file where the amino acid
                  is in lower case
   V1.3  28.08.13 Modified for new ParseResSpec()

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"


/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

typedef struct _patch
{
   struct _patch *next;
   int           resnum;
   char          chain,
                 insert,
                 aacode;
}  PATCH;


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *patchfile);
PATCH *ReadPatchFile(FILE *fp);
void Usage(void);
BOOL ApplyPatches(PDB *pdb, PATCH *patches);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for PDB number patching

   09.08.95 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   PATCH *patches = NULL;
   FILE  *in      = stdin,
         *out     = stdout,
         *patchfp = NULL;
   char  infile[MAXBUFF],
         outfile[MAXBUFF],
         patchfile[MAXBUFF];
   PDB   *pdb;
   int   natoms;
   
   if(ParseCmdLine(argc, argv, infile, outfile, patchfile))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((patchfp=fopen(patchfile,"r"))==NULL)
         {
            fprintf(stderr,"patchpdbnum: Unable to open patch file\n");
            return(1);
         }

         if((patches = ReadPatchFile(patchfp))==NULL)
         {
            fprintf(stderr,"patchpdbnum: Unable to read patch file\n");
            return(1);
         }

         if((pdb = ReadPDB(in,&natoms)) != NULL)
         {
            if(ApplyPatches(pdb, patches))
            {
               WritePDB(out, pdb);
            }
            else
            {
               fprintf(stderr,"patchpdbnum: Patching failed\n");
               return(1);
            }

            FREELIST(pdb, PDB);
         }
         else
         {
            fprintf(stderr,"patchpdbnum: No atoms read from PDB file\n");
         }

         FREELIST(patches, PATCH);
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
                     char *patchfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            char   *patchfile   Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   09.08.95 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *patchfile)
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
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 1--3 arguments left                    */
         if(argc < 1 || argc > 3)
            return(FALSE);
         
         /* Copy the first to patchfile                                 */
         strcpy(patchfile, argv[0]);
         argc--;
         argv++;
         
         if(argc)
         {
            /* Copy the next to infile                                  */
            strcpy(infile, argv[0]);
            argc--;
            argv++;
         
            /* If there's another, copy it to outfile                   */
            if(argc)
               strcpy(outfile, argv[0]);
         }
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>PATCH *ReadPatchFile(FILE *fp)
   ------------------------------
   Reads a patch file and creates a linked list of patches which is 
   returned.

   09.08.95 Original    By: ACRM
   14.06.06 Skips records where the amino acid is in lower case
            Also fixed for newer GetWord() which needs buffer size
   28.08.13 Modified for new ParseResSpec()
*/
PATCH *ReadPatchFile(FILE *fp)
{
   PATCH *patch = NULL,
         *p;
   char  buffer[MAXBUFF],
         resid[MAXBUFF],
         aacode[MAXBUFF],
         chain[8],
         insert[8],
         *chp;
   int   resnum;

   /* Get lines from the file                                           */
   while(fgets(buffer,MAXBUFF,fp))
   {
      TERMINATE(buffer);
      
      /* If it's not a blank line or a comment                          */
      if(strlen(buffer) && buffer[0] != '!' && buffer[0] != '#')
      {
         /* If it's a warning or error, echo to stderr                  */
         if(!strncmp(buffer,"WARNING:",8) ||
            !strncmp(buffer,"ERROR:",6))
         {
            fprintf(stderr,"Patch file %s\n",buffer);
         }
         else
         {
            /* Get the first two words out of the buffer                */
            chp = GetWord(buffer,resid,MAXBUFF);
            GetWord(chp,aacode,MAXBUFF);

            if(isupper(aacode[0]))   /* 14.06.06                        */
            {
               /* If the first word is a valid residue specification    */
               if(ParseResSpec(resid, chain, &resnum, insert))
               {
                  /* Allocate an item in the linked list                */
                  if(patch==NULL)
                  {
                     INIT(patch,PATCH);
                     p = patch;
                  }
                  else
                  {
                     ALLOCNEXT(p,PATCH);
                  }
                  if(p==NULL)
                  {
                     if(patch!=NULL) FREELIST(patch, PATCH);
                     return(NULL);
                  }
                  
                  /* Copy the data into the linked list entry           */
                  p->chain  = chain[0];
                  p->resnum = resnum;
                  p->insert = insert[0];
                  p->aacode = aacode[0];
               }
            }
         }
      }
   }

   return(patch);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
   Prints a usage message

   09.08.95 Original    By: ACRM
   14.06.06 V1.2
   28.08.13 V1.3
*/
void Usage(void)
{
   fprintf(stderr,"\nPatchPDBNum V1.3 (c) 1995-2013, Dr. Andrew C.R. \
Martin, UCL\n");

   fprintf(stderr,"\nUsage: patchpdbnum <patchfile> [<in.pdb> \
[<out.pdb>]]\n");
   fprintf(stderr,"PDB file I/O is through stdin/stdout if files are not \
specified.\n");

   fprintf(stderr,"\nPatchPDBNum patches the numbering of a PDB file \
from a patch file\n");
   fprintf(stderr,"containing residue numbers in the form [c]num[i] \
(where c is an optional\n");
   fprintf(stderr,"chain name, num is the residue number and i is an \
optional insert code)\n");
   fprintf(stderr,"followed by a 1-letter amino acid name. The numbering \
of the PDB file\n");
   fprintf(stderr,"is modified to match that in the patch file and the \
PDB file is \n");
   fprintf(stderr,"terminated after all specified residues.\n\n");
}


/************************************************************************/
/*>BOOL ApplyPatches(PDB *pdb, PATCH *patches)
   -------------------------------------------
   Does the real work of applying the sequence numbering patches and
   unlinking unused parts of the PDB linked list.

   Walks the patch linked list to find required residue numbers and
   walks the PDB linked list applying them. When an end of chain is
   found in the patch list, the current chain in the PDB list is 
   truncated. If a PDB chain ends before a patch list chain, then an
   error is issued.

   09.08.95 Original    By: ACRM
   12.02.97 Fixed NULL pointer reference on first entry round loop
            (prevchain was set from r which was not initialised)
*/
BOOL ApplyPatches(PDB *pdb, PATCH *patches)
{
   PDB   *p,                      /* Start of a residue                 */
         *q,                      /* Start of next residue              */
         *r = NULL,               /* For stepping through a residue     */
         *pp,                     /* Previous record to r               */
         *prev = NULL;            /* Previous record to p               */
   PATCH *a;                      /* Step through patches               */
   char  patchchain,              /* Current patch chain                */
         pdbchain,                /* Current pdb chain                  */
         prevchain;               /* PDB chain before renaming          */
   BOOL  NewPatchChain = FALSE;   /* Indicates new chain in patch list  */
   int   resnum;                  /* Current patch residue number       */

   if((patches != NULL) && (pdb != NULL))
   {
      /* Record the chains we are starting in                           */
      patchchain = patches->chain;
      resnum     = patches->resnum;
      pdbchain   = pdb->chain[0];
      prevchain  = pdbchain;      /* 12.02.97 Initialised this - shouldn't
                                     be needed, but just in case...     */

      /* Initialise p to start of PDB linked list                       */
      p = pdb;
      
      /* Step through the patch linked list                             */
      for(a=patches; a!=NULL; NEXT(a))
      {
         /* If it's not an insertion                                    */
         if(a->aacode != '-')
         {
            /* See if the chain has just changed in the patch file      */
            NewPatchChain = FALSE;
            if((a->chain != patchchain) ||(a->resnum < resnum))
            {
               /* Chain has ended in patch file, so unlink up to start of
                  next chain from PDB file
               */
               if(prev != NULL)
               {
                  /* Walk along from prev till we find the next chain,
                     keeping pp as the previous item in the bit we are 
                     junking
                  */
                  pp=NULL;
                  for(r=p; r->chain[0]==prevchain; NEXT(r)) pp = r;

                  /* If we skipped over anything, unlink and free       */
                  if(pp!=NULL)
                  {
                     prev->next=r;
                     pp->next = NULL;
                     FREELIST(p,PDB);     
                     p=r;
                  }
               }
               
               NewPatchChain = TRUE;
               patchchain    = a->chain;
               pdbchain      = p->chain[0];
            }

            /* If the patch list hasn't changed chain, see if the PDB
               list has. This is an error as the PDB list is truncated
            */
            if(!NewPatchChain)
            {
               if(p->chain[0] != pdbchain)
               {
                  fprintf(stderr,"patchpdbnum: Chain %c too short for \
patches\n",pdbchain);
                  return(FALSE);
               }
            }

            /* Check the AA code is correct                             */
            if(throne(p->resnam) != a->aacode)
            {
               fprintf(stderr,"Residue mismatch between patch file and \
PDB file.\n");
               fprintf(stderr,"Patch file expects residue %c. PDB record \
is:\n",a->aacode);
               WritePDBRecord(stderr,p);
               return(FALSE);
            }
               
            /* Find the next residue in the PDB file                    */
            q = FindNextResidue(p);

            /* Apply the new residue numbering to the current residue   */
            if(r!=NULL)                  /* 12.02.97 Added NULL check   */
               prevchain = r->chain[0];
            for(r=p; r!=q; NEXT(r))
            {
               r->chain[0]  = a->chain;
               r->resnum    = a->resnum;
               r->insert[0] = a->insert;
            }

            /* Update record of previous item in linked list so we can
               unlink chunks
            */
            for(prev=p; prev->next!=q; NEXT(prev));
            
            /* Update start of residue PDB pointer                      */
            p = q;
         }  /* It's not an insert in the patch file                     */
         /* Update current residue number                               */
         resnum = a->resnum;
      }  /* End of stepping though patches                              */
   }  /* End of check on patches and PDB                                */
   else
   {
      return(FALSE);
   }

   /* Unlink anything which remains in the PDB linked list              */
   FREELIST(prev->next, PDB);
   prev->next = NULL;
      
   return(TRUE);
}


