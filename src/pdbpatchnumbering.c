/************************************************************************/
/**

   \file       pdbpatchnumbering.c
   
   \version    V1.9
   \date       30.09.17
   \brief      Patch the numbering of a PDB file from a file of numbers
               and sequence (as created by KabatSeq, etc)
   
   \copyright  (c) Dr. Andrew C. R. Martin / UCL 1995-2017
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
-  V1.0  09.08.95 Original    By: ACRM
-  V1.1  12.02.97 Fixed NULL pointer reference which was shown up under
                  Linux
-  V1.2  14.06.06 Now skips records in the patch file where the amino acid
                  is in lower case
-  V1.3  28.08.13 Modified for new ParseResSpec()
-  V1.4  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.5  06.11.14 Renamed from patchpdbnum By: ACRM
-  V1.6  07.11.14 Initialized a variable
-  V1.7  13.02.15 Added whole PDB support
-  V1.8  12.03.15 Changed to allow multi-character chain names and
                  three-letter code in patch file
-  V1.9  30.09.17 Now allows the patch file to skip the first few residues.
                  i.e. if there is an N-terminal extension to the known
                  numbering, this will be removed from the file.

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
   char          chain[8],
                 insert[8],
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
BOOL ApplyPatches(PDB **pPDB, PATCH *patches);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for PDB number patching

-  09.08.95 Original    By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   PATCH    *patches = NULL;
   FILE     *in      = stdin,
            *out     = stdout,
            *patchfp = NULL;
   char     infile[MAXBUFF],
            outfile[MAXBUFF],
            patchfile[MAXBUFF];
   WHOLEPDB *wpdb;
   PDB      *pdb;
   
   if(ParseCmdLine(argc, argv, infile, outfile, patchfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((patchfp=fopen(patchfile,"r"))==NULL)
         {
            fprintf(stderr,"pdbpatchnumbering: Unable to open patch \
file\n");
            return(1);
         }

         if((patches = ReadPatchFile(patchfp))==NULL)
         {
            fprintf(stderr,"pdbpatchnumbering: Unable to read patch \
file\n");
            return(1);
         }

         if((wpdb = blReadWholePDB(in)) != NULL)
         {
            pdb=wpdb->pdb;
            if(ApplyPatches(&pdb, patches))
            {
               wpdb->pdb = pdb;
               blWriteWholePDB(out, wpdb);
            }
            else
            {
               fprintf(stderr,"pdbpatchnumbering: Patching failed\n");
               return(1);
            }

            FREELIST(pdb, PDB);
         }
         else
         {
            fprintf(stderr,"pdbpatchnumbering: No atoms read from PDB \
file\n");
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
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *patchfile   Output file (or blank string)
   \return                     Success?

   Parse the command line
   
-  09.08.95 Original    By: ACRM
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
*//**

   Reads a patch file and creates a linked list of patches which is 
   returned.

-  09.08.95 Original    By: ACRM
-  14.06.06 Skips records where the amino acid is in lower case
            Also fixed for newer GetWord() which needs buffer size
-  28.08.13 Modified for new ParseResSpec()
-  07.11.14 Initialized p
-  12.03.15 Changed to allow multi-character chain names and
            to allow three-letter code
*/
PATCH *ReadPatchFile(FILE *fp)
{
   PATCH *patch = NULL,
         *p     = NULL;
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
            chp = blGetWord(buffer,resid,MAXBUFF);
            blGetWord(chp,aacode,MAXBUFF);

            if(isupper(aacode[0]))   /* 14.06.06                        */
            {
               /* If the first word is a valid residue specification    */
               if(blParseResSpec(resid, chain, &resnum, insert))
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
                  strncpy(p->chain, chain, 8);
                  p->resnum = resnum;
                  strncpy(p->insert, insert, 8);
                  if(strlen(aacode) > 1)
                  {
                     p->aacode = blThrone(aacode);
                  }
                  else
                  {
                     p->aacode = aacode[0];
                  }
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
*//**

   Prints a usage message

-  09.08.95 Original    By: ACRM
-  14.06.06 V1.2
-  28.08.13 V1.3
-  22.07.14 V1.4 By: CTP
-  07.11.14 V1.6 By: ACRM
-  13.02.15 V1.7 By: ACRM
-  12.03.15 V1.8
-  30.09.17 V1.9
*/
void Usage(void)
{
   fprintf(stderr,"\npdbpatchnumbering V1.9 (c) 1995-2017, Dr. Andrew \
C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: pdbpatchnumbering patchfile [in.pdb \
[out.pdb]]\n");
   fprintf(stderr,"PDB file I/O is through stdin/stdout if files are not \
specified.\n");

   fprintf(stderr,"\npdbpatchnumbering patches the numbering of a PDB \
file from a patch file\n");
   fprintf(stderr,"containing resspec residue specifiers:\n");

   blPrintResSpecHelp(stderr);

   fprintf(stderr,"\nThe numbering of the PDB file is modified to match \
that in the patch\n");
   fprintf(stderr,"file and the PDB file is terminated after all \
specified residues.\n\n");
   fprintf(stderr,"The patch file must contain all the residues present \
in the PDB file and\n");
   fprintf(stderr,"typically comes from a program such as abnum or \
abynum which applies\n");
   fprintf(stderr,"standard numbering to a sequence.\n\n");
   fprintf(stderr,"The patch file consists of records of the form:\n");
   fprintf(stderr,"L1 ALA          L1 A\n");
   fprintf(stderr,"L2 CYS   -or-   L2 C\n");
   fprintf(stderr,"L3 ASP          L3 D\n\n");
}


/************************************************************************/
/*>BOOL ApplyPatches(PDB **pPDB, PATCH *patches)
   --------------------------------------------
*//**

   Does the real work of applying the sequence numbering patches and
   unlinking unused parts of the PDB linked list.

   Walks the patch linked list to find required residue numbers and
   walks the PDB linked list applying them. When an end of chain is
   found in the patch list, the current chain in the PDB list is 
   truncated. If a PDB chain ends before a patch list chain, then an
   error is issued.

-  09.08.95 Original    By: ACRM
-  12.02.97 Fixed NULL pointer reference on first entry round loop
            (prevchain was set from r which was not initialised)
-  12.03.15 Changed to allow multi-character chain names
*/
BOOL ApplyPatches(PDB **pPDB, PATCH *patches)
{
   PDB   *p,                      /* Start of a residue                 */
         *q,                      /* Start of next residue              */
         *r = NULL,               /* For stepping through a residue     */
         *pp,                     /* Previous record to r               */
         *prev = NULL,            /* Previous record to p               */
         *pdb  = NULL;            /* Derefernced from pPDB              */
   PATCH *a;                      /* Step through patches               */
   char  patchchain[8],           /* Current patch chain                */
         pdbchain[8],             /* Current pdb chain                  */
         prevchain[8];            /* PDB chain before renaming          */
   BOOL  NewPatchChain = FALSE,   /* Indicates new chain in patch list  */
         NewPDBChain   = TRUE;    /* Indicates new chain in PDB file    */
   int   resnum,                  /* Current patch residue number       */
         skippedResidues = 0,     /* Count Nter skipped residues        */
         maxSkippedResidues = 5;  /* Max Nter residues to skip          */
   
   pdb = *pPDB;

   if((patches != NULL) && (pdb != NULL))
   {
      /* Record the chains we are starting in                           */
      strncpy(patchchain, patches->chain, 8);
      strncpy(pdbchain,   pdb->chain,     8);
      resnum     = patches->resnum;
      strncpy(prevchain, pdbchain, 8); /* 12.02.97 Initialised this - 
                                          shouldn't be needed, but just
                                          in case...     
                                       */

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
            if(!CHAINMATCH(a->chain, patchchain) || (a->resnum < resnum))
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
                  for(r=p; CHAINMATCH(r->chain, prevchain); NEXT(r))
                     pp = r;

                  /* If we skipped over anything, unlink and free       */
                  if(pp!=NULL)
                  {
                     prev->next=r;
                     pp->next = NULL;
                     FREELIST(p,PDB);     
                     p=r;
                  }
               }
               
               skippedResidues = 0;
               NewPDBChain     = TRUE;
               NewPatchChain   = TRUE;
               strncpy(patchchain, a->chain, 8);
               strncpy(pdbchain,   p->chain, 8);
            }

            /* If the patch list hasn't changed chain, see if the PDB
               list has. This is an error as the PDB list is truncated
            */
            if(!NewPatchChain)
            {
               if(!CHAINMATCH(p->chain, pdbchain))
               {
                  fprintf(stderr,"pdbpatchnumbering: Chain %s too short \
for patches\n",pdbchain);
                  return(FALSE);
               }
            }

            /* 30.09.17 Skip up to the first maxSkippedResidues if the 
               amino acids don't match
            */
            while(NewPDBChain && blThrone(p->resnam) != a->aacode)
            {
               if(skippedResidues++ < maxSkippedResidues)
               {
                  p = blDeleteResiduePDB(pPDB, p);
                  pdb = *pPDB;
               }
            }

            NewPDBChain = FALSE;

            if(blThrone(p->resnam) != a->aacode)
            {
               fprintf(stderr,"Residue mismatch between patch file and \
PDB file.\n");
               fprintf(stderr,"Patch file expects amino acid %c. PDB \
record is:\n",a->aacode);
               blWritePDBRecord(stderr,p);
               return(FALSE);
            }
               
            /* Find the next residue in the PDB file                    */
            q = blFindNextResidue(p);

            /* Apply the new residue numbering to the current residue   */
            if(r!=NULL)                  /* 12.02.97 Added NULL check   */
               strncpy(prevchain, r->chain, 8);
            for(r=p; r!=q; NEXT(r))
            {
               strcpy(r->chain,  a->chain);
               strcpy(r->insert, a->insert);
               r->resnum = a->resnum;
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


