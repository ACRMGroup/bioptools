/************************************************************************/
/**

   \file       pdbpatchbval.c
   
   \version    V1.7
   \date       12.03.15
   \brief      Patch the b-value (or occupancy) column using values from 
               a file
   
   \copyright  (c) Dr. Andrew C. R. Martin / UCL 1996-2015
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
-  V1.0  29.05.96 Original   By: ACRM
-  V1.1  17.02.97 Fixed bug in not updating doubly linked list correctly
                  when removing items from patch list
-  V1.2  28.08.13 Modified for new ParseResSpec()
-  V1.3  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.4  06.11.14 Renamed from patchbval  By: ACRM
-  V1.5  25.11.14 Initialized a variable 
-  V1.6  13.02.15 Added whole PDB support
-  V1.7  12.03.15 Changed to allow multi-character chain names

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

typedef struct _patch
{
   struct _patch *next, 
                 *prev;
   REAL value;
   int  resnum;
   char chain[8],   /* 28.08.13 Now a string                            */
        insert[8];  /* 28.08.13 Now a string                            */
}  PATCH;


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ApplyPatches(FILE *in, FILE *out, PATCH *patchlist, BOOL occup,
                  BOOL verbose);
PATCH *ReadPatchFile(FILE *fp);
BOOL ParseCmdLine(int argc, char **argv, char *datafile, char *infile, 
                  char *outfile, BOOL *occup, BOOL *verbose);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for patching bvalues or occupancies per-residue

-  29.05.96 Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE  *in        = stdin,
         *out       = stdout,
         *data      = NULL;
   PATCH *patchlist = NULL;
   char  infile[MAXBUFF],
         outfile[MAXBUFF],
         datafile[MAXBUFF];
   BOOL  occup   = FALSE,
         verbose = FALSE;

   if(ParseCmdLine(argc, argv, datafile, infile, outfile, &occup, 
                   &verbose))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((data=fopen(datafile,"r"))!=NULL)
         {
            if((patchlist = ReadPatchFile(data))!=NULL)
            {
               if(!ApplyPatches(in, out, patchlist, occup, verbose))
               {
                  fprintf(stderr,"pdbpatchbval: Patching failed!\n");
                  return(1);
               }
               /* Normal exit from here!                                */
            }
            else
            {
               fprintf(stderr,"pdbpatchbval: Unable to read patch \
data\n");
               return(1);
            }
         }
         else
         {
            fprintf(stderr,"pdbpatchbval: Unable to open patch file: \
%s\n", datafile);
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
/*>BOOL ApplyPatches(FILE *in, FILE *out, PATCH *patch, BOOL occup,
                     BOOL verbose)
   --------------------------------------------------------------------
*//**

   Does the work of patching a PDB file, including reading it in and
   writing out the patched version

-  29.05.96 Original   By: ACRM
-  17.02.97 The linked list wasn't being updated in the prev direction
            when items were removed leading to possible corruption.
            Also, no longer treats end of list as a separate special
            case.
-  28.08.13 PATCH.chain and PATCH.insert are now strings
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  12.03.15 Changed to allow multi-character chain names  By: ACRM
*/
BOOL ApplyPatches(FILE *in, FILE *out, PATCH *patch, BOOL occup,
                  BOOL verbose)
{
   WHOLEPDB *wpdb;
   PDB      *pdb,
            *end,
            *p,
            *q;
   PATCH    *pa;
   
   
   if((wpdb = blReadWholePDB(in))==NULL)
   {
      fprintf(stderr,"Unable to read PDB file\n");
      return(FALSE);
   }

   pdb = wpdb->pdb;

   /* Set values to 0.0                                                 */
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(occup) p->occ  = (REAL)0.0;
      else      p->bval = (REAL)0.0;
   }

   /* Run through again a residue at a time                             */
   for(p=pdb; p!=NULL; p=end)
   {
      /* If all patches applied, break out of the linked list           */
      if(patch==NULL)
         break;
      
      end = blFindNextResidue(p);
      
      /* Now compare this residue with items in the patch list          */
      for(pa=patch; pa!=NULL; NEXT(pa))
      {
         if((pa->resnum == p->resnum)          &&
            CHAINMATCH(pa->chain, p->chain) &&
            (pa->insert[0] == p->insert[0]))
         {
            break;
         }
      }
      

      if(pa!=NULL)
      {
         /* Since we broke out, then we have a match so apply patch     */
         for(q=p; q!=end; NEXT(q))
         {
            if(occup) q->occ  = pa->value;
            else      q->bval = pa->value;
         }

         /* Unlink this from the patch list                             */
         if(pa==patch)                    /* First in list              */
         {
            patch = pa->next;
            if(patch!=NULL)
               patch->prev = NULL;
         }
         else
         {
            if(pa->prev != NULL)
               pa->prev->next = pa->next;
            /* 17.02.97 Added update this way round                     */
            if(pa->next != NULL)
               pa->next->prev = pa->prev;
         }
         free(pa);
      }
   }

   /* All patches should now ba applied; list any which weren't         */
   if(patch!=NULL && verbose)
   {
      fprintf(stderr,"The following patches were not applied:\n");
      for(pa=patch; pa!=NULL; NEXT(pa))
      {
         fprintf(stderr,"%s.%d%c %f\n",
                 pa->chain,pa->resnum,pa->insert[0],pa->value);
      }
   }
   
   /* Write the patched PDB file                                        */
   blWriteWholePDB(out,wpdb);

   return(TRUE);
}


/************************************************************************/
/*>PATCH *ReadPatchFile(FILE *fp)
   ------------------------------
*//**

   Reads a patch file into a linked list

-  29.05.96 Original   By: ACRM
-  28.08.13 Modified for new ParseResSpec()
-  25.11.14 Initialized p
*/
PATCH *ReadPatchFile(FILE *fp)
{
   char  buffer[MAXBUFF],
         resspec[MAXBUFF];
   REAL  value;
   PATCH *patch = NULL,
         *p     = NULL;
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      TERMINATE(buffer);
      if((sscanf(buffer,"%s %lf",resspec, &value)) == 2)
      {
         if(patch == NULL)
         {
            INITPREV(patch, PATCH);
            p = patch;
         }
         else
         {
            ALLOCNEXTPREV(p, PATCH);
         }
         if(p==NULL)
         {
            FREELIST(patch, PATCH);
            return(NULL);
         }
         blParseResSpec(resspec, p->chain, &p->resnum, p->insert);
         p->value = value;
      }
   }

   return(patch);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *datafile, char *infile, 
                     char *outfile, BOOL *occup, BOOL *verbose)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *datafile    The patch datafile
   \param[out]     *occup       Put the patches in the occupancy column
   \param[out]     *verbose     Report failed patches
   \return                     Success?

   Parse the command line
   
-  29.05.96 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *datafile, char *infile, 
                  char *outfile, BOOL *occup, BOOL *verbose)
{
   argc--;
   argv++;

   *occup = FALSE;
   *verbose = FALSE;
   
   infile[0] = outfile[0] = datafile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'v':
            *verbose = TRUE;
            break;
         case 'o':
            *occup = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are up to 3 arguments left                 */
         if(argc > 3)
            return(FALSE);
         
         /* Copy the first to datafile                                  */
         strcpy(datafile, argv[0]);
         
         /* If there's another, copy it to infile                       */
         argc--;
         argv++;
         if(argc)
         {
            strcpy(infile, argv[0]);
            
            /* If there's another, copy it to outfile                   */
            argc--;
            argv++;
            if(argc)
               strcpy(outfile, argv[0]);
         }
         
         return(TRUE);
      }
      argc--;
      argv++;
   }

   if(!datafile[0])
      return(FALSE);
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  29.05.96 Original   By: ACRM
-  28.08.13 V1.2
-  22.07.14 V1.3 By: CTP
-  06.11.14 V1.4 By: ACRM
-  25.11.14 V1.5
-  13.02.15 V1.6
-  13.03.15 V1.7
*/
void Usage(void)
{
   fprintf(stderr,"\npdbpatchbval V1.7 (c) 1996-2015, Dr. Andrew C.R. \
Martin, UCL\n");

   fprintf(stderr,"\nUsage: pdbpatchbval [-o] [-v] patchfile [in.pdb \
[out.pdb]]\n");
   fprintf(stderr,"       -o  Place the patches in the occupancy \
column\n");
   fprintf(stderr,"       -v  Verbose: report failed patches\n");

   fprintf(stderr,"\npdbpatchbval takes a patch file containing \
resspec residue specifications\n");
   fprintf(stderr,"and values one to a line and patches the B-value (or \
occupancy) for\n");
   fprintf(stderr,"that residue with the specified values.\n");
   blPrintResSpecHelp(stderr);
   fprintf(stderr,"\n");
}



 
