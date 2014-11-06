/*************************************************************************

   Program:    getresidues
   File:       getresidues.c
   
   Version:    V1.0
   Date:       15.06.10
   Function:   Extract a set of residues from a PDB file
   
   Copyright:  (c) Dr. Andrew C. R. Martin 2010
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
   V1.0  15.06.10  Original
                    
*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define MAXRESID 16
typedef struct _reslist
{
   char resid[MAXRESID],
        chain,
        insert;
   int  resnum;
   struct _reslist *next;
}  RESLIST;

   

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *resfile,
                  char *infile, char *outfile);
void Usage(void);
RESLIST *ReadResidueList(FILE *fp);
void PrintResidues(FILE *out, PDB *pdb, RESLIST *reslist);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   main routine

   22.07.96 Original    By: ACRM
   29.09.05 Modified for -l By: TL
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout,
        *rfp = NULL;
   char InFile[MAXBUFF],
        OutFile[MAXBUFF],
        ResFile[MAXBUFF];
   int  natom;
   PDB  *pdb;
   RESLIST *reslist = NULL;


   if(ParseCmdLine(argc, argv, ResFile, InFile, OutFile))
   {
      if((rfp=fopen(ResFile, "r"))!=NULL)
      {
         if(OpenStdFiles(InFile, OutFile, &in, &out))
         {
            if((pdb=ReadPDB(in, &natom))==NULL)
            {
               fprintf(stderr,"Error: getresidues - No atoms read from \
PDB file\n");
               return(1);
            }
            
            if((reslist = ReadResidueList(rfp))==NULL)
            {
               fprintf(stderr,"Error: getresidues - Failed to read \
residues from list\n");
               return(1);
            }
            
            PrintResidues(out, pdb, reslist);
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
void PrintResidues(FILE *out, PDB *pdb, RESLIST *reslist)
{
   PDB *p;
   RESLIST *r;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(r=reslist; r!=NULL; NEXT(r))
      {
         if((p->resnum == r->resnum) &&
            (p->chain[0] == r->chain) &&
            (p->insert[0] == r->insert))
         {
            WritePDBRecord(out, p);
         }
      }
   }
}


/************************************************************************/
RESLIST *ReadResidueList(FILE *fp)
{
   char buffer[MAXBUFF];
   char chain[8], insert[8];
   int  resnum;
   RESLIST *reslist = NULL, *r;
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      TERMINATE(buffer);
      if(ParseResSpec(buffer, chain, &resnum, insert))
      {
         if(reslist == NULL)
         {
            INIT(reslist, RESLIST);
            r = reslist;
         }
         else
         {
            ALLOCNEXT(r, RESLIST);
         }
         if(r==NULL)
         {
            FREELIST(reslist, RESLIST);
            return(NULL);
         }
         
         strncpy(r->resid, buffer, MAXRESID);
         r->chain = chain[0];
         r->insert = insert[0];
         r->resnum = resnum;
      }
   }
   
   return(reslist);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *resfile,
                     char *infile, char *outfile)
   ----------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *Zone1       First end of zone
            char   *Zone2       Second end of zone
            char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
            BOOL   *uppercaseresspec  Should residue spec be upcased?
                                (Default: yes)
   Returns: BOOL                Success?

   Parse the command line
   
   15.06.10 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *resfile,
                  char *infile, char *outfile)
{
   argc--;
   argv++;

   resfile[0] = infile[0] = outfile[0] = '\0';

   if(!argc)
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
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 1, 2 or 3 arguments left               */
         if(argc < 1 || argc > 3)
            return(FALSE);
         
         /* Copy the first to resfile                                   */
         strcpy(resfile, argv[0]);
         argc--;
         argv++;
         
         /* Copy the next  to infile                                    */
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
void Usage(void)
{
   fprintf(stderr,"\ngetresidues V1.0 (c) 2010, UCL, Dr. Andrew C.R. \
Martin\n");
   fprintf(stderr,"\nUsage: getresidues resfile [in.pdb [out.pdb]]\n");

   fprintf(stderr,"\nresfile is a file listing residue specifications in \
the format \n");
   fprintf(stderr,"[c][.]nnn[i]\n");
   fprintf(stderr,"where [c] is an optional chain specification which \
may be followed by a \n");
   fprintf(stderr,"full stop (required if the chain label is a number), \
nnn is a residue \n");
   fprintf(stderr,"number and [i] is an optional insert code.\n");
   fprintf(stderr,"Input is from stdin and output is to stdout if no \
files are specified. \n");

   fprintf(stderr,"\nTakes a list of residue specifications and extracts \
just those residues \n");
   fprintf(stderr,"from a PDB file\n");
}
