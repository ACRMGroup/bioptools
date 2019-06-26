/************************************************************************/
/**

   \file       pdbgetresidues.c
   
   \version    V1.5.1
   \date       26.06.19
   \brief      Extract a set of residues from a PDB file
   
   \copyright  (c) Prof. Andrew C. R. Martin 2010-2019
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
-  V1.0   15.06.10 Original
-  V1.1   22.07.14 Renamed deprecated functions with bl prefix.
                   Added doxygen annotation. By: CTP
-  V1.2   06.11.14 Renamed from getresidues  By: ACRM
-  V1.3   25.11.14 Initialized a variable
-  V1.4   12.03.15 Changed to allow multi-character chain names
-  V1.5   13.03.19 Terminate strings after strncpy()
-  V1.5.1 26.06.19 Fix to termination of resnam which was beyond array
                   boundary

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
        chain[8],
        insert[8];
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
*//**

   main routine

-  22.07.96 Original    By: ACRM
-  29.09.05 Modified for -l By: TL
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
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
         if(blOpenStdFiles(InFile, OutFile, &in, &out))
         {
            if((pdb=blReadPDB(in, &natom))==NULL)
            {
               fprintf(stderr,"Error: pdbgetresidues - No atoms read from \
PDB file\n");
               return(1);
            }
            
            if((reslist = ReadResidueList(rfp))==NULL)
            {
               fprintf(stderr,"Error: pdbgetresidues - Failed to read \
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
/*
-  12.03.15 Changed to allow multi-character chain names
*/
void PrintResidues(FILE *out, PDB *pdb, RESLIST *reslist)
{
   PDB *p;
   RESLIST *r;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      for(r=reslist; r!=NULL; NEXT(r))
      {
         if((p->resnum == r->resnum) &&
            CHAINMATCH(p->chain, r->chain) &&
            !strcmp(p->insert, r->insert))
         {
            blWritePDBRecord(out, p);
         }
      }
   }
}


/************************************************************************/
/*
-  13.03.19  Terminate string ofter strncpy()
-  26.06.19  Array boundary fix
*/
RESLIST *ReadResidueList(FILE *fp)
{
   char buffer[MAXBUFF];
   char chain[8], insert[8];
   int  resnum;
   RESLIST *reslist = NULL, 
           *r       = NULL;  /* ACRM 25.11.15 Initialize this           */
   
   while(fgets(buffer, MAXBUFF, fp))
   {
      TERMINATE(buffer);
      if(blParseResSpec(buffer, chain, &resnum, insert))
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
         r->resid[MAXRESID-1] = '\0'; /* ACRM 26.06.19 Added the -1     */

         strncpy(r->chain, chain, 8);
         r->chain[7] = '\0';

         strncpy(r->insert, insert, 8);
         r->insert[7] = '\0';

         r->resnum = resnum;
      }
   }
   
   return(reslist);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *resfile,
                     char *infile, char *outfile)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *Zone1       First end of zone
   \param[out]     *Zone2       Second end of zone
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \return                      Success?

   Parse the command line
   
-  15.06.10 Original    By: ACRM
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
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  22.07.14 V1.1 By: CTP
-  06.11.14 V1.2 By: ACRM
-  25.11.14 V1.3 By: ACRM
-  12.03.15 V1.4
-  13.03.19 V1.5
*/
void Usage(void)
{
   fprintf(stderr,"\npdbgetresidues V1.5 (c) 2010-2019, UCL, Dr. Andrew \
C.R. Martin\n");
   fprintf(stderr,"\nUsage: pdbgetresidues resfile [in.pdb [out.pdb]]\n");

   fprintf(stderr,"\nresfile is a file listing residue specifications in \
the resspec format:\n\n");
   blPrintResSpecHelp(stderr);

   fprintf(stderr,"\nInput is from stdin and output is to stdout if no \
files are specified. \n");

   fprintf(stderr,"\nTakes a list of residue specifications and extracts \
just those residues \n");
   fprintf(stderr,"from a PDB file. Note that residue specifications are \
case sensitive.\n");
}
