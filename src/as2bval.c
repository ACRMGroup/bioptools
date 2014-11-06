/*************************************************************************

   Program:    as2bval
   File:       as2bval.c
   
   Version:    V1.5
   Date:       30.09.05
   Function:   Sum B-vals over each residue and replace with the
               summed or average value
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1994-2005
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
   V1.0  05.07.94 Original    By: ACRM
   V1.1  24.08.94 Changed to call OpenStdFiles()
   V1.2  27.09.94 Corrected usage message.
   V1.3  30.05.02 Changed PDB field from 'junk' to 'record_type'
   V1.4  03.06.04 Fixed for FixAtomName() which now needs occupancy
   V1.5  30.09.05 Now writes the results itself using WritePDBRecordAtnam
                  since we aren't using ReadPDB() so record widths are
                  different

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/MathUtil.h"
#include "bioplib/general.h"
#include "bioplib/fsscanf.h"

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
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
void Usage(void);
PDB *ReadSolv(FILE *fp, int *natom);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for moving accall atom accessibilities into B-val column

   05.07.94 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE *in       = stdin,
        *out      = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        PrevChain[8];
   int  natoms;
   PDB  *pdb, 
        *p;
   
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = ReadSolv(in,&natoms)) != NULL)
         {
            for(p=pdb; p!=NULL; NEXT(p))
            {
               if(strncmp(PrevChain,p->chain,1))
               {
                  /* Chain change, insert TER card                               */
                  fprintf(out,"TER   \n");
                  strcpy(PrevChain,p->chain);
               }
               WritePDBRecordAtnam(out, p);
            }
            fprintf(out,"TER   \n");
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *infile      Input file (or blank string)
            char   *outfile     Output file (or blank string)
   Returns: BOOL                Success?

   Parse the command line
   
   05.07.94 Original    By: ACRM
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

   05.07.94 Original    By: ACRM
   24.08.94 V1.1
   27.09.94 V1.2, corrected message
   30.05.02 V1.3
   03.05.04 V1.4
   30.09.05 V1.5
*/
void Usage(void)
{
   fprintf(stderr,"\nAS2BVal V1.5 (c) 1994-2005, Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Usage: as2bval [<in.pdb>] [<out.pdb>]\n");
   fprintf(stderr,"Rewrites the output from accall solvent accessibility \
as a standard PDB\n");
   fprintf(stderr,"format file with accessibility in the B-val column \
and radius in the\n");
   fprintf(stderr,"occupancy column.\n\n");
}

/************************************************************************/
/*>PDB *ReadSolv(FILE *fp, int *natom)
   -----------------------------------
   Input:   FILE   *fp     PDB file pointer
   Output:  int    *natom  Number of atoms read
   Returns: PDB    *       PDB linked list

   Reads a PDB-like file but with the OCC and BVAL columns widened for
   use by solvent accessibility data

   05.07.94 Original based on doReadPDB()   By: ACRM
   03.06.05 FixAtomName() now needs the occupancy
*/
PDB *ReadSolv(FILE *fp, int *natom)
{
   char     record_type[8],
            atnambuff[8],
            *atnam,
            resnam[8],
            chain[4],
            insert[4],
            buffer[160];
   int      atnum,
            resnum;
   double   x,y,z,
            occ,
            bval;
   PDB      *pdb  = NULL,
            *p;

   *natom         = 0;

   while(fgets(buffer,159,fp))
   {
      fsscanf(buffer,"%6s%5d%1x%5s%4s%1s%4d%1s%3x%8lf%8lf%8lf%8lf%6lf",
              record_type,&atnum,atnambuff,resnam,chain,&resnum,insert,
              &x,&y,&z,&bval,&occ);
      if(!strncmp(record_type,"ATOM  ",6) || 
         !strncmp(record_type,"HETATM",6))
      {
         /* Fix the atom name accounting for start in column 13 or 14   */
         atnam = FixAtomName(atnambuff, occ);
         
         /* Trim the atom name to 4 characters                          */
         atnam[4] = '\0';
         
         /* Allocate space in the linked list                           */
         if(pdb == NULL)
         {
            INIT(pdb, PDB);
            p = pdb;
         }
         else
         {
            ALLOCNEXT(p, PDB);
         }
         
         /* Failed to allocate space; free up list so far & return      */
         if(p==NULL)
         {
            if(pdb != NULL) FREELIST(pdb, PDB);
            *natom = (-1);
            return(NULL);
         }
         
         /* Increment the number of atoms                               */
         (*natom)++;
         
         /* Store the information read                                  */
         p->atnum  = atnum;
         p->resnum = resnum;
         p->x      = (REAL)x;
         p->y      = (REAL)y;
         p->z      = (REAL)z;
         p->occ    = (REAL)occ;
         p->bval   = (REAL)bval;
         p->next   = NULL;
         strcpy(p->record_type, record_type);
         strcpy(p->atnam,       atnam);
         strcpy(p->resnam,      resnam);
         strcpy(p->chain,       chain);
         strcpy(p->insert,      insert);
      }
   }

   /* Return pointer to start of linked list                            */
   return(pdb);
}

