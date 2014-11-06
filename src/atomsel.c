/*************************************************************************

   Program:    atomsel
   File:       atomsel.c
   
   Version:    V1.1
   Date:       24.08.94
   Function:   Select atoms from a PDB file. Acts as filter
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1994
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   Phone:      (Home) +44 (0372) 275775
   EMail:      INTERNET: amartin@scitec.adsp.sub.org
                         martin@bsm.bioc.ucl.ac.uk
               UUCP:     ...{uunet|rutgers}!cbmehq!cbmuk!scitec!amartin
               JANET:    martin@uk.ac.ucl.bioc.bsm
               
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
   V1.0  15.07.94 Original    By: ACRM
   V1.1  24.08.94 Changed to call OpenStdFiles()

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160

typedef struct _atomtype
{
   struct _atomtype *next;
   char type[8];
}  ATOMTYPE;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, ATOMTYPE **atoms, char *infile, 
                  char *outfile);
void UpcaseTypes(ATOMTYPE *atoms);
BOOL InAtomList(ATOMTYPE *atoms, PDB *p);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for selecting atoms

   04.07.94 Original    By: ACRM
   24.08.94 Changed to call OpenStdFiles()
*/
int main(int argc, char **argv)
{
   PDB      *pdb, 
            *p;
   int      natom;
   char     infile[MAXBUFF],
            outfile[MAXBUFF],
            lastchain;
   FILE     *in  = stdin, 
            *out = stdout;
   ATOMTYPE *atoms = NULL;

   if(ParseCmdLine(argc, argv, &atoms, infile, outfile))
   {
      if(OpenStdFiles(infile, outfile, &in, &out))
      {
         /* If no atoms specified, assume CA                            */
         if(atoms == NULL)
         {
            INIT(atoms, ATOMTYPE);
            if(atoms == NULL)
            {
               fprintf(stderr,"Error: Unable to allocate memory for atom \
list.\n");
               return(1);
            }

            strcpy(atoms->type,"CA  ");
         }

         /* Upcase all atom names                                       */
         UpcaseTypes(atoms);

         /* Read in the PDB file                                        */
         if((pdb=ReadPDB(in,&natom))!=NULL)
         {
            lastchain = pdb->chain[0];
            
            /* Run through the linked list writing out atoms in the
               required atom list
            */
            for(p=pdb;p!=NULL;NEXT(p))
            {
               if(p->chain[0] != lastchain)
               {
                  lastchain = p->chain[0];
                  fprintf(out,"TER   \n");
               }
                  
               if(InAtomList(atoms, p))
                  WritePDBRecord(out,p);
            }
            fprintf(out,"TER   \n");
         }
         else
         {
            fprintf(stderr,"Warning: No atoms read from PDB file.\n");
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
   Prints a usage message

   15.07.94 Original    By: ACRM
   24.08.94 V1.1
*/
void Usage(void)
{            
   fprintf(stderr,"\nAtomSel V1.1 (c) 1994, Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: atomsel [-atom] [-atom...] [<in.pdb>] \
[<out.pdb>]\n\n");
   fprintf(stderr,"Selects specfied atom types from a PDB file. \
Assumes C-alpha if no atoms\n");
   fprintf(stderr,"are specified. I/O is through stdin/stdout if files \
are not specified.\n\n");
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, ATOMTYPE **atoms,
                     char *infile, char *outfile)
   ---------------------------------------------------------------------
   Input:   int      argc        Argument count
            char     **argv      Argument array
   Output:  ATOMTYPE **atoms     Linked list of atoms to keep
            char     *infile     Input filename (or blank string)
            char     *outfile    Output filename (or blank string)
   Returns: BOOL                 Success

   Parse the command line

   15.07.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, ATOMTYPE **atoms, char *infile, 
                  char *outfile)
{
   char     buffer[80],
            *chp,
            *chq;
   ATOMTYPE *a;
   

   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         if(!strncmp(argv[0]+1,"help",4))
            return(FALSE);
         
         if(*atoms == NULL)
         {
            INIT((*atoms),ATOMTYPE);
            a = *atoms;
         }
         else
         {
            ALLOCNEXT(a, ATOMTYPE);
         }
         
         if(a == NULL)
         {
            fprintf(stderr,"Error: Unable to allocate memory for atom \
list.\n");
            exit(1);
         }
         
         strcpy(a->type, argv[0]+1);
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
/*>void UpcaseTypes(ATOMTYPE *atoms)
   ---------------------------------
   I/O:     ATOMTYPE   *atoms    Linked list of atom types

   Upcases the atom types in the specified atom list and pads all names
   to four characters.

   15.07.94 Original   By: ACRM
*/
void UpcaseTypes(ATOMTYPE *atoms)
{
   ATOMTYPE *a;
   
   for(a=atoms; a!=NULL; NEXT(a))
   {
      padterm(a->type, 4);
      UPPER(a->type);
   }
}

/************************************************************************/
/*>BOOL InAtomList(ATOMTYPE *atoms, PDB *p)
   ----------------------------------------
   Input:   ATOMTYPE   *atoms    Linked list of required atom types
            PDB        *p        PDB item to check against list
   Returns: BOOL                 In list?

   Compares an atom name from a PDB item with the linked list of allowed
   atom names.

   15.07.94 Original   By: ACRM
*/
BOOL InAtomList(ATOMTYPE *atoms, PDB *p)
{
   ATOMTYPE *a;
   

   for(a=atoms; a!=NULL; NEXT(a))
   {
      if(!strncmp(p->atnam, a->type, 4))
         return(TRUE);
   }
   
   return(FALSE);
}

