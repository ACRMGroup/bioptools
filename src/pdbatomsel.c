/************************************************************************/
/**

   \file       pdbatomsel.c
   
   \version    V1.7
   \date       02.03.15
   \brief      Select atoms from a PDB file. Acts as filter
   
   \copyright  (c) Dr. Andrew C. R. Martin 1994-2015
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
-  V1.0  15.07.94 Original    By: ACRM
-  V1.1  24.08.94 Changed to call OpenStdFiles()
-  V1.2  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.3  19.08.14 Removed unused variables in ParseCmdLine() By: CTP
-  V1.4  06.11.14 Renamed from atomsel By: ACRM
-  V1.5  07.11.14 Initialized a variable
-  V1.6  12.02.15 Uses Whole PDB
-  V1.7  02.03.15 Major rewrite to use blSelectAtomsPDBAsCopy()

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
#include "bioplib/array.h"

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
void UpcaseAndPadAtomTypes(ATOMTYPE *atoms);
char **ConvertAtomsToArray(ATOMTYPE *atoms, int *nSelected);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for selecting atoms

-  04.07.94 Original    By: ACRM
-  24.08.94 Changed to call OpenStdFiles()
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  12.02.15 WholePDB support   By: ACRM
            Better support for TER cards
-  02.03.15 Major rewrite to use blSelectAtomsPDBAsCopy()
*/
int main(int argc, char **argv)
{
   WHOLEPDB *wpdb;
   char     infile[MAXBUFF],
            outfile[MAXBUFF],
            **selectedAtoms;
   FILE     *in  = stdin, 
            *out = stdout;
   ATOMTYPE *atoms = NULL;
   int      nSelected;

   if(ParseCmdLine(argc, argv, &atoms, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
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
         UpcaseAndPadAtomTypes(atoms);
         
         /* Convert atoms linked list to select array                   */
         selectedAtoms = ConvertAtomsToArray(atoms, &nSelected);
         FREELIST(atoms, ATOMTYPE);

         /* Read in the PDB file                                        */
         if((wpdb=blReadWholePDB(in))!=NULL)
         {
            PDB *pdb = NULL;
            int natoms;
            
            pdb = blSelectAtomsPDBAsCopy(wpdb->pdb, nSelected,
                                         selectedAtoms, &natoms);
            FREELIST(wpdb->pdb, PDB);
            wpdb->pdb    = pdb;
            wpdb->natoms = natoms;

            blWriteWholePDB(out, wpdb);
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
/*>char **ConvertAtomsToArray(ATOMTYPE *atoms, int *nSelected)
   -----------------------------------------------------------
*//**
   \param[in]   *atoms       ATOMTYPE linked list
   \param[out]  *nSelected   Number of selected atoms in the array
   \return                   2D array of atom types. Each is padded to
                             4 chars. NULL on error.

   Converts the linked list of strings to an array of strings

-  02.03.15  Original   By: ACRM
*/
char **ConvertAtomsToArray(ATOMTYPE *atoms, int *nSelected)
{
   ATOMTYPE *a;
   char     **selectedAtoms = NULL;
   int      i;

   /* Count the atom types                                              */
   *nSelected = 0;
   for(a=atoms; a!=NULL; NEXT(a))
      (*nSelected)++;

   /* Allocate a 2D array to contain the list                           */
   if((selectedAtoms = (char **)blArray2D(sizeof(char), 
                                          *nSelected, 5))==NULL)
      return(NULL);
   
   /* Populate the array                                                */
   for(a=atoms, i=0; a!=NULL; NEXT(a), i++)
   {
      strncpy(selectedAtoms[i], a->type, 4);
      selectedAtoms[i][5] = '\0';
      PADMINTERM(selectedAtoms[i], 4);
   }

   return(selectedAtoms);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  15.07.94 Original    By: ACRM
-  24.08.94 V1.1
-  22.07.14 V1.2 By: CTP
-  06.11.14 V1.4 By: ACRM
-  07.11.14 V1.5
-  12.02.15 V1.6
-  02.03.15 V1.7
*/
void Usage(void)
{            
   fprintf(stderr,"\npdbatomsel V1.7 (c) 1994-2015, Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"Usage: pdbatomsel [-atom] [-atom...] [<in.pdb> \
[<out.pdb>]]\n\n");
   fprintf(stderr,"Selects specfied atom types from a PDB file. \
Assumes C-alpha if no atoms\n");
   fprintf(stderr,"are specified. I/O is through stdin/stdout if files \
are not specified.\n\n");
   fprintf(stderr,"Note that this program does not currently support \
PDBML output\n\n");
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, ATOMTYPE **atoms,
                     char *infile, char *outfile)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     **atoms     Linked list of atoms to keep
   \param[out]     *infile     Input filename (or blank string)
   \param[out]     *outfile    Output filename (or blank string)
   \return                     Success

   Parse the command line

-  15.07.94 Original    By: ACRM
-  19.08.14 Removed unused variables. By: CTP
-  07.11.14 Initialized a
*/
BOOL ParseCmdLine(int argc, char **argv, ATOMTYPE **atoms, char *infile, 
                  char *outfile)
{
   ATOMTYPE *a = NULL;

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
/*>void UpcaseAndPadAtomTypes(ATOMTYPE *atoms)
   -------------------------------------------
*//**

   \param[in,out]  *atoms    Linked list of atom types

   Upcases the atom types in the specified atom list and pads all names
   to four characters.

-  15.07.94 Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  02.03.15 Added initial termination of string
*/
void UpcaseAndPadAtomTypes(ATOMTYPE *atoms)
{
   ATOMTYPE *a;
   
   for(a=atoms; a!=NULL; NEXT(a))
   {
      a->type[5] = '\0';
      blPadterm(a->type, 4);
      UPPER(a->type);
   }
}

