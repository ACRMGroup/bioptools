/************************************************************************/
/**

   \file       pdb2xyz.c
   
   \version    V1.1
   \date       22.07.14
   \brief      Convert PDB to Gromos XYZ
   
   \copyright  (c) Dr. Andrew C. R. Martin 1994-2014
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
-  V1.0  23.08.94 Original   By: ACRM
-  V1.1  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>

#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/general.h"

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
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  char *title);
void Usage(void);
void WriteXYZ(FILE *out, PDB *pdb, int natoms, char *title);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for format conversion

-  23.08.94 Original    By: ACRM
-  24.08.94 Changed to call OpenStdFiles()
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        title[MAXBUFF];
   PDB  *pdb;
   int  natoms;
   
   if(ParseCmdLine(argc, argv, infile, outfile, title))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDB(in,&natoms)) != NULL)
         {
            WriteXYZ(out, pdb, natoms, title);
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
                     char *title)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *title       Title for file
   \return                     Success?

   Parse the command line
   
-  23.08.94 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *title)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   title[0]  = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 't':
            argc--;
            argv++;
            strncpy(title,argv[0],MAXBUFF);
            title[MAXBUFF-1] = '\0';
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

-  23.08.94 Original    By: ACRM
-  22.07.14 V1.1 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\nPDB2XYZ V1.1 (c) 1994-2014, Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: pdb2xyz [-t title] [<in.pdb>] [<out.pdb>]\n\n");
   fprintf(stderr,"Convert PDB format to GROMOS XYZ. N.B. Does NOT \
correct atom order.\n\n");
}

/************************************************************************/
/*>void WriteXYZ(FILE *out, PDB *pdb, int natoms, char *title)
   -----------------------------------------------------------
*//**

   \param[in]      *out      Output file pointer
   \param[in]      *pdb      PDB linked list
   \param[in]      natoms    Number of atoms
   \param[in]      title     Title for XYZ file or blank string

   Write the PDB linked list out in Gromos XYZ format

-  23.08.94 Original    By: ACRM
*/
void WriteXYZ(FILE *out, PDB *pdb, int natoms, char *title)
{
   PDB *p;
   int i = 1;

   if(title[0])
      fprintf(out,"%s\n",title);
   else
      fprintf(out,"Gromos XYZ file generated by PDB2XYZ\n");

   fprintf(out,"%5d\n",natoms);
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      fprintf(out,"%5d %-4s %-4s %5d%8.3f%8.3f%8.3f\n",
              p->resnum,
              p->resnam,
              p->atnam,
              i++,
              p->x,
              p->y,
              p->z);
   }
}

