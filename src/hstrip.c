/************************************************************************/
/**

   \file       hstrip.c
   
   \version    V1.2
   \date       22.07.14
   \brief      Strip hydrogens from a PDB file. Acts as filter
   
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
-  V1.0  04.07.94 Original
-  V1.1  15.07.94 Writes TER cards and returns correctly
-  V1.2  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"


/************************************************************************/
/* Defines and macros
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for stripping hydrogens

-  04.07.94 Original    By: ACRM
-  15.07.94 Now writes TER cards and returns 0 correctly
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   PDB  *pdb, 
        *p;
   int  natom;
   FILE *in  = stdin, 
        *out = stdout;
   char lastchain;

   argc--;
   argv++;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            Usage();
            return(0);
            break;
         default:
            Usage();
            return(1);
            break;
         }
      }
      else
      {
         break;
      }
      
      argc--;
      argv++;
   }

   if(argc)
   {
      if((in=fopen(argv[0],"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input file: %s\n",argv[0]);
         return(1);
      }
      argc--;
      argv++;
   }
   
   if(argc)
   {
      if((out=fopen(argv[0],"w"))==NULL)
      {
         fprintf(stderr,"Unable to open output file: %s\n",argv[0]);
         return(1);
      }
      argc--;
      argv++;
   }

   if(argc)
   {
      Usage();
      return(1);
   }

   if((pdb=blReadPDB(in,&natom))!=NULL)
   {
      lastchain = pdb->chain[0];
      
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(p->atnam[0] != 'H')
         {
            if(p->chain[0] != lastchain)
            {
               lastchain = p->chain[0];
               fprintf(out,"TER   \n");
            }
            blWritePDBRecord(out,p);
         }
      }
      fprintf(out,"TER   \n");
   }

   return(0);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  04.07.94 Original    By: ACRM
-  22.07.14 V1.2 By: CTP
*/
void Usage(void)
{            
   fprintf(stderr,"\nHStrip V1.2 (c) 1994-2014, Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Usage: hstrip [<in.pdb>] [<out.pdb>]\n\n");
   fprintf(stderr,"Removes hydrogens from a PDB file. I/O is through \
stdin/stdout if files\n");
   fprintf(stderr,"are not specified.\n\n");
}
