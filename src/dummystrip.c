/************************************************************************/
/**

   \file       dummystrip.c
   
   \version    V1.1
   \date       22.07.14
   \brief      Remove dummy atoms from a PDB file
   
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
   
-  V1.1  22.07.14 Renamed deprecated functions with bl prefix.
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

   Main program for stripping dummy atoms

-  04.07.94 Original    By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   PDB *pdb, *p;
   int natom;
   FILE *in  = stdin, 
        *out = stdout;

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
      for(p=pdb;p!=NULL;NEXT(p))
      {
         if(!(p->x > 9998.9 && p->x < 9999.1) &&
             (p->y > 9998.9 && p->y < 9999.1) &&
             (p->z > 9998.9 && p->z < 9999.1))
            blWritePDBRecord(out,p);
      }
   }
}

/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  04.07.94 Original    By: ACRM
-  22.07.14 V1.1 By: CTP
*/
void Usage(void)
{            
   fprintf(stderr,"\nDummyStrip V1.1 (c) 1994-2014, Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Usage: dummystrip [<in.pdb>] [<out.pdb>]\n\n");
   fprintf(stderr,"Removes dummy atoms from a PDB file. I/O is through \
stdin/stdout if files\n");
   fprintf(stderr,"are not specified.\n\n");
}
