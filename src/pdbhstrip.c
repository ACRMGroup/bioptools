/************************************************************************/
/**

   \file       pdbhstrip.c
   
   \version    V1.4
   \date       13.02.15
   \brief      Strip hydrogens from a PDB file. Acts as filter
   
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
-  V1.0  04.07.94 Original
-  V1.1  15.07.94 Writes TER cards and returns correctly
-  V1.2  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.3  06.11.14 Renamed from hstrip  By: ACRM
-  V1.4  13.02.15 Added whole PDB support and re-written to use
                  blStripHPDBAsCopy()

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
-  13.02.15 Added whole PDB support and re-written to use
            blStripHPDBAsCopy()  By: ACRM
*/
int main(int argc, char **argv)
{
   WHOLEPDB *wpdb;
   FILE     *in  = stdin, 
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

   if((wpdb=blReadWholePDB(in))!=NULL)
   {
      PDB *pdbin  = NULL,
          *pdbout = NULL;
      int natoms;
      pdbin  = wpdb->pdb;
      pdbout = blStripHPDBAsCopy(pdbin, &natoms);
      FREELIST(pdbin, PDB);
      wpdb->pdb = pdbout;

      blWriteWholePDB(out, wpdb);
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
-  06.11.14 V1.3 By: ACRM
-  13.02.15 V1.4
*/
void Usage(void)
{            
   fprintf(stderr,"\npdbhstrip V1.4 (c) 1994-2015, Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Usage: pdbhstrip [in.pdb [out.pdb]]\n\n");
   fprintf(stderr,"Removes hydrogens from a PDB file. I/O is through \
stdin/stdout if files\n");
   fprintf(stderr,"are not specified.\n\n");
}
