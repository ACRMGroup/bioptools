/************************************************************************/
/**

   \file       pdbdummystrip.c
   
   \version    V1.4
   \date       09.11.21
   \brief      Strips atoms with NULL coordinates
   
   \copyright  (c) Prof. Andrew C. R. Martin 1996-2021
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
-  V1.0  13.11.96 Original    By: ACRM
-  V1.1  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.2  06.11.14 Renamed from nullstrip. This replaces an older program
                  called pdbstrip.   By: ACRM
-  V1.3  13.02.15 Added whole PDB support
-  V1.4  09.11.21 Now checks for >= 9999.0 instead of ==

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
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
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile);
PDB *StripNulls(PDB *pdb);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for stripping NULL coordinates

-  03.11.96 Original    By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE     *in  = stdin,
            *out = stdout;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];
   WHOLEPDB *wpdb;
   PDB      *pdb;
   
   if(ParseCmdLine(argc, argv, infile, outfile))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb = blReadWholePDB(in)) != NULL)
         {
            pdb=wpdb->pdb;
            pdb = StripNulls(pdb);
            blWriteWholePDB(out, wpdb);
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
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \return                     Success?

   Parse the command line
   
-  13.11.96 Original    By: ACRM
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
/*>PDB *StripNulls(PDB *pdb)
   -------------------------
*//**

   Strip NULL coordinate atoms from a PDB linked list

-  13.11.96 Original   By: ACRM
-  09.11.21 CHanged from == to >= 9999.0
*/
PDB *StripNulls(PDB *pdb)
{
   PDB *p, *prev;
   
   /* Remove from start                                                 */
   for(p=pdb; p!=NULL;)
   {
      if((p->x >= (REAL)9999.0) &&
         (p->y >= (REAL)9999.0) &&
         (p->z >= (REAL)9999.0))
      {
         pdb = p->next;
         free(p);
      }
      else
      {
         break;
      }
      p=pdb;
   }

   /* Remove rest                                                       */
   prev = pdb;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if((p->x >= (REAL)9999.0) &&
         (p->y >= (REAL)9999.0) &&
         (p->z >= (REAL)9999.0))
      {
         prev->next = p->next;
         free(p);
         p=prev;
      }
      prev = p;
   }
   
   return(pdb);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  13.11.96 Original    By: ACRM
-  22.07.14 V1.1 By: CTP
-  06.11.14 V1.2 By: ACRM
-  13.02.15 V1.4 By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\npdbdummystrip V1.4 (c) 1996-2021, Prof. Andrew \
C.R. Martin, UCL\n");

   fprintf(stderr,"\nUsage: pdbdummystrip [in.pdb [out.pdb]]\n");

   fprintf(stderr,"\nRemoves atoms from a PDB file which have NULL \
coordinates (i.e.\n");
   fprintf(stderr,"x,y,z >= 9999.0)\n");
}
