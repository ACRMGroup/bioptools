/************************************************************************/
/**

   \file       torsions.c
   
   \version    V1.4
   \date       19.08.14
   \brief      Generate a complete set of backbone torsion angles for a 
               protein.
   
   \copyright  (c) Andrew C.R. Martin 1994-5-2014
   \author     Dr. Andrew C. R. Martin
   \par
               Biomolecular Structure and Modelling Unit,
               Department of Biochemistry and Molecular Biology,
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

   Notes:
   ======
   
**************************************************************************

   Revision History:
   =================
-  V1.0  10.06.94 Original
-  V1.1  16.08.94 Added -m option for Martin's output format
-  V1.2  12.06.95 Added -c option for pseudo-CA torsions     
-  V1.3  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.4  19.08.14 Added AsCopy suffix to call to blSelectAtomsPDB() 
                  By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/angle.h"
#include "bioplib/seq.h"

/************************************************************************/
/* Defines
*/

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void ShowTorsions(FILE *out, PDB *pdb, REAL *tors, BOOL Radians,
                  BOOL MartinFormat);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for converting a PDB file to torsions.

-  10.06.94 Original   By: ACRM
-  16.08.94 Added -m option
-  12.06.95 Added -c option
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  19.08.14 Added AsCopy suffix to call to blSelectAtomsPDB() By: CTP
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   PDB  *fullpdb,
        *pdb,
        *p,
        *p1,
        *p2,
        *p3,
        *p4;
   int  natoms,
        TorNum;
   char *sel[4];
   REAL tors[3];
   BOOL Radians      = FALSE,
        MartinFormat = FALSE,
        CATorsions   = FALSE;

   argc--;
   argv++;
   
   /* Handle any switches                                               */
   if(argc)
   {
      while(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            Usage();
            return(0);
            break;
         case 'r':
            Radians = TRUE;
            break;
         case 'm':
            MartinFormat = TRUE;
            break;
         case 'c':
            CATorsions = TRUE;
            break;
         default:
            Usage();
            return(1);
            break;
         }
         argc--;
         argv++;
      }
   }

   if(argc > 2) 
   {
      Usage();
      return(1);
   }

   /* Handle any filenames                                              */
   if(argc)
   {
      if((in=fopen(argv[0],"r"))==NULL)
      {
         fprintf(stderr,"Unable to open input file: %s\n",argv[0]);
         return(1);
      }
      
      argc--;
      argv++;

      if(argc)
      {
         if((out=fopen(argv[0],"w"))==NULL)
         {
            fprintf(stderr,"Unable to open output file: %s\n",argv[0]);
            return(1);
         }
      }
   }
   
   /* Read in the structure                                             */
   if((fullpdb = blReadPDB(in, &natoms))==NULL)
   {
      fprintf(stderr,"No atoms read from PDB file\n");
      return(1);
   }

   /* Set up the atom selection and select them                         */
   SELECT(sel[0],"CA  ");
   if(!CATorsions)
   {
      SELECT(sel[1],"N   ");
      SELECT(sel[2],"C   ");
   }
   
   if((pdb = blSelectAtomsPDBAsCopy(fullpdb,(CATorsions?1:3),sel,&natoms))
      == NULL)
   {
      fprintf(stderr,"Unable to select backbone atoms from PDB \
file (no memory?)\n");
      return(1);
   }

   /* Print title                                                       */
   if(CATorsions)
   {
      fprintf(out,"Res_N    CA_N--CA_(N+1)\n");
   }
   else
   {
      
      if(MartinFormat)
         fprintf(out,"Residue    PHI      PSI      OMEGA\n");
      else
         fprintf(out,"               PHI      PSI     OMEGA\n");
   }

   fprintf(out,"--------------------------------------\n");
   

   /* Walk the linked list and calculate torsions                       */
   tors[0] = tors[1] = tors[2] = 9999.0;
   p1      = p2      = p3      = p4      = NULL;

   if(CATorsions)
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         p1 = p2;
         p2 = p3;
         p3 = p4;
         p4 = p;
         if(p1 && p2 && p3 && p4)   /* Got all 4 atoms                  */
         {
            tors[0] = blPhi(p1->x, p1->y, p1->z,
                            p2->x, p2->y, p2->z,
                            p3->x, p3->y, p3->z,
                            p4->x, p4->y, p4->z);
            if(!Radians) tors[0] *= 180.0 / PI;
            fprintf(out,"   %c    %8.3f\n",blThrone(p2->resnam),tors[0]);
         }
         else if(p1==NULL && p2==NULL && p3==NULL && p4) /* Got 1 atom  */
         {
            fprintf(out,"   %c        -\n",blThrone(p4->resnam));
         }
      }
      /* Finish off by printing the last 2 residues                     */
      fprintf(out,"   %c        -\n",blThrone(p3->resnam));
      fprintf(out,"   %c        -\n",blThrone(p4->resnam));
   }
   else
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(!strncmp(p->atnam,"C   ",4))
         {
            ShowTorsions(out, p, tors, Radians, MartinFormat);
            tors[0] = tors[1] = tors[2] = 9999.0;
         }
         
         /* Get pointers to four atoms in sequence                      */
         p1 = p;
         p2 = p->next;
         if(p2 != NULL) p3 = p2->next;
         if(p3 != NULL) p4 = p3->next;
         
         if(p1==NULL || p2==NULL || p3==NULL || p4==NULL)
         {
            ShowTorsions(out, p, tors, Radians, MartinFormat);
            break;
         }
         
         if(!strncmp(p->atnam,"N   ",4))
            TorNum = 1;
         else if(!strncmp(p->atnam,"CA  ",4))
            TorNum = 2;
         else if(!strncmp(p->atnam,"C   ",4))
            TorNum = 0;
         
         tors[TorNum] = blPhi(p1->x, p1->y, p1->z,
                              p2->x, p2->y, p2->z,
                              p3->x, p3->y, p3->z,
                              p4->x, p4->y, p4->z);
      }
   }
   return(0);
}

/************************************************************************/
/*>void ShowTorsions(FILE *out, PDB *pdb, REAL *tors, BOOL Radians,
                     BOOL MartinFormat)
   ----------------------------------------------------------------
*//**

   \param[in]      *out          Output file
   \param[in]      *pdb          PDB record pointer
   \param[in]      *tors         Array of torsion angles
   \param[in]      Radians       Should output be in radians
   \param[in]      MartinFormat  Output in Martin's required format

   Displays the torsion angles converting to degrees if Radians flag
   not set.

-  10.06.94 Original    By: ACRM
-  16.08.94 Added MartinFormat
*/
void ShowTorsions(FILE *out, PDB *pdb, REAL *tors, BOOL Radians,
                  BOOL MartinFormat)
{
   if(!Radians)
   {
      int i;

      for(i=0; i<3; i++)
      {
         if(tors[i] < (REAL)9990.0)
            tors[i] *= (REAL)180.0/PI;
      }
      
   }
   
   if(MartinFormat)
   {
      fprintf(out,"   %c    ", blThrone(pdb->resnam));

      if(tors[0] > (REAL)9998.0)
      {
         fprintf(out,"    -    %8.3f %8.3f\n",tors[1],tors[2]);
      }
      else if(tors[1] > (REAL)9998.0)
      {
         fprintf(out,"%8.3f     -        -   \n",tors[0]);
      }
      else
      {
         fprintf(out,"%8.3f %8.3f %8.3f\n",tors[0], tors[1], tors[2]);
      }
   }
   else
   {
      fprintf(out,"%5d%c %-4s %8.3f %8.3f %8.3f\n",
              pdb->resnum, pdb->insert[0], pdb->resnam,
              tors[0], tors[1], tors[2]);
   }
}
   
/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Displays a usage message

-  10.06.94 original   By: ACRM
-  16.08.94 Added -m
-  12.06.95 Added -c
-  22.07.14 V1.3 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\ntorsions V1.3. (c) 1994-2014 Andrew Martin, UCL. \
Freely Distributable\n");
   fprintf(stderr,"Generates a set of backbone torsions from a PDB \
file.\n\n");
   fprintf(stderr,"Usage: torsions [-h] [-r] [<in.pdb> [<out.tor>]]\n");
   fprintf(stderr,"       -h   This help message\n");
   fprintf(stderr,"       -r   Give results in radians\n");
   fprintf(stderr,"       -m   Output format required by Martin \
Reczko\n");
   fprintf(stderr,"       -c   Generate CA-CA pseudo-torsions\n");
   fprintf(stderr,"I/O is to stdin/stdout if unspecified.\n\n");
}

