/************************************************************************/
/**

   \file       rotate.c
   
   \version    V1.3
   \date       22.07.14
   \brief      Simple program to rotate PDB files
   
   \copyright  (c) Dr. Andrew C. R. Martin / UCL 1994-7-2014
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
-  V1.0  17.06.94 Original
-  V1.1  21.07.95 Added -m option
-  V1.2  29.09.97 Added -n option
-  V1.3  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/matrix.h"
#include "bioplib/pdb.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
typedef struct _rotlist
{
   char direction;
   REAL angle;
   struct _rotlist *next;
}  ROTLIST;

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void DoRotations(PDB *pdb, ROTLIST *rotlist);
BOOL BuildRotInstruct(ROTLIST **pRotList, char direction, char *amount);
void Usage(void);
BOOL SetupMatrix(int argc, char **argv, REAL matrix[3][3]);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for PDB rotation

-  17.06.94 Original    By: ACRM
-  21.07.95 Added -m
-  29.09.97 Added -n
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE    *in      = stdin,
           *out     = stdout;
   ROTLIST *rotlist = NULL;
   PDB     *pdb;
   int     natoms;
   BOOL    GotMatrix = FALSE,
           GotRot    = FALSE,
           DoCentre  = TRUE;
   REAL    matrix[3][3];

   argc--;
   argv++;

   /* Handle all switches                                               */
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            Usage();
            return(0);
         case 'x':
            argc--;
            argv++;
            GotRot=TRUE;
            if(GotMatrix || !BuildRotInstruct(&rotlist, 'x', argv[0]))
            {
               Usage();
               return(1);
            }
            break;
         case 'y':
            argc--;
            argv++;
            GotRot=TRUE;
            if(GotMatrix || !BuildRotInstruct(&rotlist, 'y', argv[0]))
            {
               Usage();
               return(1);
            }
            break;
         case 'z':
            argc--;
            argv++;
            GotRot=TRUE;
            if(GotMatrix || !BuildRotInstruct(&rotlist, 'z', argv[0]))
            {
               Usage();
               return(1);
            }
            break;
         case 'm':
            argc--;
            argv++;
            GotMatrix = TRUE;
            if(GotRot || !SetupMatrix(argc, argv, matrix))
            {
               Usage();
               return(1);
            }
            argv+=8;
            argc-=8;
            break;
         case 'n':
            DoCentre = FALSE;
            break;
         default:
            Usage();
            return(1);
         }
      }
      else
      {
         break;
      }
      
      argc--;
      argv++;
   }

   if(!DoCentre && !GotMatrix)
   {
      fprintf(stderr,"rotate: Error -n may may only bed used with -m\n");
      return(1);
   }
   
   /* Check for excess file name specifications                         */
   if(argc > 2)
   {
      Usage();
      return(1);
   }
   
   /* Open input file if specified                                      */
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

   /* Open output file if specified                                     */
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
   
   /* Read in the PDB file                                              */
   if((pdb = blReadPDB(in, &natoms))==NULL)
   {
      fprintf(stderr,"No atoms read from PDB file\n");
      return(1);
   }

   /* Apply the rotations                                               */
   if(GotMatrix)
   {
      if(DoCentre)
         blRotatePDB(pdb, matrix);
      else
         blApplyMatrixPDB(pdb, matrix);
   }
   else
   {
      DoRotations(pdb,rotlist);
   }
   
   /* Write the new PDB file                                            */
   blWritePDB(out,pdb);
   
   return(0);
}

   
/************************************************************************/
/*>void DoRotations(PDB *pdb, ROTLIST *rotlist)
   --------------------------------------------
*//**

   \param[in,out]  *pdb        PDB linked list for rotation
   \param[in]      *rotlist    Linked list of rotation instructions

   Applies a set of rotation instructions to a PDB linked list

-  17.06.94 Original    By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
void DoRotations(PDB *pdb, ROTLIST *rotlist)
{
   ROTLIST *p;
   REAL    rotmat[3][3];
   
   for(p=rotlist; p!=NULL; NEXT(p))
   {
      blCreateRotMat(p->direction, p->angle, rotmat);
      blRotatePDB(pdb, rotmat);
   }
}


/************************************************************************/
/*>BOOL BuildRotInstruct(ROTLIST **pRotList, char direction, char *amount)
   ----------------------------------------------------------------------
*//**

   \param[in,out]  **pRotList   Pointer to linked list of instructions
   \param[in]      direction    Direction for rotation
   \param[in]      *amount      Amount by which to rotate (degrees)

   Add an item to a linked list of rotation insructions

-  17.06.94 Original    By: ACRM
*/
#define RotList *pRotList
BOOL BuildRotInstruct(ROTLIST **pRotList, char direction, char *amount)
{
   static ROTLIST  *p;
   
   if(RotList == NULL)
   {
      INIT((RotList), ROTLIST);
      p = RotList;
   }
   else
   {
      ALLOCNEXT(p, ROTLIST);
   }
   
   if(p==NULL) return(FALSE);
   
   p->direction = direction;
   if((sscanf(amount,"%lf",&(p->angle)))==0)
      return(FALSE);
   
   p->angle *= PI/180.0;
   
   return(TRUE);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

-  17.06.94 Original    By: ACRM
-  21.07.95 V1.1 Added -m
-  29.09.97 V1.2 Added -n
-  22.07.14 V1.3 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\nRotate V1.3  (c) 1994-2014 Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"Freely distributable if no profit is made\n\n");
   fprintf(stderr,"Usage: rotate [-m 11 12 13 21 22 23 31 32 33] [-h]\n");
   fprintf(stderr,"              [-n] [<input.pdb> [<output.pdb>]]\n");
   fprintf(stderr,"       --or--\n");
   fprintf(stderr,"       rotate [-x <ang>] [-y <ang>] [-z <ang>] \
[-h]\n");
   fprintf(stderr,"              [<input.pdb> [<output.pdb>]]\n\n");
   fprintf(stderr,"       -m           Specify rotation matrix\n");
   fprintf(stderr,"       -n           Do not move to CofG before applying\
 matrix\n");
   fprintf(stderr,"       -x, -y, -z   Specify rotations (in degrees)\n");
   fprintf(stderr,"       -h           This help message\n");
   fprintf(stderr,"I/O is to stdin/stdout if not specified\n\n");
   fprintf(stderr,"Rotates a PDB file using the given rotation matrix \
or using the sequence\n");
   fprintf(stderr,"of specified rotations. All rotations are performed \
around the centre\n");
   fprintf(stderr,"of geometry of the molecule. -x, -y and -z rotations \
are applied in\n");
   fprintf(stderr,"sequence and as many rotations are are required may \
be given.\n\n");
}


/************************************************************************/
/*>BOOL SetupMatrix(int argc, char **argv, REAL matrix[3][3])
   ----------------------------------------------------------
*//**

   Set up a rotation matrix from values on the command line

-  21.07.95 Original    By: ACRM
*/
BOOL SetupMatrix(int argc, char **argv, REAL matrix[3][3])
{
   int i, j;
   
   if(argc < 8)
      return(FALSE);
   
   for(i=0; i<3; i++)
   {
      for(j=0; j<3; j++)
      {
         if(sscanf(argv[0],"%lf",&(matrix[i][j]))==0)
         {
            return(FALSE);
         }
         argv++;
         if(--argc < 0)
         {
            return(FALSE);
         }
      }
   }
   
   return(TRUE);
}

