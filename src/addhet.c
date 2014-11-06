/***********************************************************************/
/* Includes */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"

/***********************************************************************/
/* Prototypes */

int main (int argc, char *argv[]);
void DetermineBoundingBox (PDB *pdb, 
                           REAL *xmin, REAL *xmax, 
                           REAL *ymin, REAL *ymax, 
                           REAL *zmin, REAL *zmax);
PDB *ReadPDBHetAtoms(FILE *fp1, int *natoms);
void PrintBoundedHets(FILE *fp3, PDB *pdb, PDB *pdb2, 
                      REAL xmin, REAL xmax, 
                      REAL ymin, REAL ymax, 
                      REAL zmin, REAL zmax);

/***********************************************************************/
#define MAXDISTSQ 36

/***********************************************************************/
/* OK */
int main (int argc, char *argv[])
{
  
   FILE *fp1;
   FILE *fp2;
   FILE *fp3;
   
   PDB *pdbDomain;
   PDB *pdbHetatm;
   int natoms;
   
   REAL xmin, xmax, ymin, ymax, zmin, zmax;
   
   /* check correct number of files are specified on command line */
   if(argc !=4)
   {
      fprintf(stderr, "Usage: addhet complete.pdb domain.pdb \
domhet.pdb\n");
      exit(1);
   }
   
   if((fp1= fopen(argv[1], "r")) == NULL)
   {
      fprintf(stderr, "Error opening pdb file\n");
      exit(1);
   }
   
   if((fp2= fopen(argv[2], "r")) == NULL)
   {
      fprintf(stderr, "Error opening domain file\n");
      exit(1);
   }
   
   if((fp3= fopen(argv[3], "w")) == NULL)
   {
      fprintf(stderr, "Error opening output file\n");
      exit(1);
   }
   
   
   if((pdbDomain =  ReadPDB(fp2, &natoms))!=NULL)
   {
      DetermineBoundingBox(pdbDomain, &xmin, &xmax, &ymin, &ymax,  
                           &zmin, &zmax);
      WritePDB(fp3, pdbDomain);
      if((pdbHetatm=ReadPDBHetAtoms(fp1, &natoms))!=NULL)
      {
         PrintBoundedHets(fp3, pdbDomain, pdbHetatm, xmin, xmax, 
                          ymin, ymax, zmin, zmax);
      }
      /* This isn't always an error (there may have been no HETATMs)   */
/* ACRM--- 03.07.02
//    else
//    {
//       fprintf(stderr, "Error reading hetatoms\n");
//    }
*/
   }
   else
   {
      fprintf(stderr, "Error reading PDB\n");
   }
   return(0);
   
}

/***********************************************************************/
/* OK */
void DetermineBoundingBox (PDB *pdb, REAL *xmin, REAL *xmax, 
                           REAL *ymin, REAL *ymax, 
                           REAL *zmin, REAL *zmax)
{
   PDB *p;
   *xmin = *xmax = pdb->x;
   *ymin = *ymax = pdb->y;
   *zmin = *zmax = pdb->z;
   
   
   for (p=pdb; p!=NULL; NEXT(p))
   {
      if(p->x < *xmin)
      {
         *xmin = p->x;
      }
      if(p->x > *xmax)
      {
         *xmax = p->x;
      }
      if(p->y < *ymin)
      {
         *ymin = p->y;
      }
      if(p->y > *ymax)
      {
         *ymax = p->y;
      }
      if(p->z < *zmin)
      {
         *zmin = p->z;
      }
      if(p->z > *zmax)
      {
         *zmax = p->z;
      }
   }  
}

/***********************************************************************/
/* OK */
PDB *ReadPDBHetAtoms(FILE *fp1, int *natoms)
{
   PDB *pdb = NULL, 
       *pdbHetatm = NULL, 
       *p, *q;
   int natoms2;
       *natoms = 0;
   
   if((pdb = ReadPDB(fp1, &natoms2))!=NULL)
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(!strncmp(p->record_type, "HETATM", 6))
         {
            if(pdbHetatm==NULL)
            {
               INIT(pdbHetatm, PDB);
               q = pdbHetatm;
            }
            else
            {
               ALLOCNEXT(q, PDB);
            }
            if(q==NULL)
            {
               FREELIST(pdb, PDB);
               FREELIST(pdbHetatm, PDB);
               /* ACRM+++ 03.07.02 Report error here and exit(1) as 
                  NULL may be a valid return
//                return(NULL);
               */
               fprintf(stderr,"No memory for HETATM list\n");
               exit(1);
               /* ACRM-END 03.07.02                                    */
            }
            
            CopyPDB(q, p);
            (*natoms)++;
         }
      }
      FREELIST(pdb, PDB);
      return(pdbHetatm);
   }
   return(NULL);
}

/***********************************************************************/
/*>void PrintBoundedHets(FILE *fp3, PDB *pdb, PDB *pdbHetatm, 
                         REAL xmin, REAL xmax, 
                         REAL ymin, REAL ymax, 
                         REAL zmin, REAL zmax)
   -----------------------------------------------------------
   Input:   FILE   *fp3       File pointer for output
            PDB    *pdb       PDB linked list of protein atoms
            PDB    *pdbHetatm PDB linked list of HET atoms

   blah blah blah

   08.07.02 Original By: ALC
*/
void PrintBoundedHets(FILE *fp3, PDB *pdb, PDB *pdbHetatm, 
                      REAL xmin, REAL xmax, 
                      REAL ymin, REAL ymax, 
                      REAL zmin, REAL zmax)
{
   PDB  *start, 
        *stop, 
        *p, *q;
   BOOL isNeighbour,
        isNonStdAA,
        inBounds,
        clash;
   REAL distsq;
   int  bbcount;
   
   
   for(start=pdbHetatm; start!=NULL; start=stop)
   {
      stop = FindNextResidue(start);
      isNeighbour = FALSE;
      clash       = FALSE;
      isNonStdAA  = FALSE;
      bbcount     = 0;
      inBounds    = FALSE;
      
      /* See if any of the atoms of this HET group are in the bounds   */
      for(p=start; p!=stop; NEXT(p))
      {
         if((p->x>=xmin) && (p->x<=xmax) && 
            (p->y>=ymin) && (p->y<=ymax) && 
            (p->z>=zmin) && (p->z<=zmax))
         {
            inBounds = TRUE;
            break;
         }
      }
      
      if(inBounds)
      {
         /* See if this is actually a non-standard amino acid rather than
            a true HET group
         */
         for(p=start; p!=stop; NEXT(p))
         {
            if(!strncmp(p->atnam, "N   ", 4) ||
               !strncmp(p->atnam, "CA  ", 4) ||
               !strncmp(p->atnam, "C   ", 4) ||
               !strncmp(p->atnam, "O   ", 4))
            {
               if(++bbcount == 4)
               {
                  isNonStdAA = TRUE;
                  break;
               }
            }
         }
         
         /* If it's a standard AA then check if any atoms are in range */
         if(!isNonStdAA)
         {
            for(p=start; p!=stop; NEXT(p))
            {
               /* Look and see if any atom is close enough but not
                  clashing
               */
               for(q = pdb; q!=NULL; NEXT(q))
               {
                  distsq = DISTSQ(p,q);
                  if(distsq < MAXDISTSQ)
                  {                  
                     isNeighbour = TRUE;
                     /* It has been picked up already as it clashes    */
                     if(distsq < 0.1)
                     {
                        clash = TRUE;
                        break;
                     }
                  }
               }
               if(clash)
               {
                  break;
               }
            }

            /* If it's a neighbour and it doesn't clash, print it      */
            if(isNeighbour && !clash)
            {
               for(p=start; p!=stop; NEXT(p))
               {
                  WritePDBRecord(fp3, p);
               }
            }
         }  /* Is not a non-standard AA group                          */
      }  /* Is within boundaries                                       */
   }  /* loop over HET groups                                          */
}

/***********************************************************************/

  





