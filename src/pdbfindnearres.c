/************************************************************************/
/**

   Program:    pdbfindnearres
   \file       pdbfindnearres.c
   
   \version    V1.0.2
   \date       29.06.19       
   \brief      Finds residues of a specified type near to the given
               zones   
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2019
   \author     Prof. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

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
   V1.0   05.06.19  Original
   V1.0.1 18.06.19  Fixed buffer size for fussy compiler
   V1.0.2 19.06.19  Fully documented and default radius in help message

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bioplib/macros.h"
#include "bioplib/pdb.h"
#include "bioplib/array.h"


/************************************************************************/
/* Defines and macros
*/
#define SMALLBUFF 32
#define MAXBUFF   160
#define DEFRAD    8.0

typedef struct _zone
{
   char start[SMALLBUFF],
        stop[SMALLBUFF];
   struct _zone *next;
} ZONE;

#define STRNCPY(a,b,c) \
do {  strncpy(a,b,c);  \
   a[c-1] = '\0';      \
} while(0)
#define ISBACKBONE(z)  (!strncmp((z)->atnam,"N   ",4) || \
                        !strncmp((z)->atnam,"CA  ",4) || \
                        !strncmp((z)->atnam,"C   ",4) || \
                        !strncmp((z)->atnam,"O   ",4))
#define ISSIDECHAIN(z) !ISBACKBONE(z)



/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void WriteFlaggedAtoms(FILE *out, PDB *pdb);
void ListFlaggedResidues(FILE *out, PDB *pdb);
BOOL ResInZone(PDB *res, ZONE *zones);
void ClearOccup(PDB *pdb);
void FlagNearRes(PDB *pdb, ZONE *zones, char *restype, REAL radiusSq);
void SetOccup(PDB *start, PDB *stop);
BOOL CheckInRange(PDB *res1, PDB *nextRes1, PDB *res2, PDB *nextRes2, 
                  REAL radiusSq);
ZONE *ParseZoneSpec(char *zonespec);
void PopulateZone(ZONE *z, char *zoneDescription);
char **blSplitStringOnCharacter(char *string, char charac, 
                                int minItemLen);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *radius, char *zonespec, char *restype, 
                  BOOL *listOnly);
void Usage(void);

/************************************************************************/
int main(int argc, char **argv)
{
   FILE *in      = stdin,
        *out     = stdout;
   REAL radius   = DEFRAD;
   PDB  *pdb     = NULL;
   ZONE *zones   = NULL;
   BOOL listOnly = FALSE;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        zonespec[MAXBUFF],
        restype[MAXBUFF];
   
   
   if(ParseCmdLine(argc, argv, infile, outfile, &radius, zonespec, 
                   restype, &listOnly))
   {
      if((zones = ParseZoneSpec(zonespec))!=NULL)
      {
#ifdef DEBUG
         ZONE *z;
         for(z=zones; z!=NULL; NEXT(z))
         {
            fprintf(stderr, "ZONE: %s to %s\n", z->start, z->stop);
         }
#endif

         /* Square the radius to save on distance sqrt()s               */
         radius *= radius;
         
         if(blOpenStdFiles(infile, outfile, &in, &out))
         {
            WHOLEPDB *wpdb;
            if((wpdb = blReadWholePDB(in)) != NULL)
            {
               pdb = wpdb->pdb;
               FlagNearRes(pdb, zones, restype, radius);
               if(listOnly)
               {
                  ListFlaggedResidues(out, pdb);
               }
               else
               {
                  WriteFlaggedAtoms(out, pdb);
               }
            }
            else
            {
               fprintf(stderr,"No atoms read from PDB file\n");
            }
         }
      }
      else
      {
         fprintf(stderr,"Unable to parse zones\n");
      }
   }
   else
   {
      Usage();
   }

   return(0);
}


/************************************************************************/
/*>void WriteFlaggedAtoms(FILE *out, PDB *pdb)
   -------------------------------------------
*//**
   \param[in]    out   Output file pointer
   \param[in]    pdb   PDB linked list

   Writes the ATOM records of a PDB file only writing those with the
   occupancy set to > 0.001 (used as a flag)

-  05.06.19 Original    By: ACRM
*/
void WriteFlaggedAtoms(FILE *out, PDB *pdb)
{
   PDB *p;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->occ > 0.001)
      {
         blWritePDBRecord(out, p);
      }
   }
}


/************************************************************************/
/*>void ListFlaggedResidues(FILE *out, PDB *pdb)
   ---------------------------------------------
*//**
   \param[in]    out   Output file pointer
   \param[in]    pdb   PDB linked list

   Writes a list of residues where the occupancy of the first atom in the
   residues is > 0.001 (used as a flag)

-  05.06.19 Original    By: ACRM
*/
void ListFlaggedResidues(FILE *out, PDB *pdb)
{
   PDB *res, *nextRes;
   for(res=pdb; res!=NULL; res=nextRes)
   {
      nextRes = blFindNextResidue(res);
      if(res->occ > 0.01)
      {
         char resid[SMALLBUFF];
         MAKERESID(resid, res);
         fprintf(out, "%s\n", resid);
      }
   }
}


/************************************************************************/
/*>BOOL ResInZone(PDB *res, ZONE *zones)
   -------------------------------------
*//**
   \param[in]    res     Pointer to the start of a residue
   \param[in]    zones   Linked list of zones
   \return               Is the residue in any of the zones

   Determine whether a residue is within any the zones specified in the
   linked list

-  05.06.19 Original    By: ACRM
*/
BOOL ResInZone(PDB *res, ZONE *zones)
{
   ZONE *z;
   
   for(z=zones; z!=NULL; NEXT(z))
   {
      if(blInPDBZoneSpec(res, z->start, z->stop))
      {
         return(TRUE);
      }
   }
   return(FALSE);
}


/************************************************************************/
/*>void ClearOccup(PDB *pdb)
   -------------------------
*//**
   \param[in]     pdb    PDB linked list

   Sets the occupancy to zero across the whole PDB file. Used as a flag

-  05.06.19 Original    By: ACRM
*/
void ClearOccup(PDB *pdb)
{
   PDB *p;
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->occ = 0.0;
   }
}


/************************************************************************/
/*>void FlagNearRes(PDB *pdb, ZONE *zones, char *restype, REAL radiusSq)
   ---------------------------------------------------------------------
*//**
   \param[in]    pdb       PDB linked list
   \param[in]    zones     Linked list of zone specifications
   \param[in]    restype   Amino acid type we are looking for
   \param[in]    radiusSq  Squared distance for a residue to be in range

   Uses occupancy as a flag - first clears this for all atoms, then
   looks for residues of the specified type (restype) within the given
   distance of any of the given zones. If a residue has any s/c atom in
   range then the occupancy for the residue is set to 1

-  05.06.19 Original    By: ACRM
*/
void FlagNearRes(PDB *pdb, ZONE *zones, char *restype, REAL radiusSq)
{
   PDB  *nextRes1, *nextRes2, *res1, *res2;

   /* Set occupancies to zero - we will use this as a flag              */
   ClearOccup(pdb);
   
   /* First step through the residues                                   */
   for(res1=pdb; res1!=NULL; res1=nextRes1)
   {
      nextRes1 = blFindNextResidue(res1);

      /* If this residue is in a zone                                   */
      if(ResInZone(res1, zones))
      {
         /* Second step through the residues                            */
         for(res2=pdb; res2!=NULL; res2=nextRes2)
         {
            nextRes2 = blFindNextResidue(res2);

            /* If this residue is of the correct type, but not in any
               zones
            */
            if(!strncmp(res2->resnam, restype, 3) &&
               !ResInZone(res2, zones))
            {
               /* See if any atoms in res1 are within distance of any
                  atoms in res2
               */
               if(CheckInRange(res1, nextRes1, res2, nextRes2, radiusSq))
               {
                  SetOccup(res2, nextRes2);
               }
            }
         }
      }
   }
}


/************************************************************************/
/*>void SetOccup(PDB *start, PDB *stop)
   ------------------------------------
*//**
   \param[in]   start   Start of a residue in PDB linked list
   \param[in]   stop    Start of next residue

   Sets the occupancy to 1 for all atoms from start to before stop

-  05.06.19 Original    By: ACRM
*/
void SetOccup(PDB *start, PDB *stop)
{
   PDB *p;
   for(p=start; p!=stop; NEXT(p))
   {
      p->occ = 1.0;
   }
}


/************************************************************************/
/*>BOOL CheckInRange(PDB *res1, PDB *nextRes1, PDB *res2, PDB *nextRes2, 
                     REAL radiusSq)
   ----------------------------------------------------------------------
*//**
   \param[in]    res1        Start of residue range (all atoms)
   \param[in]    nextRes1    End of residue range (all atoms)
   \param[in]    res2        Start of residue range (sidechain atoms)
   \param[in]    nextRes2    End of residue range (sidechain atoms)
   \param[in]    radiusSq    Distance squared
   \return                   BOOL: In range?

   Looks to see if sidechain atoms of res2 are within the specified
   distance of all atoms of res1. Returns true if they are, false 
   otherwise

-  05.06.19 Original    By: ACRM
*/
BOOL CheckInRange(PDB *res1, PDB *nextRes1, PDB *res2, PDB *nextRes2, 
                  REAL radiusSq)
{
   PDB *p, *q;
   
   for(p=res1; p!=nextRes1; NEXT(p))
   {
      for(q=res2; q!=nextRes2; NEXT(q))
      {
         if(ISSIDECHAIN(q))
         {
            if(DISTSQ(p,q) <= radiusSq)
            {
               return(TRUE);
            }
         }
      }
   }
   return(FALSE);
}


/************************************************************************/
/*>ZONE *ParseZoneSpec(char *zonespec)
   -----------------------------------
*//**
   \param[in]    zonespec   Zone specification
   \return                  Linked list of ZONEs

   Takes a comma separated list of zones each of which may be a single
   residue label or a range of labels (i.e. [c[.]]NNN[i]-[c[.]]NNN[i])
   and creates a linked list of ZONEs

-  05.06.19 Original    By: ACRM
*/
ZONE *ParseZoneSpec(char *zonespec)
{
   ZONE *zones = NULL, 
        *z     = NULL;
   char **zoneArray;
   int  i;
   
   if((zoneArray = blSplitStringOnCommas(zonespec, 8))==NULL)
      return(NULL);

   for(i=0; strlen(zoneArray[i]); i++)
   {
      if(zones == NULL)
      {
         INIT(zones, ZONE);
         z = zones;
      }
      else
      {
         ALLOCNEXT(z, ZONE);
      }
      
      if(z==NULL)
         return(NULL);
   
      PopulateZone(z, zoneArray[i]);
      free(zoneArray[i]);
   }
   free(zoneArray);

   return(zones);
}


/************************************************************************/
/*>void PopulateZone(ZONE *z, char *zoneDescription)
   -------------------------------------------------
*//**
   \param[in]    z                 Pointer to a ZONE structure
   \param[in]    zoneDescription   Zone specification

   Takes a zone specification which may be a single residue label or a 
   range of labels (i.e. [c[.]]NNN[i]-[c[.]]NNN[i]) and creates a ZONE
   structure

-  05.06.19 Original    By: ACRM
*/
void PopulateZone(ZONE *z, char *zoneDescription)
{
   char *dash, *z2;

   if((dash=strchr(zoneDescription, '-'))!=NULL)
   {
      z2=dash+1;
      *dash = '\0';
   }
   else
   {
      z2 = zoneDescription;
   }

   STRNCPY(z->start, zoneDescription, SMALLBUFF);
   STRNCPY(z->stop,  z2,              SMALLBUFF);
}


/************************************************************************/
/*>char **blSplitStringOnCharacter(char *string, char charac, 
                                   int minItemLen)
   -----------------------------------------------------------------------
*//**
   \param[in]   string       A string
   \param[in]   charac       A character
   \param[in]   minItemLen   Minimum length to allocate for each item
   \return                   Array of strings

   Takes a string and splits it into an array of strings at the specified
   character. This is done as a 2D array where each string has the same
   length.

-  05.06.19 Original    By: ACRM
*/
char **blSplitStringOnCharacter(char *string, char charac, int minItemLen)
{
   int  nitems = 0;
   char **items = NULL;
   char *c, 
        *buffer;
   int  maxItemLen = minItemLen-1,
        itemLen,
        i;

   /* Count the number of comma-separated items in the string. Also record
      the length of the longest item
   */
   itemLen = 0;
   for(c=string; *c; c++)
   {
      if(*c == charac)
      {
         if(itemLen > maxItemLen)
            maxItemLen = itemLen;
         nitems++;
         itemLen = 0;
      }
   }
   
   if(itemLen > maxItemLen)
      maxItemLen = itemLen;
   nitems++;
   maxItemLen++;

   /* Allocate space for the items                                      */
   if((items = (char **)blArray2D(sizeof(char), nitems+1, 
                                  maxItemLen))==NULL)
      return(NULL);
   
   /* And copy in the data                                              */
   buffer = string;
   for(i=0; i<nitems; i++)
   {
      if((c = strchr(buffer, charac))!=NULL)
         *c = '\0';
      strncpy(items[i], buffer, maxItemLen);
      buffer=c+1;
   }
   items[nitems][0] = '\0';

   return(items);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     REAL *radius, char *zonespec, char *restype,
                     char *listOnly)
   ---------------------------------------------------------------------
*//**
   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *radius      Neighbour radius
   \param[out]     *zonespec    Zone specification
   \param[out]     *restype     Res type to look for
   \param[out]     *listOnly    Only list the near residues rather than
                                PDB output
   \return                      Success?

   Parse the command line
   
-  05.06.19 Original    By: ACRM
-  18.06.19 Improved help message
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  REAL *radius, char *zonespec, char *restype, 
                  BOOL *listOnly)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = zonespec[0] = restype[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'r':
            argc--;
            argv++;
            sscanf(argv[0],"%lf",radius);
            break;
         case 'l':
            *listOnly = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2, 3 or 4 arguments left               */
         if((argc > 4) || (argc < 2))
            return(FALSE);
         
         /* Copy the first to zonespec                                  */
         STRNCPY(zonespec, argv[0], MAXBUFF);
         
         /* Copy the next to restype                                    */
         argc--;
         argv++;
         STRNCPY(restype, argv[0], MAXBUFF);
         UPPER(restype);
         
         /* If there's another, copy it to infile                       */
         argc--;
         argv++;
         if(argc)
         {
            STRNCPY(infile, argv[0], MAXBUFF);

            /* If there's another, copy it to outfile                   */
            argc--;
            argv++;
            if(argc)
            {
               STRNCPY(outfile, argv[0], MAXBUFF);
            }
         }
            
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
   Usage message

-  05.06.19 Original    By: ACRM
-  18.06.19 V1.0.1
-  19.06.19 V1.0.2

*/
void Usage(void)
{
   printf("\npdbfindneares V1.0.2 (c) 2019 UCL, Prof. Andrew C.R. \
Martin\n");

   printf("\nUsage: pdbfindnearres [-r nnn][-l] zone[,zone...] resnam \
[in.pdb [out.pdb]]\n");
   printf("       -r   Specify the radius used to look for nearby \
residues [%.3f]\n", DEFRAD);
   printf("       -l   Simply list residues instead of PDB output\n");

   printf("\nFinds occurrences of residues of type 'resnam' that have \
sidechains\n");
   printf("within the specified distance of any atoms in the specified \
residue\n");
   printf("range(s).\n");

   printf("\nzone is specified as a single residue specification or two \
residue\n");
   printf("        specifications separated by a dash (-).\n");

   printf("\n        A residue specification is of the form \
[c[.]]nnn[i]\n");
   printf("        where c   is an optional chain specification followed \
by a . if numeric\n");
   printf("              nnn is a residue number\n");
   printf("              i   is an optional insert code\n");
   printf("\n        Multiple zones may be listed separated by commas \
(,)\n");
   printf("\nresnam is a three-letter code amino acid name (upper or \
lower case)\n");
   printf("\nFor example:\n");
   printf("        pdbfindnearres L24-L34 tyr test.pdb\n");
   printf("        pdbfindnearres -l L50-L56 tyr test.pdb\n");
   printf("        pdbfindnearres -l L24-L34,L50-L56 tyr test.pdb\n");
   printf("        pdbfindnearres -l L24,L34,L50,L56 tyr test.pdb\n");
   printf("        pdbfindnearres -l -r 16 L50 lys test.pdb\n");
   printf("        pdbfindnearres -l -r 16 L50,L24-L34 lys test.pdb\n\n");
}

