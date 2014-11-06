/*************************************************************************

   Program:    pdbsphere
   File:       pdbsphere.c
   
   Version:    V1.8
   Date:       27.07.12
   Function:   Output all aminoacids within range from central aminoacid 
               in a PDB file
   
   Copyright:  (c) UCL/Anja Baresic/Dr. Andrew C.R. Martin
   Author:     Anja Baresic/Dr. Andrew C.R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
 
   Email:      anya@biochem.ucl.ac.uk; andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   Inputs ID of a central residue and a PDB file. 
   Outputs all residues from that PDB file, within 8 angstroms (default 
   range) of central residue's coordinates.

**************************************************************************

   Usage:
   ======
   
**************************************************************************

   Revision History:
   =================
   V1.0  23.10.07  Original
   V1.1  05.11.07  Added -r, -s, -h command line options
   V1.2  17.12.07  Summary output format changed to 
                   [chain]:resnum:[insert]
                   By: Anja
   V1.5  17.05.11  Default summary output format changed to 
                   [chain].resnum[insert] (as used by ProFit et al)
                   Use -c to get colon separated format. Cleaned up
                   usage message
                   Various code cleanup
                   Changed to use PDB->extras rather than PDB->occ
   V1.6  26.10.11  Modified to take -H flag which causes the program to
                   assume the residue spec is for a HETATM instead of an
                   ATOM
   V1.7  14.05.12  Added -a option
   V1.8  27.07.12  Fixed bug in handling (lack of) -a

**************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "bioplib/pdb.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/macros.h"

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
void FlagResiduesInRange(PDB *pdb, PDB *central, REAL radiusSq);
void WriteAtoms(PDB *pdb, FILE *out);
void WriteResidues(PDB *pdb, FILE *out, BOOL colons, BOOL compact);
BOOL ParseCmdLine(int argc, char **argv,char *resspec, char *InFile, 
                  char *OutFile, BOOL *summary, REAL *radiusSq,
                  BOOL *colons, BOOL *isHet, BOOL *doAuto);
void ClearExtras(PDB *pdb);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Takes a PDB file and a central residue ID in [chain[.]]num[insert] 
   format and writes a PDB file containing only those residues within
   a specified radius (default 8A, override with -r). Summary output
   (just the residue list) can be generated with -s and -c provides an
   alternative output format.
*/
int main(int argc, char **argv)
{
   FILE   *in  = stdin,
          *out = stdout;
   PDB    *pdb,
          *central;         
   int    natom;
   REAL   radiusSq;
   char   resspec[MAXBUFF],
          InFile[MAXBUFF],
          OutFile[MAXBUFF];
   BOOL   summary,
          colons = FALSE,
          isHet = FALSE,
          doAuto = FALSE;
   
   if (ParseCmdLine(argc, argv, resspec, InFile, OutFile, &summary, 
                    &radiusSq, &colons, &isHet, &doAuto))
   {
      if (OpenStdFiles(InFile, OutFile, &in, &out))
      {
         if((pdb=ReadPDB(in, &natom))==NULL)
         {
            fprintf(stderr,"Error: (pdbsphere) No atoms read from PDB \
file\n");
            return(1);
         }
         else
         {        
            /* Clear the ->extras field                                 */
            ClearExtras(pdb);

            if(doAuto)
            {
               for(central=pdb; central!=NULL; central=FindNextResidue(central))
               {
                  FlagResiduesInRange(pdb, central, radiusSq);
                  fprintf(out, "%s %s%d%s:",central->resnam, 
                          central->chain,
                          central->resnum,
                          central->insert);
                  WriteResidues(pdb, out, colons, TRUE);
                  ClearExtras(pdb);
               }
            }
            else
            {
               if(isHet)
               {
                  central=FindHetatmResidueSpec(pdb, resspec);
               }
               else
               {
                  central=FindResidueSpec(pdb, resspec);
               }
               
               if(central==NULL)
               {
                  fprintf(stderr,"Error: (pdbsphere) Residue %s not \
found in %s\n", resspec, InFile);
                  return(1);
               }
               else                  
               {
                  FlagResiduesInRange(pdb, central, radiusSq);
               
                  if (summary)
                  {
                     WriteResidues(pdb, out, colons, FALSE);
                  }
                  else
                  {
                     WriteAtoms(pdb, out);
                  }      
               }
            }
         }
      }      
   }
   else 
   {
      Usage();
   }
   return 0;      
}

/**********************************************************************/
/*>void FlagResiduesInRange(PDB *pdb, PDB *central, REAL *radiusSq)
   ----------------------------------------------------------------
   Input:   PDB    *pdb     
            PDB    *central    Pointer to the first atom of a central 
                               residue
            REAL   radiusSq    To be flagged, atom has to be within that 
                               range (radius is squared for speed)

   If any atom in a residue is within range from *central, marks all 
   atoms in that residue (sets extras field to 1) 

   17.05.11 Changed double to REAL and use extras field rather than occ
            By: ACRM
*/
void FlagResiduesInRange(PDB *pdb, PDB *central, REAL radiusSq)
{
   PDB  *p,
        *q,
        *current,
        *nextPRes,
        *nextCurrentRes;
   BOOL aaInRange;
      

   nextPRes=FindNextResidue(central);
  
   for (current=pdb; current!=NULL; current=nextCurrentRes)
   {
      aaInRange=FALSE;
      nextCurrentRes=FindNextResidue(current);
      
      for (q=current; q!=nextCurrentRes; NEXT(q))
      {            
         for (p=central; p!=nextPRes; NEXT(p))
         {
            if (DISTSQ(p, q)<radiusSq)
            {
               aaInRange=TRUE;
               /*q=nextCurrentRes;
                 break;*/
            }
         }
      }
      
      if (aaInRange)
      {
         for (q=current; q!=nextCurrentRes; NEXT(q))
         {
            q->extras=(APTR)1;
         }    
      }
   }
} 


/************************************************************************/
/*>void WriteAtoms(PDB *pdb, FILE *out)
   ------------------------------------
   Input:    PDB  *pdb   pointer to beginning of a linked list
             FILE *out   output file

   Writes atom's node in *out if extras is set

   17.05.11 Uses PDB->extras instead of PDB->occ
*/
void WriteAtoms(PDB *pdb, FILE *out)
{  
   PDB *p;

   for (p=pdb; p!=NULL; NEXT(p))
   {
      if (p->extras)
      {
         WritePDBRecord(out, p);
      }
   }
}


/************************************************************************/
/*>void WriteResidues(PDB *pdb, FILE *out, BOOL colons, BOOL compact)
   ------------------------------------------------------------------
   Input:    PDB  *pdb     pointer to beginning of a linked list
             FILE *out     output file
             BOOL *colons  Include colons in output format

   Writes a list of residues' IDs for residues marked with 
   occ>1.90 (returns occ to original value)

   17.05.11 New default output format and takes colon parameter to get
            old format and uses extras rather than occ
   14.05.12 Added compact
*/
void WriteResidues(PDB *pdb, FILE *out, BOOL colons, BOOL compact)
{
   PDB *p;
   
   for (p=pdb; p!=NULL; p=FindNextResidue(p))
   {
      if (p->extras)
      {
         p->extras=(APTR)0;
         if(compact)
         {
            if(isdigit(p->chain[0]))
            {
               fprintf(out, " %c.%d%c", p->chain[0], p->resnum, 
                       p->insert[0]);
            }
            else
            {
               fprintf(out, " %c%d%c", p->chain[0], p->resnum, 
                       p->insert[0]);
            }
         }
         else
         {
            if(colons)
            {
               fprintf(out, "%s:%d:%s\n", p->chain, p->resnum, p->insert);
            }
            else
            {
               if(isdigit(p->chain[0]))
               {
                  fprintf(out, "%c.%d%c\n", p->chain[0], p->resnum, 
                          p->insert[0]);
               }
               else
               {
                  fprintf(out, "%c%d%c\n", p->chain[0], p->resnum, 
                          p->insert[0]);
               }
            }
         }
      }
   }

   if(compact)
   {
      fprintf(out,"\n");
   }
   
}  


/***********************************************************************/
/*>void Usage(void)
   ----------------
   Print Usage message

   17.05.11 V1.5 By: ACRM
   26.10.11 V1.6 By: ACRM
   14.05.12 V1.7 By: ACRM
   27.07.12 V1.8 By: ACRM
*/

void Usage(void)
{
   fprintf(stderr,"\n");
   fprintf(stderr,"PDBsphere V1.8 (c) 2011-2012 UCL, Anja Baresic, \
Andrew Martin.\n");
   fprintf(stderr,"\nUsage: \
pdbsphere [-s] [-c] [-r radius] [-h] [-H] resID\n                 \
[in.pdb [out.pdb/out.txt]]\n");
   fprintf(stderr,"-or-   \
pdbsphere -a [-r radius] [in.pdb [out.txt]]\n");
   fprintf(stderr,"       -s  Output summary: only list of residue \
IDs.\n");
   fprintf(stderr,"       -c  Colon separated summary format.\n");
   fprintf(stderr,"       -H  Residue spec is for a HETATM\n");
   
   fprintf(stderr,"       -r  Set your own allowed range to radius.\n");
   fprintf(stderr,"       -a  'Auto' mode - analyses all residues \
producing\n           summary for each.\n");

   fprintf(stderr,"\npdbsphere identifies residues within a specified \
radius of a specified\n");
   fprintf(stderr,"residue. All atoms of any residue containing at least \
one atom within \n");
   fprintf(stderr,"range (default 8A) are output. The -s option provides \
a summary format\n");
   fprintf(stderr,"listing the residues in range instead of providing \
PDB output.\n");

   fprintf(stderr,"\nResID is in form [c[.]]num[i]where [c] is an \
optional chain specification\n");
   fprintf(stderr,"with an optional '.' for numeric chain IDs, num is \
the residue number and\n");
   fprintf(stderr,"[i] is an optional insertion code.\n");
   fprintf(stderr,"I/O is through standard input/output if files not \
specified.\n\n");   
   
}

/**********************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *resspec, char *InFile, 
                    char *OutFile, BOOL *summary, REAL *radiusSq,
                    BOOL *colons, BOOL *isHet, BOOL *doAuto)
   ----------------------------------------------------------------------
   Input:   int    argc         Argument count
            char   **argv       Argument array
   Output:  char   *resspec     Central residue ID in [chain]num[insert]
                                format
            char   *InFile      Input file (or blank string)
            char   *OutFile     Output file (or blank string)
            BOOL   *summary     Should output be summarised?
                                (Default: no)
            REAL *radiusSq      Maximum allowed distance of residues' 
                                coordinates - squared
                                (Default:64, max range:8 angstroms)
            BOOL   *colons      Colon separated output format
            BOOL   *isHet       Residue spec is for a HETATM (-H)
   Returns: BOOL                Success?

   Parse the command line
   
   26.10.07 Original    By: Anya
   17.05.11 Changed double to REAL; Added *colons     By: ACRM
   26.10.11 Added -H
   14.05.12 Added -a
   27.07.12 Fixed bug in checking of doAuto
*/
BOOL ParseCmdLine(int argc, char **argv,char *resspec, char *InFile, 
                  char *OutFile, BOOL *summary, REAL *radiusSq, 
                  BOOL *colons, BOOL *isHet, BOOL *doAuto)
{
   argc--;
   argv++;

   InFile[0] = OutFile[0] = '\0';
   *summary = FALSE;
   *radiusSq = 64.00;
   *colons   = FALSE;
   *isHet    = FALSE;
   *doAuto   = FALSE;
 
   if (argc<1)
   {      
      return (FALSE);
   }
   
   while(argc)
   {     
      if(argv[0][0] == '-')
      {
         if(argv[0][2]!='\0')
         {
            return(FALSE);
         }
         else
         {
            switch(argv[0][1])
            {
            case 'h':
               return(FALSE);
               break;
            case 's':
               *summary = TRUE;
               break;
            case 'r':
               argc--;
               argv++;
               if(argc < 0)
                  return(FALSE);
               if(!sscanf(argv[0],"%lf",radiusSq))
                  return(FALSE);
               else
                  *radiusSq *= *radiusSq;  
               break;
            case 'c':
               *colons = TRUE;
               break;
            case 'H':
               *isHet = TRUE;
               break;
            case 'a':
               *doAuto = TRUE;
               break;
            default:
               return(FALSE);
               break;
            }
         }
      }
      else
      {
         /* Check that there are 1, 2 or 3 arguments left               */
         if(argc<1 || argc > 3)
            return(FALSE);

         /* Copy the first to resspec                                   */
         if(!(*doAuto))
         {
            if(argc)
            {
               strncpy(resspec, argv[0], MAXBUFF);
               argc--;
               argv++;
            }
            else
            {
               return(FALSE);
            }
         }

         /* Copy the second to InFile                                   */
         if(argc)
         {
            strncpy(InFile, argv[0], MAXBUFF);
            argc--;
            argv++;
         }         

         /* If there's another, copy it to OutFile                      */
         if(argc)
         {
            strncpy(OutFile, argv[0], MAXBUFF);            
            argc--;
            argv++;
         }
            
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>void ClearExtras(PDB *pdb)
   --------------------------
   Clear the extras flag in the PDB linked list

   17.05.11 Original   By: ACRM
*/
void ClearExtras(PDB *pdb)
{
   PDB *p;

   for(p=pdb; p!=NULL; NEXT(p))
   {
      p->extras = (APTR)0;
   }
}
