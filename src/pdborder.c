/************************************************************************/
/**

   \file       pdborder.c
   
   \version    V1.8
   \date       13.03.19
   \brief      Correct the atom order in a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL 1994-2019
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
-  V1.0  24.08.94 Original
-  V1.1  15.01.97 Tidied a few bits
-  V1.2  31.05.02 Changed PDB field from 'junk' to 'record_type'
-  V1.3  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.4  07.11.14 Initialized variables
-  V1.5  13.02.15 Added whole PDB support and fixed some core dumps
-  V1.6  05.03.15 Replaced blFindEndPDB() with blFindNextResidue()
-  V1.7  12.03.15 Changed to allow multi-character chain names
-  V1.8  13.03.19 Fixed possible unterminated string

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
#include "bioplib/array.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define MAXATOMS 19
#define MAXATNAM 8
#define TERM_MIDCHAIN 0
#define TERM_NTER     1
#define TERM_CTER     2

/************************************************************************/
/* Globals
*/
char ***gAtomLists = NULL;
BOOL gVerbose = FALSE,
     gWarnH   = FALSE;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *COLast, BOOL *GromosILE);
void Usage(void);
PDB *CorrectOrder(PDB *pdb, BOOL COLast);
PDB *CorrectResidue(PDB *start, PDB *end, BOOL COLast, int terminus);
void FixGromosILE(PDB *pdb,BOOL Gromos);
void SpliceNTerHs(PDB **from, PDB **to);
void SpliceCTerOs(PDB **from, PDB **to);
BOOL GotSomeHydrogens(PDB *pdb);
BOOL InitializeGlobalAtomLists(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for format conversion

-  23.08.94 Original    By: ACRM
-  24.08.94 Changed to call OpenStdFiles()
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  13.02.15 Added whole PDB support and initialize atom lists
            dynamically  By: ACRM
*/
int main(int argc, char **argv)
{
   FILE     *in  = stdin,
            *out = stdout;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];
   WHOLEPDB *wpdb;
   PDB      *pdb;
   BOOL     COLast    = FALSE,
            GromosILE = FALSE;
   
   if(ParseCmdLine(argc, argv, infile, outfile, &COLast, &GromosILE))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb = blReadWholePDB(in)) != NULL)
         {
            pdb = wpdb->pdb;

            if(!InitializeGlobalAtomLists())
            {
               fprintf(stderr,"pdborder (error): No memory for atom \
lists\n");
               return(1);
            }

            FixGromosILE(pdb,GromosILE);
            
            if((pdb = CorrectOrder(pdb, COLast)) != NULL)
            {
               if(gWarnH)
               {
                  if(GotSomeHydrogens(pdb))
                  {
                     fprintf(stderr,"pdborder (warning): There were \
hydrogens missing.\n");
                  }
               }
               blWriteWholePDB(out, wpdb);
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
      Usage();
   }
   
   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *COLast, BOOL *GromosILE)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *COLast      Put the C and O after the s/c
   \param[out]     *GromosILE   Change ILE CD1 to CD for GROMOS
   \return                      Success?

   Parse the command line

   -v sets global gVerbose
   
-  23.08.94 Original    By: ACRM
-  24.08.94 Added verbose option
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *COLast, BOOL *GromosILE)
{
   argc--;
   argv++;

   infile[0]  = outfile[0] = '\0';
   *COLast    = FALSE;
   *GromosILE = FALSE;
   gVerbose   = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'g':
            *COLast    = TRUE;
            *GromosILE = TRUE;
            break;
         case 'c':
            *COLast = TRUE;
            break;
         case 'i':
            *GromosILE = TRUE;
            break;
         case 'v':
            gVerbose = TRUE;
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
-  24.08.94 Added -v
-  15.01.97 V1.1
-  22.07.14 V1.2 By: CTP
-  07.11.14 V1.4 By: ACRM
-  13.02.15 V1.5 
-  05.03.15 V1.6
-  12.03.15 V1.7
-  13.03.19 V1.8
*/
void Usage(void)
{
   fprintf(stderr,"\npdborder V1.8 (c) 1994-2019, Andrew C.R. Martin, \
UCL\n\n");
   fprintf(stderr,"Usage: pdborder [-c] [-i] [-g] [in.pdb \
[out.pdb]]\n");
   fprintf(stderr,"       -c = N CA s/c C O order\n");
   fprintf(stderr,"       -i = ILE has CD instead of CD1\n");
   fprintf(stderr,"       -g = -c -i (i.e. for Gromos)\n");
   fprintf(stderr,"       -v = Report all missing Hs\n\n");
   fprintf(stderr,"Correct atom order of a PDB file.\n");
   fprintf(stderr,"By default, fixes ILE CD to CD1 and imposes standard \
N,CA,C,O,s/c atom\n");
   fprintf(stderr,"ordering.\n\n");
}

/************************************************************************/
/*>PDB *CorrectOrder(PDB *pdb, BOOL COLast)
   ----------------------------------------
*//**

   \param[in]      *pdb     PDB linked list
   \param[in]      COLast   Put the C,O after the s/c
   \return                  Fixed order linked list or NULL if error

-  23.08.94 Original    By: ACRM
-  24.08.94 Added code to handle additional Hs at NTER and O at CTER
            Only tries to re-order ATOM records.
-  07.11.14 Initialized variables
-  05.03.15 Replaced blFindEndPDB() with blFindNextResidue()
-  12.03.15 Changed to allow multi-character chain names
*/
PDB *CorrectOrder(PDB *pdb, BOOL COLast)
{
   PDB  *start = NULL,
        *prev  = NULL,
        *ret   = NULL,
        *end   = NULL,
        *p;
   char chain[8];
   int  terminus;
   BOOL GotNTER = FALSE;
   
   strcpy(chain, "-");
   start = pdb;
   ret   = NULL;
   
   while(start != NULL)
   {
      end = blFindNextResidue(start);
      
      /* Only correct the order if this is an amino acid residue        */
      if(!strncmp(start->record_type,"ATOM  ",6))
      {
         /* Find out where we are in the chain                          */
         terminus = TERM_MIDCHAIN;

         if(!CHAINMATCH(start->chain, chain))
         {
            terminus = TERM_NTER;
            strncpy(chain, start->chain, 8);
         }
         else if(end==NULL                      || 
                 !CHAINMATCH(end->chain, chain) || 
                 !strncmp(end->resnam, "CTER",4))
         {
            terminus = TERM_CTER;
         }

         /* Set flag so residue after NTER will be treated as Nter      */
         if(!strncmp(start->resnam,"NTER",4))
            GotNTER = TRUE;
         if(GotNTER && strncmp(start->resnam,"NTER",4))
         {
            GotNTER = FALSE;
            terminus = TERM_NTER;
         }
         
         /* Correct the order of this residue and add to linked list    */
         if((p = CorrectResidue(start, end, COLast, terminus)) != NULL) 
         {
            if(ret == NULL)
            {
               prev = p;
               ret = p;
            }
            else
            {
               prev->next = p;
            }
         }
         else
         {
            return(NULL);
         }
      }
         
      if(prev != NULL)            /* 13.02.15                           */
      {
         /* Walk to the last atom in this residue                       
            13.02.15 Added check on p!=NULL
          */
         for(p = prev->next; 
             p != NULL && p->next != NULL && p->next != end; 
             NEXT(p)) ;
         prev = p;
      }
      
      /* Step to the next residue                                       */
      start = end;
   }

   return(ret);
}
   
/************************************************************************/
/*>PDB *CorrectResidue(PDB *start, PDB *end, BOOL COLast, int terminus)
   --------------------------------------------------------------------
*//**

   \param[in]      *start      Start of residue linked list
   \param[in]      *end        Start of next residue (or NULL)
   \param[in]      COLast      Put the carboxyl at the end
   \param[in]      terminus    Flag to indicate position in chain
                               (TERM_MIDCHAIN, TERM_NTER, or TERM_CTER)
   \return                     Fixed order PDB linked list

   Do the actual order correction for this residue

-  23.08.94 Original    By: ACRM
-  24.08.94 Added terminus handling code
-  12.03.15 Changed to allow multi-character chain names
*/
PDB *CorrectResidue(PDB *start, PDB *end, BOOL COLast, int terminus)
{
   int  offset,
        i;
   PDB  *p,
        *ret    = NULL;
   BOOL special = FALSE;

   if(!strncmp(start->resnam,"NTER",4) || 
      !strncmp(start->resnam,"CTER",4))
      special = TRUE;
   
   /* Terminate the PDB linked list at the end of this residue          */
   for(p=start; p->next != end; NEXT(p)) ;
   p->next = NULL;

   /* Search the atom order list for this residue type                  */
   for(offset = 0; gAtomLists[offset][0] != NULL; offset++)
   {
      if(!strncmp(p->resnam,gAtomLists[offset][0],4))
         break;
   }

   /* Issue warning if residue type not found                           */
   if(gAtomLists[offset][0] == NULL)
   {
      fprintf(stderr,"Warning: Residue type `%s' unknown. Atom order \
unchanged.\n",p->resnam);
      return(start);
   }

   /* Find each atom in order and build into the ret linked list        */
   for(i=1; gAtomLists[offset][i] != NULL; i++)
   {
      /* If not CO last or not C or not O then move it                  */
      if(!COLast || 
         (strncmp(gAtomLists[offset][i],"C   ",4) &&
          strncmp(gAtomLists[offset][i],"O   ",4)))
      {
         /* Find this PDB record                                        */
         for(p=start; p!=NULL; NEXT(p))
         {
            if(!strncmp(p->atnam,gAtomLists[offset][i],4))
               break;
         }
         
         if(p==NULL)
         {
            if((gAtomLists[offset][i][0] != 'H') ||
               (gVerbose && 
                (terminus != TERM_NTER || 
                 strncmp(gAtomLists[offset][i],"H   ",4))))
            {
               if(terminus != TERM_NTER || 
                  strncmp(gAtomLists[offset][i],"N   ",4))
               {
                  if(terminus != TERM_CTER ||
                     strncmp(gAtomLists[offset][i],"O   ",4))
                  {
                     if(start != NULL)
                     {
                        fprintf(stderr,"Warning: Missing atom `%s' in \
residue %s %s.%d%c\n",          gAtomLists[offset][i],
                                start->resnam,
                                start->chain,
                                start->resnum,
                                start->insert[0]);
                     }
                     else if(ret != NULL)
                     {
                        fprintf(stderr,"Warning: Missing atom `%s' in \
residue %s %s.%d%c\n",          gAtomLists[offset][i],
                                ret->resnam,
                                ret->chain,
                                ret->resnum,
                                ret->insert[0]);
                     }
                  }
               }
            }
            else if(!gVerbose && gAtomLists[offset][i][0] == 'H' &&
                    (terminus != TERM_NTER ||
                     strncmp(gAtomLists[offset][i],"H   ",4)))
            {
               gWarnH = TRUE;
            }
         }
         else
         {
            blMovePDB(p, &start, &ret);
         }
      }
   }

   /* If the C and O are last, move them now                            */
   if(COLast && !special)
   {
      for(i=0; i<2; i++)
      {
         char atom[8];
         
         if(i==0)
            strcpy(atom,"C   ");
         else
            strcpy(atom,"O   ");
         
         /* Find this PDB record                                        */
         for(p=start; p!=NULL; NEXT(p))
         {
            if(!strncmp(p->atnam,atom,4))
               break;
         }
         
         if(p==NULL)
         {
            if(terminus != TERM_CTER ||
               strncmp(atom,"O   ",4))
            {
               if(start != NULL)
               {
                  fprintf(stderr,"Warning: Missing atom `%s' in residue \
%s %s.%d%c\n",            atom,
                          start->resnam,
                          start->chain,
                          start->resnum,
                          start->insert[0]);
               }
               else if(ret != NULL)
               {
                  fprintf(stderr,"Warning: Missing atom `%s' in residue \
%s %s.%d%c\n",            atom,
                          ret->resnam,
                          ret->chain,
                          ret->resnum,
                          ret->insert[0]);
               }
            }               
         }
         else
         {
            blMovePDB(p, &start, &ret);
         }
      }
   }

   
   if(terminus == TERM_NTER)       /* Handle extra hydrogens at NTER    */
      SpliceNTerHs(&start, &ret);
   else if(terminus == TERM_CTER)  /* Handle extra oxygen at CTER       */
      SpliceCTerOs(&start, &ret);
   
   /* Check to see if any records haven't been moved                    */
   for(p=start; p!=NULL; NEXT(p))
   {
      fprintf(stderr,"Warning: Extra atom `%s' in residue %s %s.%d%c\n",
              p->atnam,
              p->resnam,
              p->chain,
              p->resnum,
              p->insert[0]);
   }
   
   return(ret);
}
   
/************************************************************************/
/*>void FixGromosILE(PDB *pdb, BOOL Gromos)
   ----------------------------------------
*//**

   \param[in,out]  *pdb      PDB linked list
   \param[in]      Gromos    True: Change to CD
                             False: Change to CD1

   Fix the PDB linked list to have CD rather than CD1 for ILE and
   fix the atom lists table.

-  23.08.94 Original    By: ACRM
-  15.01.97 Added the Gromos parameter (previously always -> Gromos form)
-  13.03.19 Fixed possible non-terminated string
*/
void FixGromosILE(PDB *pdb, BOOL Gromos)
{
   PDB *p;
   int offset, i;
   
   if(Gromos)
   {
      /* Fix PDB linked list ILE CD1->CD                                */
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(!strncmp(p->resnam,"ILE ",4) &&
            !strncmp(p->atnam, "CD1 ",4))
         {
            strcpy(p->atnam,     "CD  ");
            strcpy(p->atnam_raw, " CD ");
         }
      }
      
      /* Fix the atom order list ILE CD1->CD                            */
      for(offset = 0; gAtomLists[offset][0] != NULL; offset++)
      {
         if(!strncmp("ILE ",gAtomLists[offset][0],4))
         {
            for(i=0; gAtomLists[offset][i] != NULL; i++)
            {
               if(!strncmp(gAtomLists[offset][i],"CD1 ",4))
               {
                  strncpy(gAtomLists[offset][i],"CD  ", 5);
                  break;
               }
            }
            break;
         }
      }
   }
   else
   {
      /* Fix PDB linked list ILE CD->CD1                                */
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if(!strncmp(p->resnam,"ILE ",4) &&
            !strncmp(p->atnam, "CD  ",4))
            strcpy(p->atnam,"CD1 ");
      }
   }
}

/************************************************************************/
/*>void SpliceNTerHs(PDB **from, PDB **to)
   ---------------------------------------
*//**

   Moves Nterminal hydrogens to the output PDB linked list

-  24.08.94 Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
void SpliceNTerHs(PDB **from, PDB **to)
{
   PDB *ret    = NULL,
       *p      = NULL,
       *H1     = NULL,
       *H2     = NULL,
       *H3     = NULL,
       *NT     = NULL,
       *N      = NULL;

   if(*from == NULL)
      return;
   
   /* Search the from list for the 3 hydrogens and NT                   */
   for(p = (*from); p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"H1  ",4) || !strncmp(p->atnam,"HT1 ",4))
         H1 = p;
      if(!strncmp(p->atnam,"H2  ",4) || !strncmp(p->atnam,"HT2 ",4))
         H2 = p;
      if(!strncmp(p->atnam,"H3  ",4) || !strncmp(p->atnam,"HT3 ",4))
         H3 = p;
      if(!strncmp(p->atnam,"NT  ",4))
         NT = p;
   }

   /* Search the to list for N                                          */
   for(p = (*to); p!= NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"N   ",4))
         N = p;
   }

   /* Move H1 and H2 to the ret list                                    */
   blMovePDB(H1, from, &ret);
   blMovePDB(H2, from, &ret);

   /* Try to move N and if this fails, move NT                          */
   if(!blMovePDB(N, to, &ret))
      blMovePDB(NT, from, &ret);
   
   /* Move H3 to the ret list                                           */
   blMovePDB(H3, from, &ret);

   if(ret != NULL)
   {
      /* Now link the remains of the to list to ret                     */
      p = ret;
      LAST(p);
      p->next = *to;
      
      /* Now reset the to list to be the ret list                       */
      *to = ret;
   }
}


/************************************************************************/
/*>void SpliceCTerOs(PDB **from, PDB **to)
   ---------------------------------------
*//**

   Moves Cterminal oxygens to the output PDB linked list

-  24.08.94 Original   By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
void SpliceCTerOs(PDB **from, PDB **to)
{
   PDB  *p    = NULL,
        *O1   = NULL,
        *O2   = NULL,
        *O    = NULL,
        *OXT  = NULL,
        *C    = NULL,
        *tail = NULL;

   /* Search the from list for the 3 oxygens                            */
   for(p = (*from); p!=NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"O1  ",4) || 
         !strncmp(p->atnam,"OT1 ",4))
         O1 = p;
      if(!strncmp(p->atnam,"O2  ",4) || 
         !strncmp(p->atnam,"OT2 ",4))
         O2 = p;
      if(!strncmp(p->atnam,"OXT ",4))
         OXT = p;
   }

   /* Search the to list for O and C                                    */
   for(p = (*to); p!= NULL; NEXT(p))
   {
      if(!strncmp(p->atnam,"O   ",4))
         O = p;
      if(!strncmp(p->atnam,"C   ",4))
         C = p;
   }

   /* Split the to list after the C                                     */
   if(C != NULL)
   {
      tail = C->next;
      C->next = NULL;
   }
   
   /* Try to move O from tail to the to list, if we fail, move O1.
      The first move will be quite happy even if tail hasn't been
      defined.
   */
   if(!blMovePDB(O, &tail, to))
      blMovePDB(O1, from, to);
   
   /* Now move O2 from the from list to the to list                     */
   blMovePDB(O2, from, to);

   /* Now join tail onto the to list                                    */
   p = *to;
   LAST(p);
   p->next = tail;

   /* Last, try to move OXT to the to list                              */
   blMovePDB(OXT, from, to);
}


/************************************************************************/
/*>BOOL GotSomeHydrogens(PDB *pdb)
   -------------------------------
*//**

   \param[in]      *pdb     PDB linked list to search for hydrogens
   \return                  Were any hydrogens found?

   Tests whether there are any hydrogens in the PDB linked list.

-  15.01.97 Original   By: ACRM
*/
BOOL GotSomeHydrogens(PDB *pdb)
{
   PDB *p;
   for(p=pdb; p!=NULL; NEXT(p))
   {
      if(p->atnam[0] == 'H')
         return(TRUE);
   }
   return(FALSE);
}   

/************************************************************************/
/*>BOOL InitializeGlobalAtomLists(void)
   ------------------------------------
*//**

   \return                    Success?

   Initializes the global atom name array. We can't just use the static
   hard-coded version since we need to write to it for Gromos ILE CD

-  13.02.15  Original   By: ACRM
*/
BOOL InitializeGlobalAtomLists(void)
{
   char *atomLists[][MAXATOMS] =   /* Residue/atom table                  */
      {  {"ALA ","N   ","H   ","CA  ","C   ","O   ","CB  ",NULL  ,NULL  ,
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"CYS ","N   ","H   ","CA  ","C   ","O   ","CB  ","SG  ",NULL  ,
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"CYS1","N   ","H   ","CA  ","C   ","O   ","CB  ","SG  ",NULL  ,
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"CYS2","N   ","H   ","CA  ","C   ","O   ","CB  ","SG  ",NULL  ,
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"CYSH","N   ","H   ","CA  ","C   ","O   ","CB  ","SG  ","HG  ",
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"ASP ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","OD1 ",
          "OD2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"GLU ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ",
          "OE1 ","OE2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"PHE ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ",
          "CD2 ","CE1 ","CE2 ","CZ  ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"GLY ","N   ","H   ","CA  ","C   ","O   ",NULL  ,NULL  ,NULL  ,
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"HIS ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","ND1 ",
          "HD1 ","CD2 ","CE1 ","NE2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"HIS1","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","ND1 ",
          "HD1 ","CD2 ","CE1 ","NE2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"HISB","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","ND1 ",
          "CD2 ","CE1 ","NE2 ","HE2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"HISH","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","ND1 ",
          "HD1 ","CD2 ","CE1 ","NE2 ","HE2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"ILE ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG1 ","CG2 ",
          "CD1 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"LYS ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ",
          "CE  ","NZ  ","HZ1 ","HZ2 ","HZ3 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"LEU ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ",
          "CD2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"MET ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","SD  ",
          "CE  ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"ASN ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","OD1 ",
          "ND2 ","HD21","HD22",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"PRO ","N   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ",NULL  ,
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"GLN ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ",
          "OE1 ","NE2 ","HE21","HE22",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"ARG ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","CD  ",
          "NE  ","HE  ","CZ  ","NH1 ","HH11","HH12","NH2 ","HH21","HH22",NULL},
         {"SER ","N   ","H   ","CA  ","C   ","O   ","CB  ","OG  ","HG  ",
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"THR ","N   ","H   ","CA  ","C   ","O   ","CB  ","OG1 ","HG1 ",
          "CG2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"VAL ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG1 ","CG2 ",
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"TRP ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ",
          "CD2 ","NE1 ","HE1 ","CE2 ","CE3 ","CZ2 ","CZ3 ","CH2 ",NULL  ,NULL},
         {"TYR ","N   ","H   ","CA  ","C   ","O   ","CB  ","CG  ","CD1 ",
          "CD2 ","CE1 ","CE2 ","CZ  ","OH  ","HH  ",NULL  ,NULL  ,NULL  ,NULL},
         {"CTER","OT2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {"NTER","HT1 ","HT2 ",NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL},
         {NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,
          NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL  ,NULL}
      }  ;
   int nres, atom, res;

   /* Count how many residues we have specified                         */           
   for(nres=0; atomLists[nres][0] != NULL; nres++);
   nres++;

   /* Creat the malloc'd 3D array                                       */
   if((gAtomLists = (char ***)blArray3D(sizeof(char), nres, 
                                        MAXATOMS, MAXATNAM))==NULL)
   {
      return(FALSE);
   }

   /* Copy in the names                                                 */
   for(res=0; res<nres; res++)
   {
      for(atom=0; atom<MAXATOMS; atom++)
      {
         /* Free up the memory if there was no atom name                */
         if(atomLists[res][atom] == NULL)
         {
            free(gAtomLists[res][atom]);
            gAtomLists[res][atom] = NULL;
         }
         else
         {
            strcpy(gAtomLists[res][atom], atomLists[res][atom]);
         }
      }
   }

   return(TRUE);
}
