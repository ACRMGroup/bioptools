/*************************************************************************

   Program:    chaincontacts
   File:       chaincontacts.c
   
   Version:    V1.3
   Date:       28.10.15
   Function:   Calculate details of contacts between chains
   
   Copyright:  (c) Dr. Andrew C. R. Martin 1995-2015
   Author:     Dr. Andrew C. R. Martin
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk
               
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

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================
   V1.0  17.10.95 Original
   V1.1  07.04.06 Added ability to group chains
   V1.2  04.02.15 Now only reads ATOM records
   V1.3  28.10.15 Now takes a -H option to allow analysis of contacts 
                  with HETATOMs

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
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
#define DEF_RAD 3.0

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *radsq, char *chainsx, char *chainsy, BOOL *doHet,
                  BOOL *verbose, BOOL *stripWater);
void DoProteinProteinAnalysis(FILE *out, PDB *pdb, REAL RadSq, 
                              char *filename, char *chainsx, 
                              char *chainsy, BOOL verbose);
void DoProteinHetAnalysis(FILE *out, PDB *pdb, REAL RadSq, 
                          char *filename, char *chainsx, 
                          char *chainsy, BOOL verbose);
void PrintContacts(FILE *out, PDB *p, PDB *pe, PDB *q, PDB *qe, 
                   REAL RadSq, BOOL verbose);
BOOL InChainList(PDB *p, char *chains);
void PrintHeader(FILE *out, char *filename, REAL RadSq);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
   Main program for performing contact analysis

   17.10.95 Original    By: ACRM
   07.04.06 Added chainsx and chainsy checking
   04.03.15 Now just reads PDB atoms. Updated for new BiopLib
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        chainsx[MAXBUFF],
        chainsy[MAXBUFF];
   PDB  *pdb;
   int  natom;
   FILE *in        = stdin,
        *out       = stdout;
   REAL radsq      = DEF_RAD * DEF_RAD;
   BOOL doHet      = FALSE,
        verbose    = FALSE,
        keepWater  = FALSE;

   chainsx[0] = '\0';
   chainsy[0] = '\0';
   
   if(ParseCmdLine(argc, argv, infile, outfile, &radsq, chainsx, chainsy,
                   &doHet, &verbose, &keepWater))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if(doHet)
         {
            pdb = blReadPDB(in, &natom);
            if(!keepWater)
            {
               PDB *pdb2 = blStripWatersPDBAsCopy(pdb, &natom);
               FREELIST(pdb, PDB);
               pdb = pdb2;
            }
         }
         else
         {
            pdb = blReadPDBAtoms(in, &natom);
         }
         
         if(pdb != NULL)
         {
            if(doHet)
            {
               DoProteinHetAnalysis(out, pdb, radsq, infile, 
                                    chainsx, chainsy, verbose);
            }
            else
            {
               DoProteinProteinAnalysis(out, pdb, radsq, infile, 
                                        chainsx, chainsy, verbose);
            }
         }
         else
         {
            fprintf(stderr,"Warning: No atoms read from PDB file\n");
         }
      }
   }
   else
   {
      Usage();
      return(1);
   }
   
   return(0);
}
            
/************************************************************************/
/*>void Usage(void)
   ----------------
   Print usage message.

   17.10.95 Original    By: ACRM
   04.03.15 V1.2
   28.10.15 V1.3
*/
void Usage(void)
{
   fprintf(stderr,"\nChainContacts V1.3 (c) 1995-2015, Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"Usage: chaincontacts [-r radius] [-x CCC] \
[-y CCC] [-H [-w]] [in.pdb [out.dat]]\n");
   fprintf(stderr,"       -r Specify contact radius (Default: %.3f)\n\n",
           DEF_RAD);
   fprintf(stderr,"       -x/-y Specifiy one or more chains that form \
groups\n");
   fprintf(stderr,"       -H Group Y atoms are HETATOMs\n");
   fprintf(stderr,"       -w Include waters in Group Y HETATOMs\n");

   fprintf(stderr,"I/O is through stdin/stdout if files are not \
specified.\n\n");
   fprintf(stderr,"Performs a contact analysis at atom and residue \
level\n\n");
   fprintf(stderr,"If chains are specified for groups then only contacts \
between residues\n");
   fprintf(stderr,"in the -x and -y groups will be considered. So if you \
have an antibody\n");
   fprintf(stderr,"with chains L and H and antigen with chain C, then you \
can do -x LH -y C\n");
   fprintf(stderr,"to get only contacts between chain C with chain L or H. \
If you specify\n");
   fprintf(stderr,"just -x or -y then you will get contacts between that \
chain (or chains)\n");
   fprintf(stderr,"and every other chain.\n");
}
      

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     REAL *radsq, char *chainsx, char *chainsy, 
                     BOOL *doHet, BOOL *verbose, BOOL *keepWater)
   ---------------------------------------------------------------------
   Input:   int      argc        Argument count
            char     **argv      Argument array
   Output:  char     *infile     Input filename (or blank string)
            char     *outfile    Output filename (or blank string)
            REAL     *radsq      Contact radius squared
            BOOL     *doHet      Do contacts with HETATMs
            BOOL     *verbose    Print more information on contacts
            BOOL     *keepWater Strip waters when using -H
   Returns: BOOL                 Success

   Parse the command line

   17.10.95 Original    By: ACRM
   07.04.06 Added -x and -y
   28.10.15 Added -H and -v and -w
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  REAL *radsq, char *chainsx, char *chainsy,
                  BOOL *doHet, BOOL *verbose, BOOL *keepWater)
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
         case 'r':
            argv++;
            argc--;
            sscanf(argv[0],"%lf",radsq);
            (*radsq) *= (*radsq);
            break;
         case 'x':
            argv++;
            argc--;
            strncpy(chainsx, argv[0], MAXBUFF);
            break;
         case 'y':
            argv++;
            argc--;
            strncpy(chainsy, argv[0], MAXBUFF);
            break;
         case 'H':
            *doHet = TRUE;
            break;
         case 'v':
            *verbose = TRUE;
            break;
         case 'w':
            *keepWater = TRUE;
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
/*>void DoProteinProteinAnalysis(FILE *out, PDB *pdb, REAL RadSq, 
                                 char *filename, char *chainsx,
                                 char *chainsy, BOOL verbose)
   ----------------------------------------------------------------
   Input:   FILE    *out       Output file pointer
            PDB     *pdb       PDB linked list
            REAL    RadSq      Squared radius for contact
            char    *filename  Input PDB filename or blank string
            char    *chainsx   Group1 chains
            char    *chainsy   Group2 chains
            BOOL    verbose    Print more info on contacts

   Main routine to do the contacts analysis between protein chains

   17.10.95 Original (based on code from contacts.c)   By: ACRM
   07.04.06 Added chainsx and chainsy checking
   04.03.15 Updated for new BiopLib
   28.10.15 Renamed from DoAnalysis(). Refactored to take InChainList()
            and PrintContacts() into separate subroutines. Added verbose
*/   
void DoProteinProteinAnalysis(FILE *out, PDB *pdb, REAL RadSq, 
                              char *filename, char *chainsx, 
                              char *chainsy, BOOL verbose)
{
   PDB  *p,
        *pe,
        *q,
        *qe;
   BOOL ok1, ok2;

   PrintHeader(out, filename, RadSq);

   /* Run through the linked list a residue at a time                   */
   for(p=pdb; p!=NULL; p = pe)
   {
      pe = blFindNextResidue(p);
      
      /* Test if this residue is in a group X chain                     */
      ok1 = InChainList(p, chainsx);
      
      /* It is in a group X chain                                       */
      if(ok1)
      {
         /* Run through the PDB linked list again a residue at a time   */
         for(q=pdb; q!=NULL; q = qe)
         {
            qe = blFindNextResidue(q);
         
            /* Check it's a different chain                             */
            if(p->chain[0]  != q->chain[0])
            {
               /* If it is, check if it's a group Y chain               */
               ok2 = InChainList(q, chainsy);
            
               /* If it is, then run through the first residue          */
               if(ok2)
               {
                  PrintContacts(out, p, pe, q, qe, RadSq, verbose);
               }
            }
         }
      }
   }
}

/************************************************************************/
/*>void PrintContacts(FILE *out, PDB *p, PDB *pe, PDB *q, PDB *qe, 
                      REAL RadSq, BOOL verbose)
   ---------------------------------------------------------------
   Input:   FILE    *out       Output file pointer
            PDB     *pdb       PDB linked list
            REAL    RadSq      Squared radius for contact
            char    *filename  Input PDB filename or blank string
            char    *chainsx   Group1 chains
            char    *chainsy   Group2 chains
            BOOL    doHet      Do contacts with HET groups
            BOOL    verbose    Print more information
*/
void PrintContacts(FILE *out, PDB *p, PDB *pe, PDB *q, PDB *qe, 
                   REAL RadSq, BOOL verbose)
{
   int  NContacts = 0;
   PDB  *p1, *q1;
   REAL DistSq;

   for(p1 = p; p1 != pe; NEXT(p1))
   {
      /* Run through the second residue                                 */
      for(q1 = q; q1 != qe; NEXT(q1))
      {
         /* Calculate distance and display if in range                  */
         DistSq = DISTSQ(p1,q1);
         if(DistSq <= RadSq)
         {
            NContacts++;
         }
      }
   }
   /* Print information for residue contacts                            */
   if(NContacts)
   {
      if(verbose)
      {
         fprintf(out,"Chain: %c Res:%4d%c %4s - \
Chain: %c Res:%4d%c %4s Contacts: %2d %s\n",
                 p->chain[0], p->resnum, p->insert[0], p->resnam,
                 q->chain[0], q->resnum, q->insert[0], q->resnam,
                 NContacts,
                 !strncmp(q->record_type, "HETATM", 6)?"(HET)":"");
      }
      else
      {
         fprintf(out,"Chain: %c Res:%4d%c - \
Chain: %c Res:%4d%c Contacts: %2d %s\n",
                 p->chain[0], p->resnum, p->insert[0],
                 q->chain[0], q->resnum, q->insert[0],
                 NContacts,
                 !strncmp(q->record_type, "HETATM", 6)?"(HET)":"");
      }
   }
}

/************************************************************************/
/*>BOOL InChainList(PDB *p, char *chains)
   --------------------------------------
   Input:   PDB  *p       PDB entry
            char *chains  Array of chain labels
   Returns: BOOL          Checks if the specified residue is in the list
                          of chains. The list of chains is blank, it also
                          counts as being in the list

   28.10.15 Refactored from DoProteinProteinAnalysis()
*/
BOOL InChainList(PDB *p, char *chains)
{
   BOOL ok = TRUE;
   char *chp;
   
   if(chains[0])
   {
      ok = FALSE;
      for(chp = chains; *chp; chp++)
      {
         if(p->chain[0] == *chp)
         {
            ok = TRUE;
            break;
         }
      }
   }
   return(ok);
}


/************************************************************************/
/*>void PrintHeader(FILE *out, char *filename, REAL RadSq)
   -------------------------------------------------------
   Input:  FILE *out        Output file
           char *filename   Input filename
           REAL RadSq       Radius squared for contacts

   Prints the output header

   28.10.15 Refactored out of DoProteinProteinAnalysis()
*/
void PrintHeader(FILE *out, char *filename, REAL RadSq)
{
   /* Print atom level header                                           */
   fprintf(out,"Contact Analysis\n");
   fprintf(out,"================\n\n");

   fprintf(out,"File:   %s\n",(filename[0] ? filename : "stdin"));
   fprintf(out,"Radius: %f\n",sqrt(RadSq));

   /* Now do residue level analysis                                     */
   fprintf(out,"Residue level contacts\n");
   fprintf(out,"----------------------\n\n");
   
}


/************************************************************************/
/*>void DoProteinHetAnalysis(FILE *out, PDB *pdb, REAL RadSq, 
                             char *filename, char *chainsx, 
                             char *chainsy, BOOL verbose)
   ----------------------------------------------------------------
   Input:   FILE    *out       Output file pointer
            PDB     *pdb       PDB linked list
            REAL    RadSq      Squared radius for contact
            char    *filename  Input PDB filename or blank string
            char    *chainsx   Group1 chains
            char    *chainsy   Group2 chains
            BOOL    verbose    Print more information on contacts

   Main routine to do the contacts analysis between protein chains
   and HET groups

   28.10.15 Original
*/   
void DoProteinHetAnalysis(FILE *out, PDB *pdb, REAL RadSq, 
                          char *filename, char *chainsx, 
                          char *chainsy, BOOL verbose)
{
   PDB  *p,
        *pe,
        *q,
        *qe;
   BOOL ok1, ok2;

   PrintHeader(out, filename, RadSq);

   /* Run through the linked list a residue at a time                   */
   for(p=pdb; p!=NULL; p=pe)
   {
      pe = blFindNextResidue(p);

      /* If this is a protein atom                                      */
      if(!strncmp(p->record_type, "ATOM  ", 6))
      {
         /* Test if this residue is in a group X chain                  */
         ok1 = InChainList(p, chainsx);
         
         /* It is in a group X chain                                    */
         if(ok1)
         {
            /* Run through the PDB linked list again a residue at a time*/
            for(q=pdb; q!=NULL; q = qe)
            {
               qe = blFindNextResidue(q);

               /* Check it's a HETATM group                             */
               if(!strncmp(q->record_type, "HETATM", 6))
               {
                  /* If it is, check if it's a group Y chain            */
                  ok2 = InChainList(q, chainsy);
                  
                  /* If it is, then run through the first residue       */
                  if(ok2)
                  {
                     PrintContacts(out, p, pe, q, qe, RadSq, verbose);
                  }
               }
            }
         }
      }
   }
}

