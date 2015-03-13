/************************************************************************/
/**

   \file       pdbrenum.c
   
   \version    V2.0
   \date       10.03.15
   \brief      Renumber a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin / UCL 1994-2015
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
-  V1.0  29.06.94 Original
-  V1.1  04.07.94 Added -r and -a and option for skipping in -c
-  V1.2  16.08.94 Wasn't setting insert to blank (!)
-  V1.3  24.08.94 Changed to call OpenStdFiles()
-  V1.4  14.12.94 Added check on chain count
-  V1.5  04.01.95 Added -d option
-  V1.6  06.10.95 Was blanking insert codes with -d
-  V1.7  13.05.96 Fixed potential segmentation error if command line
                  specified wrongly
-  V1.8  14.11.96 -r wasn't working with blank chain names
-  V1.9  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.10 06.11.14 Renamed from renumpdb  By: ACRM
-  V1.11 13.02.15 Added whole PDB support
-  V1.12 23.02.15 Modified for new blWritePDBTrailer
                  Uses blRenumberAtomsPDB() to do the atoms
-  V1.13 02.03.15 Deals better with header and trailer
-  V2.0  10.03.15 Chains specified with -c are now comma separated

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/pdb.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF  160
#define MAXCHAIN 160
#define MAXCHAINLABEL 8

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *DoSeq, BOOL *KeepChain, BOOL *DoAtoms,
                  char ***chains, int *ResStart, int *AtomStart, 
                  BOOL *DoRes);
void Usage(void);
void DoRenumber(PDB *pdb, BOOL DoSequential, BOOL KeepChain, 
                BOOL DoAtoms, BOOL DoRes, char **chains, int *ResStart, 
                int AtomStart);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for PDB renumbering

-  29.06.94 Original    By: ACRM
-  04.07.94 Added ResStart & AtomStart
-  24.08.94 Changed to call OpenStdFiles()
-  04.01.95 Added DoRes handling
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  13.02.15 Added whole PDB support   By: ACRM
-  23.02.15 Modified for new blWriteWholePDBTrailer()
-  02.03.15 Now uses blWriteWholePDBHeaderNoRes() if residues have been
            renumbered and always does blWriteWholePDBTrailer() since
            this now deals properly with renumbered atoms.
-  10.03.15 Chains now an array of strings
*/
int main(int argc, char **argv)
{
   FILE     *in  = stdin,
            *out = stdout;
   char     infile[MAXBUFF],
            outfile[MAXBUFF],
            **chains;
   BOOL     DoSequential = FALSE,
            KeepChain    = FALSE,
            DoAtoms      = TRUE,
            DoRes        = TRUE;
   WHOLEPDB *wpdb;
   PDB      *pdb;
   int      ResStart[MAXCHAIN],
            AtomStart;
        
   if(ParseCmdLine(argc, argv, infile, outfile, &DoSequential, &KeepChain,
                   &DoAtoms, &chains, ResStart, &AtomStart, &DoRes))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((wpdb=blReadWholePDB(in))==NULL)
         {
            fprintf(stderr,"pdbrenum: Unable to read input PDB file\n");
         }
         else
         {
            int numTer;
            
            pdb=wpdb->pdb;
            DoRenumber(pdb,DoSequential,KeepChain,DoAtoms,DoRes,chains,
                       ResStart, AtomStart);
            if(DoRes)
            {
               blWriteWholePDBHeaderNoRes(out, wpdb);
            }
            else
            {
               blWriteWholePDBHeader(out, wpdb);
            }
            numTer = blWritePDB(out, pdb);
            blWriteWholePDBTrailer(out, wpdb, numTer);
         }
      }
      else
      {
         Usage();
         return(1);
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     BOOL *DoSeq, BOOL *KeepChain, BOOL *DoAtoms, 
                     char ***chains, int *ResStart, int *AtomStart,
                     BOOL *DoRes)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \param[out]     *outfile    Output filename (or blank string)
   \param[out]     *DoSeq      Number residues sequentially throughout
   \param[out]     *KeepChain  Keep chain labels when DoSeq TRUE
   \param[out]     *DoAtoms    Renumber atoms
   \param[out]     ***chains     Chain labels
   \param[out]     *ResStart   Chain residue start numbers
   \param[out]     *AtomStart  First chain atom start number
   \param[out]     *DoRes      Renumber residues
   \return                     Success

   Parse the command line

-  29.06.94 Original    By: ACRM
-  04.07.94 Added Residue and atom start specification
            Initialise chains array.
-  04.01.95 Added -d
-  13.05.96 Checks argc after -c/-r/-a options
-  13.02.15 Added -x
-  10.03.15 -c now takes comma separated chain list
            Chains now an array of strings
            Removed redundant ForceHeader
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *DoSeq, BOOL *KeepChain, BOOL *DoAtoms, 
                  char ***chains, int *ResStart, int *AtomStart, 
                  BOOL *DoRes)
{
   int  ChainCount;
   char *chp,
        *chq;
   BOOL oldStyle;


   oldStyle = blCheckProgName(argv[0], "renumpdb");
   *chains = NULL;

   for(ChainCount = 0; ChainCount < MAXCHAIN; ChainCount++)
      ResStart[ChainCount] = (-9999);

   *AtomStart = 1;
   ChainCount = 0;
   
   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 's':
            *DoSeq = TRUE;
            break;
         case 'c':
            if(!(--argc))
            {
               Usage();
               exit(1);
            }
            argv++;
            if(oldStyle)
            {
               if((*chains = blSplitStringOnChars(argv[0]))
                  ==NULL)
               {
                  fprintf(stderr,"No memory for storing chain labels: %s\n",
                          argv[0]);
                  exit(1);
               }
            }
            else
            {
               if((*chains = blSplitStringOnCommas(argv[0], MAXCHAINLABEL))
                  ==NULL)
               {
                  fprintf(stderr,"No memory for storing chain labels: %s\n",
                          argv[0]);
                  exit(1);
               }
            }
            
            break;
         case 'k':
            *KeepChain = TRUE;
            break;
         case 'n':
            *DoAtoms = FALSE;
            break;
         case 'd':
            *DoRes = FALSE;
            break;
         case 'f':
            fprintf(stderr,"-f is now deprecated\n");
            break;
         case 'r':
            if(!(--argc))
            {
               Usage();
               exit(1);
            }
            argv++;

            for(chp=chq=argv[0]; chq!=NULL; )
            {
               if((chq = strchr(chp,','))!=NULL)
                  *chq = '\0';

               if(!strcmp(chp,"-"))
               {
                  ResStart[ChainCount++] = (-9999);
               }
               else
               {
                  if(!sscanf(chp,"%d",&(ResStart[ChainCount++])))
                     return(FALSE);
               }

               if(chq!=NULL)
                  chp = chq+1;
            }

            break;
         case 'a':
            if(!(--argc))
            {
               Usage();
               exit(1);
            }
            argv++;
            if(!sscanf(argv[0],"%d",AtomStart))
               return(FALSE);
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
            
         ResStart[ChainCount]  = (-9999);
         if(*chains==NULL)
         {
            if((*chains = blSplitStringOnCommas("", MAXCHAINLABEL))
               ==NULL)
            {
               fprintf(stderr,"No memory for storing chain labels: %s\n",
                       argv[0]);
               exit(1);
            }
         }
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   ResStart[ChainCount]  = (-9999);
   if(*chains == NULL)
   {
      if((*chains = blSplitStringOnCommas("", MAXCHAINLABEL))
         ==NULL)
      {
         fprintf(stderr,"No memory for storing chain labels: %s\n",
                 argv[0]);
         exit(1);
      }
   }
   return(TRUE);
}

/************************************************************************/
/*>void Usage(void)
   ----------------
*//**

   Print a usage message

-  29.06.94 Original    By: ACRM
-  04.07.94 V1.1
-  16.08.94 V1.2  
-  24.08.94 V1.3
-  04.01.95 V1.5
-  06.10.95 V1.6  
-  13.05.96 V1.7
-  14.11.96 V1.8
-  22.07.14 V1.9 By: CTP
-  06.11.14 V1.10 By: ACRM
-  13.02.15 V1.11
-  23.02.15 V1.12
-  02.03.15 V1.13
-  10.03.15 V2.0
*/
void Usage(void)
{
   fprintf(stderr,"\npdbrenum V2.0 (c) 1994-2015 Dr. Andrew C.R. \
Martin, UCL\n");
   fprintf(stderr,"Usage: pdbrenum [-s][-k][-c chain[,chain[...]]]\
[-n][-d]\n");
   fprintf(stderr,"                [-r num[,num][...]]][-a num]\
[in.pdb [out.pdb]]\n");
   fprintf(stderr,"       -s Renumber sequentially throughout \
structure\n");
   fprintf(stderr,"       -k Keep chain names when using -s\n");
   fprintf(stderr,"       -c Specify chain names to use\n");
   fprintf(stderr,"       -n Do not renumber atoms\n");
   fprintf(stderr,"       -d Do not renumber residues\n");
   fprintf(stderr,"       -r Specify resnum for start of each chain\n");
   fprintf(stderr,"       -a Specify first atom number\n\n");

   fprintf(stderr,"\nRenumbers the residues and atoms of a PDB file \
allowing start residues\n");
   fprintf(stderr,"and chain labels to be specified and sequential \
numbering throughout\n");
   fprintf(stderr,"multiple chains.\n");

   fprintf(stderr,"If files are not specified, stdin and stdout are \
used.\n");
   fprintf(stderr,"If a chain is to be skipped with -c or -r, use a - \
instead of the label or\nnumber.\n\n");

   fprintf(stderr,"If called as 'renumpdb' instead of 'pdbrenum', the \
old behaviour with\n");
   fprintf(stderr,"-c is used of only allowing 1-letter chain labels \
with no separating\n");
   fprintf(stderr,"comma. e.g. Chains L and H would be specified as LH \
instead of L,H\n");
}
   

/************************************************************************/
/*>void DoRenumber(PDB *pdb, BOOL DoSequential, BOOL KeepChain, 
                   BOOL DoAtoms, BOOL DoRes, char **chains, int *ResStart,
                   int AtomStart)
   -----------------------------------------------------------------------
*//**

   \param[in,out]  *pdb            PDB linked list
   \param[in]      DoSequential    Number sequentially throughout
   \param[in]      KeepChain       Keep chain labels when DoSequential
   \param[in]      DoAtoms         Renumber atoms
   \param[in]      DoRes           Renumber residues
   \param[in]      **chains        Chain labels
   \param[in]      *ResStart       Residue start numbers
   \param[in]      AtomStart       Atom number for start of first chain

   Do the actual renumbering.

-  29.06.94 Original    By: ACRM
-  04.07.94 Modified to handle ResStart and AtomStart
-  16.08.94 Sets the insert to a blank
-  14.12.94 Added check on chain count
-  04.01.95 Added DoRes parameter
-  06.10.95 If dores was false, inserts were getting blanked.
-  14.11.96 Initialise LastChain to something other than ' ' as
            -r option wasn't working with blank chain names
-  23.02.15 Now uses blRenumberAtomsPDB() for atom numbering
-  10.03.15 chains now an array of strings
*/
void DoRenumber(PDB *pdb, BOOL DoSequential, BOOL KeepChain, 
                BOOL DoAtoms, BOOL DoRes, char **chains, int *ResStart, 
                int AtomStart)
{
   PDB  *p;
   int  resnum   = 0,
        ChainNum = 0,
        ChainIndex = 0,
        LastRes;
   char LastIns,
        LastChain[MAXCHAINLABEL];
   BOOL Incremented;

   LastRes      = (-1);
   LastIns      = ' ';
   LastChain[0] = '\0';
   
   for(p=pdb; p!=NULL; NEXT(p))
   {
      Incremented = FALSE;
      
      /* Increment resnum if we have changed residue                    */
      if(p->resnum    != LastRes || 
         p->insert[0] != LastIns)
      {
         LastRes = p->resnum;
         LastIns = p->insert[0];
         
         resnum++;
         Incremented = TRUE;
      }

      /* See if we've changed chain                                     */
      if(!CHAINMATCH(p->chain, LastChain))
      {
         if(DoSequential)
         {
            if(!Incremented)
            {
               resnum++;
               Incremented = TRUE;
            }
         }
         else
         {
            if(ResStart[ChainNum] == (-9999))
               resnum = 1;
            else
               resnum = ResStart[ChainNum];
            if(++ChainNum >= MAXCHAIN)
            {
               fprintf(stderr,"Maximum number of chains (%d) exceeded. \
Try -s option or increase MAXCHAIN.\n",MAXCHAIN);
               exit(1);
            }
         }
         strcpy(LastChain, p->chain);

         /* If it's not the first atom, step to the next chain name     */
         if((p != pdb) && chains[ChainIndex][0])
            ChainIndex++;
      }
      
      /* Set the chain name if specified                                */
      if(chains[ChainIndex][0] && (chains[ChainIndex][0] != '-'))
         strcpy(p->chain, chains[ChainIndex]);
         
      /* If we're numbering sequentially and not keeping chain names,
         set the chain name to 'A'
      */
      if(DoSequential && !KeepChain)
         strcpy(p->chain, "A");
         
      if(DoRes) 
      {
         /* Set the residue number                                      */
         p->resnum = resnum;

         /* Set the insert code to a blank                              */
         p->insert[0] = ' ';
      }
   }
   /* Renumber atoms if required                                        */
   if(DoAtoms)
      blRenumAtomsPDB(pdb, AtomStart);
}
