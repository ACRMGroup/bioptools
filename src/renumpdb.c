/************************************************************************/
/**

   \file       renumpdb.c
   
   \version    V1.9
   \date       22.07.14
   \brief      Renumber a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin / UCL 1994-6-2014
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
#define MAXBUFF 160
#define MAXCHAIN 64

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *DoSeq, BOOL *KeepChain, BOOL *DoAtoms,
                  char *chains, int *ResStart, int *AtomStart, 
                  BOOL *DoRes);
void Usage(void);
void DoRenumber(PDB *pdb, BOOL DoSequential, BOOL KeepChain, 
                BOOL DoAtoms, BOOL DoRes, char *chains, int *ResStart, 
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
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char infile[MAXBUFF],
        outfile[MAXBUFF],
        chains[MAXCHAIN];
   BOOL DoSequential = FALSE,
        KeepChain    = FALSE,
        DoAtoms      = TRUE,
        DoRes        = TRUE;
   PDB  *pdb;
   int  natoms,
        ResStart[MAXCHAIN],
        AtomStart;
        
   if(ParseCmdLine(argc, argv, infile, outfile, &DoSequential, &KeepChain,
                   &DoAtoms, chains, ResStart, &AtomStart, &DoRes))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb=blReadPDB(in, &natoms))==NULL)
         {
            fprintf(stderr,"No atoms read from input file\n");
         }
         else
         {
            DoRenumber(pdb,DoSequential,KeepChain,DoAtoms,DoRes,chains,
                       ResStart, AtomStart);
            blWritePDB(out, pdb);
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
                     char *chains, int *ResStart, int *AtomStart,
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
   \param[out]     *chains     Chain labels
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
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *DoSeq, BOOL *KeepChain, BOOL *DoAtoms, 
                  char *chains, int *ResStart, int *AtomStart, 
                  BOOL *DoRes)
{
   int  ChainCount;
   char *chp,
        *chq;

   for(ChainCount = 0; ChainCount < MAXCHAIN; ChainCount++)
      ResStart[ChainCount] = (-9999);

   *AtomStart = 1;
   ChainCount = 0;
   
   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   chains[0] = '\0';
   
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
            strncpy(chains,argv[0],MAXCHAIN);
            chains[MAXCHAIN-1] = '\0';
            UPPER(chains);
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
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   ResStart[ChainCount]  = (-9999);
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
*/
void Usage(void)
{
   fprintf(stderr,"\nRenumPDB V1.9 (c) 1994-2014 Dr. Andrew C.R. Martin, \
UCL\n");
   fprintf(stderr,"Freely distributable if no profit is made\n\n");
   fprintf(stderr,"Usage: renumpdb [-s] [-k] [-c <chains>] [-n] [-d] \
[-r <num>,<num>...] [-a <num> ]\n                \
[<infile>] [<outfile>]\n");
   fprintf(stderr,"       -s Renumber sequentially throughout \
structure\n");
   fprintf(stderr,"       -k Keep chain names when using -s\n");
   fprintf(stderr,"       -c Specify chain names to use\n");
   fprintf(stderr,"       -n Do not renumber atoms\n");
   fprintf(stderr,"       -d Do not renumber residues\n");
   fprintf(stderr,"       -r Specify resnum for start of each chain\n");
   fprintf(stderr,"       -a Specify first atom number\n\n");
   fprintf(stderr,"If files are not specified, stdin and stdout are \
used.\n");
   fprintf(stderr,"If a chain is to be skipped with -c or -r, use a - \
instead of the label or\nnumber.\n\n");
}
   

/************************************************************************/
/*>void DoRenumber(PDB *pdb, BOOL DoSequential, BOOL KeepChain, 
                   BOOL DoAtoms, BOOL DoRes, char *chains, int *ResStart,
                   int AtomStart)
   ------------------------------------------------------------
*//**

   \param[in,out]  *pdb            PDB linked list
   \param[in]      DoSequential    Number sequentially throughout
   \param[in]      KeepChain       Keep chain labels when DoSequential
   \param[in]      DoAtoms         Renumber atoms
   \param[in]      DoRes           Renumber residues
   \param[in]      *chains         Chain labels (or blank string)
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
*/
void DoRenumber(PDB *pdb, BOOL DoSequential, BOOL KeepChain, 
                BOOL DoAtoms, BOOL DoRes, char *chains, int *ResStart, 
                int AtomStart)
{
   PDB  *p;
   int  resnum   = 0,
        atnum    = AtomStart,
        ChainNum = 0,
        LastRes;
   char LastIns,
        LastChain,
        *ch;
   BOOL Incremented;
   
   ch = chains;
   
   LastRes     = (-1);
   LastIns     = ' ';
   LastChain   = '\0';
   
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
      if(p->chain[0] != LastChain)
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
Try -s option.\n",MAXCHAIN);
               exit(1);
            }
         }
         LastChain = p->chain[0];

         /* If it's not the first atom, step to the next chain name     */
         if(p != pdb && *ch)
            ch++;
      }
      
      /* Set the chain name if specified                                */
      if(*ch && *ch != '-') 
         p->chain[0] = *ch;
         
      /* If we're numbering sequentially and not keeping chain names,
         set the chain name to a blank
      */
      if(DoSequential && !KeepChain)
         p->chain[0] = ' ';
         
      if(DoRes) 
      {
         /* Set the residue number                                      */
         p->resnum = resnum;

         /* Set the insert code to a blank                              */
         p->insert[0] = ' ';
      }
      
      /* Set the atoms number                                           */
      if(DoAtoms)
         p->atnum = atnum++;
   }
}
