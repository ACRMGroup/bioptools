/************************************************************************/
/**

   \file       checkforres.c
   
   \version    V1.3
   \date       22.07.14
   \brief      Checks whether a specified residue exists in a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin 2011-2014
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
-  V1.0   12.01.11  Original
-  V1.1   07.03.12  Now ignores HETATMs by default - reads them with
                    -H flag. Important for files like 1c26 which has 
                    waters are apparently valid residue numbers
-  V1.2   28.08.13  Modified for new ParseResSpec()
-  V1.3   22.07.14  Renamed deprecated functions with bl prefix.
                    Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "bioplib/pdb.h"
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
BOOL ParseCmdLine(int argc, char **argv, char *resid,
                  char *infile, char *outfile, BOOL *uppercaseresspec,
                  BOOL *readHet);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   main routine

-  22.07.96 Original    By: ACRM
-  29.09.05 Modified for -l By: TL
-  28.08.13 Modified for new ParseResSpec()
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char InFile[MAXBUFF],
        OutFile[MAXBUFF],
        resid[MAXBUFF],
        chain[8],
        insert[8];
   int  res,
        natom;
   PDB  *pdb;
   BOOL uppercaseresspec, 
        readHet;
   

   if(ParseCmdLine(argc, argv, resid, InFile, OutFile, 
                   &uppercaseresspec, &readHet))
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         BOOL ParseResSpecResult;

         if(readHet)
         {
            pdb=blReadPDB(in, &natom);
         }
         else
         {
            pdb=blReadPDBAtoms(in, &natom);
         }
         
         if(pdb==NULL)
         {
            fprintf(stderr,"checkforres: No atoms read from PDB file\n");
            return(1);
         }
         
         if (uppercaseresspec == TRUE) 
         {
            ParseResSpecResult = blParseResSpec(resid, chain, 
                                                &res, insert);
         }
         else 
         {
            ParseResSpecResult = blParseResSpecNoUpper(resid, chain, 
                                                       &res, insert);
         }
         if(!ParseResSpecResult)
         {
            fprintf(stderr,"checkforres: Illegal residue specification (%s)\n",
                    resid);
            return(1);
         }

         if((pdb = blFindResidue(pdb, chain, res, insert))!=NULL)
         {
            fprintf(out,"YES\n");
         }
         else
         {
            fprintf(out,"NO\n");
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
/*>BOOL ParseCmdLine(int argc, char **argv, char *resid,
                     char *infile, char *outfile, BOOL *uppercaseresspec)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *resid       Residue specifier
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *uppercaseresspec  Should residue spec be upcased?
                                (Default: yes)
   \param[out]     *readHet     Should we read hetatoms (Default: no)
   \param[out]     Success?

   Parse the command line
   
-  12.01.11 Original    By: ACRM
-  07.03.12 Added -H and *readHet
*/
BOOL ParseCmdLine(int argc, char **argv, char *resid,
                  char *infile, char *outfile, BOOL *uppercaseresspec,
                  BOOL *readHet)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   resid[0]               = '\0';
   *uppercaseresspec = TRUE;
   *readHet          = FALSE;

   if(!argc)               /* 05.11.07 Added this                       */
   {
      return(FALSE);
   }
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         case 'l':
            *uppercaseresspec = FALSE;
            break;
         case 'H':
            *readHet = TRUE;
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 1, 2 or 3 arguments left               */
         if(argc < 1 || argc > 3)
            return(FALSE);
         
         /* Copy the first to resid                                     */
         strcpy(resid, argv[0]);
         argc--;
         argv++;
         
         /* Copy the next to infile                                     */
         if(argc)
         {
            strcpy(infile, argv[0]);
            argc--;
            argv++;
         }
         
         /* If there's another, copy it to outfile                      */
         if(argc)
         {
            strcpy(outfile, argv[0]);
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
/*>void Usage(void)
   ----------------
*//**

   Prints a usage message

-  22.07.14 V1.3 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\ncheckforres V1.3 (c) 2011-2014, UCL, Dr. Andrew C.R. \
Martin\n");
   fprintf(stderr,"Usage: checkforres [-l][-H] resid [in.pdb [out.txt]]\n");
   fprintf(stderr,"       -l  Do not uppercase residue \
specifications.\n");
   fprintf(stderr,"           Default behaviour is to uppercase \
specifications.\n");
   fprintf(stderr,"       -H  Read HETATM records = i.e. allow residues \
that are HETATMs only\n");

   fprintf(stderr,"\nChecks whether a specified residue exists in a PDB \
file. The resid is in\n");
   fprintf(stderr,"the form [c]nnn[i] where [c] is an optional chain \
name, nnn is a residue\n");
   fprintf(stderr,"number and [i] is an optional insert code. By default \
the chain label is\n");
   fprintf(stderr,"up-cased unless the -l flag is used.\n\n");

}
