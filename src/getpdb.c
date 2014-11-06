/************************************************************************/
/**

   \file       GetPDB.c
   
   \version    V1.5
   \date       19.08.14
   \brief      Extract a numbered zone from a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin 1996-2014
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
-  V1.0   22.07.96  Original
-  V1.1   29.09.05  Added a command line option to prevent the residue
                    specifications being uppercased. This is required as 
                    some PDBs (eg 1rxsm) use upper and lower case letters 
                    as chain codes. (Tony Lewis) By: TL
-  V1.2   05.11.07  Improved error checking of command line
   V1.2.1 29.10.10  Fixed Bioplib bug where end of zone matched end of 
                    chain
-  V1.3   28.08.13  Modified for new ParseResSpec()
-  V1.4   22.07.14  Renamed deprecated functions with bl prefix.
                    Added doxygen annotation. By: CTP
-  V1.5   19.08.14  Fixed call to renamed function: 
                    blExtractZonePDBAsCopy() By: CTP
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
BOOL ParseCmdLine(int argc, char **argv, char *Zone1, char *Zone2,
                  char *infile, char *outfile, BOOL *uppercaseresspec);
void Usage(void);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   main routine

-  22.07.96 Original    By: ACRM
-  29.09.05 Modified for -l By: TL
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  19.08.14 Added AsCopy suffix for call to blExtractZonePDBAsCopy() 
            By: CTP
*/
int main(int argc, char **argv)
{
   FILE *in  = stdin,
        *out = stdout;
   char InFile[MAXBUFF],
        OutFile[MAXBUFF],
        Zone1[MAXBUFF],
        Zone2[MAXBUFF],
        chain1[8],  chain2[8],
        insert1[8], insert2[8];
   int  res1,    res2,
        natom;
   PDB  *pdb;
   BOOL uppercaseresspec;
   

   if(ParseCmdLine(argc, argv, Zone1, Zone2, InFile, OutFile, 
                   &uppercaseresspec))
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         BOOL ParseResSpec1Result, ParseResSpec2Result;

         if((pdb=blReadPDB(in, &natom))==NULL)
         {
            fprintf(stderr,"getpdb: No atoms read from PDB file\n");
            return(1);
         }
         
         if (uppercaseresspec == TRUE) 
         {
            ParseResSpec1Result = blParseResSpec(Zone1, chain1, 
                                                 &res1, insert1);
            ParseResSpec2Result = blParseResSpec(Zone2, chain2, 
                                                 &res2, insert2);
         }
         else 
         {
            ParseResSpec1Result = blParseResSpecNoUpper(Zone1, chain1, 
                                                        &res1, insert1);
            ParseResSpec2Result = blParseResSpecNoUpper(Zone2, chain2, 
                                                        &res2, insert2);
         }
         if(!ParseResSpec1Result)
         {
            fprintf(stderr,"getpdb: Illegal residue specification (%s)\n",
                    Zone1);
            return(1);
         }
         if(!ParseResSpec2Result)
         {
            fprintf(stderr,"getpdb: Illegal residue specification (%s)\n",
                    Zone2);
            return(1);
         }
         if((pdb = blExtractZonePDBAsCopy(pdb, chain1, res1, insert1,
                                          chain2, res2, insert2))==NULL)
         {
            fprintf(stderr,"getpdb: Zone not found (%s or %s)\n",
                    Zone1, Zone2);
            return(1);
         }
         
         blWritePDB(out,pdb);
      }
   }
   else
   {
      Usage();
   }
   return(0);
}

/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *Zone1, char *Zone2,
                     char *infile, char *outfile, BOOL *uppercaseresspec)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *Zone1       First end of zone
   \param[out]     *Zone2       Second end of zone
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *uppercaseresspec  Should residue spec be upcased?
                                (Default: yes)
   \return                     Success?

   Parse the command line
   
-  22.07.96 Original    By: ACRM
-  29.09.05 Added uppercaseresspec param and handling of -l  By: TL
-  05.11.07 Added first check that at least one parameter is on the
            command line
*/
BOOL ParseCmdLine(int argc, char **argv, char *Zone1, char *Zone2,
                     char *infile, char *outfile, BOOL *uppercaseresspec)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   Zone1[0]  = Zone2[0]   = '\0';
   *uppercaseresspec = TRUE;

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
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2, 3 or 4 arguments left               */
         if(argc < 2 || argc > 4)
            return(FALSE);
         
         /* Copy the first to Zone1 and second to Zone2                 */
         strcpy(Zone1, argv[0]);
         argc--;
         argv++;
         strcpy(Zone2, argv[0]);
         argc--;
         argv++;
         
         /* Copy the first to infile                                    */
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

-  22.07.14 V1.4 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\n");
   fprintf(stderr,"GetPDB V1.4 (c) 1996-2014, Dr. Andrew C.R. Martin, \
UCL.\n");
   fprintf(stderr,"            Modified by Tony Lewis, UCL, 2005\n");

   fprintf(stderr,"\nUsage: getpdb [-l] start end [in.pdb [out.pdb]]\n");
   fprintf(stderr,"       -l  Do not uppercase residue \
specifications.\n");
   fprintf(stderr,"           Default behaviour is to uppercase \
specifications.\n");
   fprintf(stderr,"\n");
   fprintf(stderr,"Start and end are residue specifications in the \
form [c]nnn[i]\n");
   fprintf(stderr,"where [c] is an optional chain specification, nnn is \
a residue number\n");
   fprintf(stderr,"and [i] is an optional insertion code.\n");

   fprintf(stderr,"\nGetPDB extracts a specified zone from a PDB file \
writing it out in\n");
   fprintf(stderr,"PDB format. I/O is through standard input/output if \
filenames are\n");
   fprintf(stderr,"not specified.\n\n");
}
