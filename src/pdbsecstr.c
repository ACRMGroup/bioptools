/************************************************************************/
/**

   \File       pdbsecstr.c
   
   \version    V1.2
   \date       06.08.18
   \brief      Secondary structure calculation program
   
   \copyright  (c) Dr. Andrew C. R. Martin, UCL, 1999-2018
   \author     Dr. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College London,
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
   V1.0   19.05.99 Original, written while at Inpharmatica   By: ACRM
   V1.1   11.08.16 Rewritten to use PDB files rather than XMAS files
                   and to use blCalcSecStrucPDB() in Bioplib
   V1.2   06.08.18 Updated Usage message

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bioplib/pdb.h"
#include "bioplib/general.h"
#include "bioplib/macros.h"
#include "bioplib/secstr.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 160
   
/************************************************************************/
/* Globals
*/
BOOL sNewChain = FALSE;

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *debug);
void Usage(void);
void WriteResults(FILE *out, PDB *pdbStart, PDB *pdbStop);



/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for secondary structure calculation.

-  19.05.99 Original   By: ACRM
-  27.05.99 Added error return if blCalcSS out of memory
-  11.08.16 Updated for using Bioplib
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   FILE *in = stdin,
        *out = stdout;
   PDB  *pdb;
   int  natoms;
   BOOL debug = FALSE;
   
   
   if(!ParseCmdLine(argc, argv, infile, outfile, &debug))
   {
      Usage();
      return(0);
   }
   else
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if((pdb = blReadPDBAtoms(in, &natoms))!=NULL)
         {
            PDB *start, *stop;
            
            for(start=pdb; start!=NULL; start=stop)
            {
               stop=blFindNextChain(start);

               if(blCalcSecStrucPDB(start, stop, debug) != 0)
               {
                  return(1);
               }
            
               WriteResults(out, start, stop);
            }
            
            FREELIST(pdb, PDB);

            if(in  != stdin)  fclose(in);
            if(out != stdout) fclose(out);
         }
      }
      else
      {
         return(1);
      }
   }
         
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *debug)
   ---------------------------------------------------------------------
*//**
   \param[in]   argc              Argument count
   \param[in]   **argv            Argument array
   \param[out]  *infile           Input filename (or blank string)
   \param[out]  *outfile          Output filename (or blank string)
   \param[out]  *debug            Debug?
   \return                        Success

   Parse the command line

-   19.05.99 Original    By: ACRM
-   11.08.16 Updated for PDB version
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *debug)
{
   argc--;
   argv++;
   
   infile[0] = outfile[0] = '\0';
   *debug    = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd':
            *debug = TRUE;
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

   Print a usage message

-  19.05.99 Original   By: ACRM
-  21.05.99 Added flags
-  11.08.16 Updated for non-xmas version
-  06.08.18 V1.2
*/
void Usage(void)
{
   fprintf(stderr,"\npdbsecstr V1.2 (c) 1999-2018, UCL, \
Dr. Andrew C.R. Martin\n");

   fprintf(stderr,"\nUsage: pdbsecstr [-d] [in.pdb [out.pdb]]\n");
   fprintf(stderr,"          -d Debug mode - reports information on\
dropped 3rd Hbonds, etc.\n");

   fprintf(stderr,"\nCalculates secondary structure assignments \
according to the method of\n");
   fprintf(stderr,"Kabsch and Sander. Reads a PDB file and writes \
a simple summary text\n");
   fprintf(stderr,"file.\n");
   fprintf(stderr,"\nInput/output is to standard input/output if \
files are not specified.\n\n");
}


/************************************************************************/
/*>void WriteResults(FILE *out, PDB *pdbStart, PDB *pdbStop)
   ---------------------------------------------------------
*//**
   \param[in]   *out        Output file pointer
   \param[in]   *pdbStart   Start of PDB linked list
   \param[in]   *pdbStop    End of PDB linked list (points to the item
                            after the stop point)

   Write a summary file with the residue names and secondary structure

   21.05.99 Original   By: ACRM
   11.08.16 Changed to use blBuildResSpec()
*/
void WriteResults(FILE *out, PDB *pdbStart, PDB *pdbStop)
{
   PDB *p;
   
   for(p=pdbStart; p!=pdbStop; p=blFindNextResidue(p))
   {
      char resspec[24];
      blBuildResSpec(p, resspec);

      fprintf(out, "%-6s %s %c\n",
              resspec,
              p->resnam,
              p->secstr);
   }
}

