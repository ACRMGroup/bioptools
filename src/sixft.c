/************************************************************************/
/**

   Program:    sixft
   \file       sexft.c
   
   \version    V1.0
   \date       10.11.17   
   \brief      Six-frame translation
   
   \copyright  (c) UCL / Dr. Andrew C. R. Martin 2017
   \author     Dr. Andrew C. R. Martin
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

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include "bioplib/macros.h"
#include "bioplib/general.h"
#include "bioplib/sequtil.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 256

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL DoTranslate(FILE *in, FILE *out, BOOL showDNA, BOOL showRF);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *showDNA, BOOL *showRF);
void Usage(void);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program

-  10.11.17 Original   By: ACRM
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   BOOL showDNA = FALSE,
        showRF  = FALSE;
   FILE *in     = stdin,
        *out    = stdout;

   if(ParseCmdLine(argc, argv, infile, outfile, &showDNA, &showRF))
   {
      if(showDNA && showRF)
         showDNA = FALSE;

      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if(!DoTranslate(in, out, showDNA, showRF))
         {
            fprintf(stderr, "Error: No memory for translation\n");
            return(1);
         }
      }
      else
      {
         fprintf(stderr, "Error: Unable to open input or output file\n");
         return(1);
      }
   }
   else
   {
      Usage();
   }
   
   return(0);
}

/************************************************************************/
/*>BOOL DoTranslate(FILE *in, FILE *out, BOOL showDNA, BOOL showRF)
   ----------------------------------------------------------------
*//**
   \param[in]   in       Input file pointer
   \param[in]   out      Output file pointer
   \param[in]   showDNA  Display the DNA sequence
   \param[in]   showRF   Display the reading frame sequence
   \return               Success?

   Does the actual work of performing the 6FT and finding the best
   translation

-  10.11.17 Original   By: ACRM
*/
BOOL DoTranslate(FILE *in, FILE *out, BOOL showDNA, BOOL showRF)
{
   char header[MAXBUFF];
   char *dna = NULL;
   BOOL returnStatus = TRUE;
   char *orf         = NULL,
        *protein     = NULL;
   
   while((dna = blReadFASTA(in, header, MAXBUFF))!=NULL)
   {
      LOWER(dna);
      
      if((orf = (char *)malloc((strlen(dna)+1) * sizeof(char))) == NULL)
      {
         returnStatus = FALSE;
         break;
      }
      
      if((protein = blSixFTBest(dna, orf)) == NULL)
      {
         returnStatus = FALSE;
         break;
      }

      if(showDNA)
      {
         blWriteFASTA(out, header, dna, 60, FALSE);
      }
      else if(showRF)
      {
         blWriteFASTA(out, header, orf, 60, FALSE);
      }
      
      blWriteFASTA(out, header, protein, 60, FALSE);
      
      FREE(orf);
      FREE(protein);
      FREE(dna);
   }
   
   FREE(orf);
   FREE(protein);
   FREE(dna);
   
   return(returnStatus);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     BOOL *showDNA, BOOL *showRF)
   ---------------------------------------------------------------------
*//**
   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *showDNA     Display the original DNA sequence
   \param[out]     *showRF      Display the reading frame DNA
   \return                      Success?

   Parse the command line
   
-  03.11.17 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *showDNA, BOOL *showRF)
{
   argc--;
   argv++;

   *showDNA  = *showRF    = FALSE;
   infile[0] = outfile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'd':
            *showDNA = TRUE;
            break;
         case 'r':
            *showRF = TRUE;
            break;
         case 'h':
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there are 2 arguments left                    */
         if(argc > 2)
            return(FALSE);
         
         /* Copy the next to infile                                  */
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

-  10.11.17  Original   By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\nsixft V1.0 (c) 2017, UCL, Dr. Andrew C.R. \
Martin\n");
   fprintf(stderr,"\nUsage: sixft [-d|-r] [dna.faa [protein.faa]]\n");
   fprintf(stderr,"       -d Output the original DNA as well\n");
   fprintf(stderr,"       -r Output the DNA, but only the reading \
frame\n");
           
   fprintf(stderr,"\nPerform simple six-frame translation displaying \
only the longest\n");
   fprintf(stderr,"translation\n\n");
}

