/************************************************************************/
/**

   Program:    sixft
   \par
   \file       sexft.c
   
   \version    V1.1
   \date       02.12.22
   \brief      Six-frame translation
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2017-22
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
   V1.0    10.11.17 Original version
   V1.1    02.12.22 Added -a option to show all translations

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
BOOL DoTranslate(FILE *in, FILE *out, BOOL showDNA, BOOL showRF,
                 BOOL showAll);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *showDNA, BOOL *showRF, BOOL *showAll);
void Usage(void);
void DisplayResults(FILE *out, BOOL showDNA, BOOL showRF, char *header,
                    char *protein, char *dna, char *orf, int width);
char *translateAnyFrame(char *inDna, int frame);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**
   Main program

-  10.11.17 Original   By: ACRM
-  02.12.22 Added -a handling
*/
int main(int argc, char **argv)
{
   char infile[MAXBUFF],
        outfile[MAXBUFF];
   BOOL showDNA = FALSE,
        showRF  = FALSE,
        showAll = FALSE;
   FILE *in     = stdin,
        *out    = stdout;

   if(ParseCmdLine(argc, argv, infile, outfile, &showDNA, &showRF,
                   &showAll))
   {
      if(showDNA && showRF)
         showDNA = FALSE;
      if(showAll && showRF)
      {
         showRF  = FALSE;
         showDNA = TRUE;
      }
      

      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         if(!DoTranslate(in, out, showDNA, showRF, showAll))
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
/*>BOOL DoTranslate(FILE *in, FILE *out, BOOL showDNA, BOOL showRF,
                    BOOL showAll)
   ----------------------------------------------------------------
*//**
   \param[in]   in       Input file pointer
   \param[in]   out      Output file pointer
   \param[in]   showDNA  Display the DNA sequence
   \param[in]   showRF   Display the reading frame sequence
   \param[in]   showAll  Display all reading frame translations
   \return               Success?

   Does the actual work of performing the 6FT and finding the best
   translation

-  10.11.17 Original   By: ACRM
*/
BOOL DoTranslate(FILE *in, FILE *out, BOOL showDNA, BOOL showRF,
                 BOOL showAll)
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

      if(showAll)
      {
         int frame;
         for(frame=0; frame<6; frame++)
         {
            if((protein = translateAnyFrame(dna, frame)) == NULL)
            {
               returnStatus = FALSE;
               break;
            }
            DisplayResults(out, showDNA, showRF, header,
                           protein, dna, NULL, 60);
            FREE(protein);
         }
      }
      else
      {
         if((protein = blSixFTBest(dna, orf)) == NULL)
         {
            returnStatus = FALSE;
            break;
         }

         DisplayResults(out, showDNA, showRF, header,
                        protein, dna, orf, 60);
         FREE(protein);
      }
   }
   
   FREE(orf);
   FREE(protein);
   FREE(dna);
   
   return(returnStatus);
}


void DisplayResults(FILE *out, BOOL showDNA, BOOL showRF, char *header,
                    char *protein, char *dna, char *orf, int width)
{
   if(showDNA)
   {
      blWriteFASTA(out, header, dna, width, FALSE);
   }
   else if(showRF)
   {
      blWriteFASTA(out, header, orf, width, FALSE);
   }
   
   blWriteFASTA(out, header, protein, width, FALSE);

}



/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     BOOL *showDNA, BOOL *showRF, BOOL *showAll)
   ---------------------------------------------------------------------
*//**
   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *showDNA     Display the original DNA sequence
   \param[out]     *showRF      Display the reading frame DNA
   \param[out]     *showAll     Show all translations
   \return                      Success?

   Parse the command line
   
-  03.11.17 Original    By: ACRM
-  02.12.22 Added -a / showAll
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  BOOL *showDNA, BOOL *showRF, BOOL *showAll)
{
   argc--;
   argv++;

   *showDNA  = *showRF    = *showAll = FALSE;
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
         case 'a':
            *showAll = TRUE;
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
/*>char *translateAnyFrame(char *inDna, int frame)
   -----------------------------------------------
*//**
   \param[in]    inDna    DNA sequence
   \param[out]   orf      DNA of the ORF (if NULL, doesn't provide this)
   \param[in]    frame    Reading frame (0-2 forward, 3-5 reverse)
   \return                Malloc'd protein sequence

   Performs a 6FT and find the longest translation starting from the
   beginning of the DNA or a Met. Returns malloc'd memory for the 
   protein sequence. Optionally outputs the ORF that was translated.

-  02.12.22 Original   By: ACRM
*/
char *translateAnyFrame(char *inDna, int frame)
{
   char *protSeq     = NULL,
        *dna;

   /* For working copy of protein sequence                              */
   if((protSeq = (char *)malloc((2+strlen(inDna)/3) * sizeof(char)))
      == NULL)
      return(NULL);

   /* Put the DNA and RC DNA into an array                              */
   if(frame < 3)
   {
      dna = inDna;
   }
   else
   {
      dna = blReverseComplement(inDna);
   }
   
   if(dna == NULL)
   {
      free(protSeq);
      return(NULL);
   }
   
   blTranslateFrame(dna, frame%3, protSeq);

   if(frame >= 3)
   {
      FREE(dna);
   }
   
   return(protSeq);
}


/************************************************************************/
/*>void Usage(void)
   ----------------
*//**
   Prints a usage message

-  10.11.17  Original   By: ACRM
-  02.12.22  V1.1 added -a
*/
void Usage(void)
{
   fprintf(stderr,"\nsixft V1.1 (c) 2017-22, UCL, Prof. Andrew C.R. \
Martin\n");
   fprintf(stderr,"\nUsage: sixft [-d|-r][-a] [dna.faa [protein.faa]]\n");
   fprintf(stderr,"       -d Output the original DNA as well\n");
   fprintf(stderr,"       -r Output the DNA, but only the reading \
frame\n");
   fprintf(stderr,"       -a Show translations from all frames\n");
           
   fprintf(stderr,"\nPerform simple six-frame translation displaying \
only the longest\n");
   fprintf(stderr,"translation by default.\n\n");
}

