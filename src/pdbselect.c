/************************************************************************/
/**

   \file       pdbselect.c
   
   \version    V1.0
   \date       30.02.15
   \brief      Select alterbative occupancies of models from a PDB file
   
   \copyright  (c) Dr. Andrew C. R. Martin 2015
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
-  V1.0  30.02.15 Original

*************************************************************************/
/* Includes
*/
#include <stdio.h>

#include "bioplib/SysDefs.h"
#include "bioplib/pdb.h"

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
void Usage(void);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *OccRank, int *ModelNum, BOOL *getInfo);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program

-  26.02.15 Original    By: ACRM
*/
int main(int argc, char **argv)
{
   FILE     *in      = stdin,
            *out     = stdout;
   WHOLEPDB *wpdb;
   int      OccRank  = 1,
            ModelNum = 1;
   char     infile[MAXBUFF],
            outfile[MAXBUFF];
   BOOL     getInfo = FALSE;

   if(ParseCmdLine(argc, argv, infile, outfile, &OccRank, &ModelNum, &getInfo))
   {
      if(blOpenStdFiles(infile, outfile, &in, &out))
      {
         BOOL DoWhole  = TRUE,
              AllAtoms = TRUE;
         
         if(((wpdb=blDoReadPDB(in, AllAtoms, OccRank, ModelNum, DoWhole))!=NULL) &&
            (wpdb->pdb != NULL))
         {
            if(getInfo)
            {
               BOOL printed = FALSE;
               
               if(gPDBMultiNMR)
               {
                  printed=TRUE;
                  fprintf(stderr,"PDB file contains %d models\n", gPDBMultiNMR);
               }
               
               if(gPDBPartialOcc)
               {
                  printed=TRUE;
                  fprintf(stderr,"PDB file contains partial occupancies\n");
               }
               
               if(!printed)
               {
                  fprintf(stderr,"PDB file does not contain partial \
occupancies or multiple models\n");
               }
            }
            else
            {
               blWriteWholePDB(out, wpdb);
            }
            
         }
         else
         {
            if(gPDBModelNotFound)
            {
               fprintf(stderr,"Requested model number not found: %d\n", ModelNum);
            }
            else
            {
               fprintf(stderr,"No atoms read from PDB file\n");
            }
            
            return(1);
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
/*>void Usage(void)
   ----------------
*//**

-  30.02.15 Original    By: ACRM
*/
void Usage(void)
{
   fprintf(stderr,"\npdbselect V1.0  (c) 2015 UCL, Andrew C.R. \
Martin\n");
   fprintf(stderr,"Usage: pdbselect [-i] [-o occupancy] [-m model] \
[<in.pdb> [<out.pdb>]]\n");
   fprintf(stderr,"       -i Print information on partial occupancy \
and models\n");
   fprintf(stderr,"       -o Specify the occupancy rank [Default: 1]\n");
   fprintf(stderr,"       -m Specify the model number [Default: 1]\n");

   fprintf(stderr,"With no command line options, this program simply \
reads a PDB file\n");
   fprintf(stderr,"and writes it out with just the first model and \
highest occupancy atoms\n");
   fprintf(stderr,"(as all bioptools programs do). The -o and -m flags \
allow different\n");
   fprintf(stderr,"occupancy ranks (2, 3, etc) and different models to \
be extracted\n");
   fprintf(stderr,"I/O is to stdin/stdout if not specified\n\n");
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     int *OccRank, int *ModelNum, BOOL *getInfo)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *OccRank     Occupancy rank
   \param[out]     *ModelNum    Model number
   \param[out]     *getInfo     Just get information
   \return                      Success?

   Parse the command line
   
-  30.02.15 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  int *OccRank, int *ModelNum, BOOL *getInfo)

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
         case 'o':
            argc--; argv++;
            if(!argc) return(FALSE);
            if(!sscanf(argv[0], "%d", OccRank)) return(FALSE);
            break;
         case 'm':
            argc--; argv++;
            if(!argc) return(FALSE);
            if(!sscanf(argv[0], "%d", ModelNum)) return(FALSE);
            break;
         case 'i':
            *getInfo = TRUE;
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

