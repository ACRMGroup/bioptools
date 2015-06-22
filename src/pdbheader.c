/************************************************************************/
/**

   \file       pdbheader.c
   
   \version    V1.3
   \date       22.06.15
   \brief      Get header info from a PDB file
   
   \copyright  (c) UCL / Dr. Andrew C.R. Martin, 2015
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

   See documentation for details

**************************************************************************

   Revision History:
   =================
-  V1.0  28.04.15 Original
-  V1.1  29.04.15 Added -p and fixed bug in -m
-  V1.2  04.06.15 Fixed bug in -c
-  V1.3  22.06.15 Added resolution info to header. With -r gives only this

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF     160
#define SMALLBUFF   16
#define TYPE_INT    0
#define TYPE_STRING 1
#define TYPE_REAL   2

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *chain, BOOL *doAll, BOOL *doSpecies,
                  BOOL *doMolecule, BOOL *noChains, BOOL *showPDB,
                  BOOL *resolOnly);
void Usage(void);
void PrintValue(FILE *fp, char *label, int width, int type, 
                char *string, int intval, REAL realval);
void ProcessChain(FILE *out, WHOLEPDB *wpdb, char *chain, 
                  BOOL doAll, BOOL doSpecies, BOOL doMolecule, 
                  BOOL showPDB);

/************************************************************************/
int main(int argc, char **argv)
{
   WHOLEPDB  *wpdb;
   FILE      *in = stdin,
             *out = stdout;
   char      header[MAXBUFF],
             infile[MAXBUFF],
             outfile[MAXBUFF],
             date[SMALLBUFF],
             chain[SMALLBUFF],
             pdbcode[SMALLBUFF],
             *title,
             **chainLabels;
   int       nChains,
             StrucType,
             i;
   BOOL      doAll      = TRUE,
             doSpecies  = FALSE,
             doMolecule = FALSE,
             noChains   = FALSE,
             showPDB    = FALSE,
             resolOnly  = FALSE;
   REAL      resolution,
             RFactor,
             FreeR;

   if(!ParseCmdLine(argc, argv, infile, outfile, chain, &doAll,
                    &doSpecies, &doMolecule, &noChains, &showPDB,
                    &resolOnly))
   {
      Usage();
      return(0);
   }
   
   if(!blOpenStdFiles(infile, outfile, &in, &out))
   {
      fprintf(stderr,"Error (pdbheader): Unable to open input or output \
file.\n");
      return(1);
   }
   
   if((wpdb = blReadWholePDB(in))!=NULL)
   {
      if(doAll)
      {
         if(!resolOnly)
         {
            if(blGetHeaderWholePDB(wpdb, 
                                   header, MAXBUFF,
                                   date,   SMALLBUFF,
                                   pdbcode, SMALLBUFF))
            {
               PrintValue(out, "PDB code:", 17, TYPE_STRING, pdbcode, 
                          0, 0.0);
               PrintValue(out, "Header:",   17, TYPE_STRING, header,  
                          0, 0.0);
               PrintValue(out, "Date:",     17, TYPE_STRING, date,    
                          0, 0.0);
            }
            
            if((title = blGetTitleWholePDB(wpdb))!=NULL)
            {
               PrintValue(out, "Title:",    17, TYPE_STRING, title,   
                          0, 0.0);
               free(title);
            }
         }

      }

      if(doAll || resolOnly)
      {
         if(blGetExptlWholePDB(wpdb, &resolution, &RFactor, &FreeR, 
                               &StrucType))
         {
            char *sTypeStr = blReportStructureType(StrucType);
            
            PrintValue(out, "Type:",       17, TYPE_STRING, sTypeStr, 
                       0, 0.0);
            PrintValue(out, "Resolution:", 17, TYPE_REAL, NULL,
                       0, resolution);
            PrintValue(out, "R-Factor:",   17, TYPE_REAL, NULL,
                       0, RFactor);
            PrintValue(out, "R-Free:",     17, TYPE_REAL, NULL,    
                       0, FreeR);
         }
      }

      if(!noChains)
      {
         if(chain[0])
         {
            ProcessChain(out, wpdb, chain, doAll, doSpecies, doMolecule,
                         showPDB);
         }
         else
         {
            chainLabels = blGetPDBChainLabels(wpdb->pdb, &nChains);
            for(i=0; i<nChains; i++)
            {
               ProcessChain(out, wpdb, chainLabels[i], 
                            doAll, doSpecies, doMolecule, showPDB);
               free(chainLabels[i]);
            }
            free(chainLabels);
         }
      }
      
   }

   return(0);
}


/************************************************************************/
/*>void ProcessChain(FILE *out, WHOLEPDB *wpdb, char *chain, 
                     BOOL doAll, BOOL doSpecies, BOOL doMolecule,
                     BOOL showPDB)
   --------------------------------------------------------------
*//**

   \param[in]   *out        File pointer
   \param[in]   *wpdb       Whole PDB structure
   \param[in]   *chain      Chain of interest
   \param[in]   doAll       Show all information for the chain
   \param[in]   doSpecies   Show the species information
   \param[in]   doMolecule  Show the molecule information
   \param[in]   showPDB     Show the PDB information

   Displays the information for a specified chain

-  28.04.15  Original   By: ACRM
-  04.06.15  If neither doSpecies nor doMolecule is specified then 
             set doAll
*/
void ProcessChain(FILE *out, WHOLEPDB *wpdb, char *chain, 
                  BOOL doAll, BOOL doSpecies, BOOL doMolecule,
                  BOOL showPDB)
{
   COMPND    compound;
   PDBSOURCE species;
   char      pdbcode[SMALLBUFF];

   pdbcode[0] = '\0';
   
   if(doAll)
   {
      PrintValue(out, "\n>Chain:",   18, TYPE_STRING, chain, 0, 0.0);
   }

   if(!doSpecies && !doMolecule)
   {
      doAll = TRUE;
   }
   
   if(showPDB)
   {
      char header[MAXBUFF],
           date[SMALLBUFF];
      
      blGetHeaderWholePDB(wpdb, header, MAXBUFF,
                          date, SMALLBUFF,
                          pdbcode, SMALLBUFF);
   }
   
   if(doAll || doMolecule)
   {
      if(blGetCompoundWholePDBChain(wpdb, chain, &compound))
      {
         if(doMolecule)
         {
            if(showPDB)
               fprintf(out,"%s : ", pdbcode);
            
            fprintf(out, "MOLECULE : %s : %s\n", 
                    chain, compound.molecule);
         }
         else if(doAll)
         {
            PrintValue(out,"MolID:",      17,  TYPE_INT,    
                       NULL,                compound.molid, 0.0);
            PrintValue(out,"Molecule:",   17,  TYPE_STRING, 
                       compound.molecule,   0, 0.0);
            PrintValue(out,"Fragment:",   17,  TYPE_STRING, 
                       compound.fragment,   0, 0.0);
            PrintValue(out,"Synonym:",    17, TYPE_STRING, 
                       compound.synonym,    0, 0.0);
            PrintValue(out,"EC:",         17, TYPE_STRING, 
                       compound.ec,         0, 0.0);
            PrintValue(out,"Engineered:", 17, TYPE_STRING, 
                       compound.engineered, 0, 0.0);
            PrintValue(out,"Mutation:",   17, TYPE_STRING, 
                       compound.mutation,   0, 0.0);
            PrintValue(out,"Other:",      17, TYPE_STRING, 
                       compound.other,      0, 0.0);
         }
      }
   }
   
   if(doAll || doSpecies)
   {
      if(blGetSpeciesWholePDBChain(wpdb, chain, &species))
      {
         if(doSpecies)
         {
            if(showPDB)
               fprintf(out,"%s : ", pdbcode);
            fprintf(out, "SPECIES  : %s : %s\n", 
                    chain, species.scientificName);
         }
         else
         {
            PrintValue(out,"Scientific name:", 17, TYPE_STRING, 
                       species.scientificName, 0, 0.0);
            PrintValue(out,"Common name:",     17, TYPE_STRING, 
                       species.commonName,     0, 0.0);
            PrintValue(out,"Strain:",          17, TYPE_STRING, 
                       species.strain,         0, 0.0);
            PrintValue(out,"Tax ID:",          17, TYPE_INT,    
                       NULL,                   species.taxid, 0.0);
         }
      }
   }
}


/************************************************************************/
/*>void PrintValue(FILE *fp, char *label, int width, int type, 
                   char *string, int intval, REAL realval)
   -----------------------------------------------------------
*//**
   \param[in]  *fp       File pointer
   \param[in]  *label    Label for printed output
   \param[in]  width     Width for output
   \param[in]  type      Type (TYPE_INT, TYPE_STRING, or TYPE_REAL)
   \param[in]  *string   String to print
   \param[in]  intval    Int to print
   \param[in]  realval   Real to print

   Prints a value with a label before it 

-  28.04.15  Original   By: ACRM
-  22.06.15  Added realval
*/
void PrintValue(FILE *fp, char *label, int width, int type, 
                char *string, int intval, REAL realval)
{
   char format[MAXBUFF];
   switch(type)
   {
   case TYPE_INT:
      sprintf(format, "%%-%ds%%d\n",width);
      fprintf(fp, format, label, intval);
      break;
   case TYPE_STRING:
      if((string!=NULL) && *string)
      {
         sprintf(format, "%%-%ds%%s\n",width);
         fprintf(fp, format, label, string);
      }
      break;
   case TYPE_REAL:
      sprintf(format, "%%-%ds%%.3f\n",width);
      fprintf(fp, format, label, realval);
      break;
   }
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                     char *chain, BOOL *doAll, BOOL *doSpecies,
                     BOOL *doMolecule, BOOL *noChains, BOOL *showPDB,
                     BOOL *resolOnly)
   ----------------------------------------------------------------------
*//**
   \param[in]   int    argc              Argument count
   \param[in]   char   **argv            Argument array
   \param[out]  char   *infile           Input filename (or blank string)
   \param[out]  char   *outfile          Output filename (or blank string)
   \param[out]  char   *chain            Chain label
   \param[out]  BOOL   *doAll            Print all the normal header 
                                         info
   \param[out]  BOOL   *doSpecies        Print the species concisely
   \param[out]  BOOL   *doMolecule       Print the molecule concisely
   \param[out]  BOOL   *noChains         Do not print chain info
   \param[out]  BOOL   *showPDB          Show PDB code with -m or -s
   \param[out]  BOOL   *resolOnly        Show only resolution info
   \return      BOOL                     Success

   Parse the command line

   28.04.15 Original    By: ACRM
   29.04.15 Added -s and showPDB parameter
   22.06.15 Added -r and resolOnly parameter
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile, 
                  char *chain, BOOL *doAll, BOOL *doSpecies,
                  BOOL *doMolecule, BOOL *noChains, BOOL *showPDB,
                  BOOL *resolOnly)
{
   argc--;
   argv++;
   
   infile[0]   = outfile[0] = chain[0] = '\0';
   *doSpecies  = FALSE;
   *doMolecule = FALSE;
   *noChains   = FALSE;
   *doAll      = TRUE;
   *showPDB    = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 's':
            *doSpecies = TRUE;
            *doAll     = FALSE;
            break;
         case 'm':
            *doMolecule = TRUE;
            *doAll      = FALSE;
            break;
         case 'c':
            argc--; argv++;
            if(!argc) return(FALSE);
            strncpy(chain, argv[0], SMALLBUFF);
            *doAll      = FALSE;
            break;
         case 'n':
            *noChains   = TRUE;
            break;
         case 'p':
            *showPDB    = TRUE;
            break;
         case 'r':
            *resolOnly  = TRUE;
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

-   28.04.15 Original   By: ACRM
-   29.04.15 Added -p
*/
void Usage(void)
{
   fprintf(stderr,"\npdbheader V1.3 (c) 2015 UCL, Dr. Andrew C.R. \
Martin\n");
   fprintf(stderr,"Usage: pdbheader [-s][-m][-p][-c chain][-n][-r] \
[in.pdb [out.pdb]]\n");
   fprintf(stderr,"       -s Show species information rather than \
everything\n");
   fprintf(stderr,"       -m Show molecule information rather than \
everything\n");
   fprintf(stderr,"       -p Show PDB code with -m or -s\n");
   fprintf(stderr,"       -c Only do the specified chain\n");
   fprintf(stderr,"       -n Do not show chain information - just the \
main header\n");
   fprintf(stderr,"       -r Only show resolution information\n");

   fprintf(stderr,"\nParses and displays the header information from a \
PDB file. The default\n");
   fprintf(stderr,"is to show the HEADER and TITLE information for the \
file followed by\n");
   fprintf(stderr,"COMPND and SOURCE information for each chain. With \
-n, no chain\n");
   fprintf(stderr,"information is shown. With -c, only the information \
for the specified\n");
   fprintf(stderr,"chain is shown. The -s and -m options lead to a more \
compact view of\n");
   fprintf(stderr,"the species and molecule information for the chains \
with no general\n");
   fprintf(stderr,"header information.\n\n");
}

