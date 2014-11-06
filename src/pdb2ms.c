/************************************************************************/
/**

   \file       pdb2ms.c
   
   \version    V1.3
   \date       22.07.14
   \brief      Create input file for Connoly MS program
   
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
   Converts a PDB file to input format for Connoly's MS program
   Optionally writes standard data files as well.

   All data are hard coded for ease of distribution.

**************************************************************************

   Usage:
   ======

   NOTE

   You need the BiopLib library to compile this program. Contact
   the author to obtain Bioplib, or look on the Web page:
   http://www.biochem.ucl.ac.uk/~martin

**************************************************************************

   Revision History:
   =================
-  V1.0  24.01.96 Original
-  V1.1  29.01.96 Added -q flag
-  V1.2  01.02.96 Has default types if not standard AAs
                  Can take atom types or radii from the PDB file
-  V1.3  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include "bioplib/pdb.h"
#include "bioplib/macros.h"
#include "bioplib/SysDefs.h"
#include "bioplib/MathType.h"
#include "bioplib/general.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF      160
#define ALLOCQUANTUM  10

typedef struct
{
   char *resnam,
        *atnam;
   int  attype1, attype2;
}  ATOMTYPES;


/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
int  main(int argc, char **argv);
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *DoStd, BOOL *Quiet, BOOL *Alt, BOOL *GotRad,
                  BOOL *GotType);
void Usage(void);
BOOL ConvertPDB2MS(FILE *out, PDB *pdb, BOOL Alt, BOOL GotRad, 
                   BOOL GotType);
int  AtomType(char *resnam, char *atnam, BOOL Alt);
void WriteStdDataFiles(BOOL Quiet, BOOL GotRad);
void WriteRadii(REAL *typerad);
REAL *StoreRadius(REAL *typerad, REAL radius, int *type);
int RadiusSeen(REAL *typerad, REAL radius);

/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for converting PDB format to MS 

-  24.01.96 Original   By: ACRM
-  29.01.96 Added -a and -q handling
-  01.02.95 Added -t and -r handling
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   char InFile[MAXBUFF],
        OutFile[MAXBUFF];
   FILE *in  = stdin,
        *out = stdout;
   PDB  *pdb;
   int  natoms;
   BOOL DoStd, Quiet, Alt, GotRad, GotType;

   if(ParseCmdLine(argc, argv, InFile, OutFile, &DoStd, &Quiet, &Alt,
                   &GotRad, &GotType))
   {
      if(blOpenStdFiles(InFile, OutFile, &in, &out))
      {
         if((pdb = blReadPDB(in, &natoms))==NULL)
         {
            fprintf(stderr,"No atoms read from PDB file\n");
            return(1);
         }

         if(!ConvertPDB2MS(out, pdb, Alt, GotRad, GotType))
            return(1);

         if(DoStd)
            WriteStdDataFiles(Quiet, GotRad);
      }
   }
   else
   {
      Usage();
   }

   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                     BOOL *DoStd, BOOL *Quiet, BOOL *Alt, BOOL *GotRad,
                     BOOL *GotType)
   ---------------------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \param[out]     *outfile     Output file (or blank string)
   \param[out]     *DoStd       Write standard data files?
   \param[out]     *Quiet       Operate quietly?
   \param[out]     *Alt         Use alternative atom type set 
   \param[out]     *GotRad      Got radii in BVal column
   \param[out]     *GotType     Got atom type numbers in BVal column
   \return                     Success?

   Parse the command line
   
-  25.01.96 Original    By: ACRM
-  29.01.96 Added -a and -q
-  01.02.95 Added -r and -t
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile, char *outfile,
                  BOOL *DoStd, BOOL *Quiet, BOOL *Alt, BOOL *GotRad,
                  BOOL *GotType)
{
   argc--;
   argv++;

   infile[0] = outfile[0] = '\0';
   *DoStd    = FALSE;
   *Quiet    = FALSE;
   *Alt      = FALSE;
   *GotRad   = FALSE;
   *GotType  = FALSE;
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 's':
            *DoStd = TRUE;
            break;
         case 'q':
            *Quiet = TRUE;
            break;
         case 'a':
            *Alt = TRUE;
            break;
         case 'r':
            *GotRad = TRUE;
            break;
         case 't':
            *GotType = TRUE;
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

-  24.01.96 Original   By: ACRM
-  29.01.95 V1.1
-  01.02.96 V1.2
-  22.07.14 V1.3 By: CTP
*/
void Usage(void)
{
   fprintf(stderr,"\npdb2ms V1.3 (c)1996-2014, Dr. Andrew C.R. Martin, UCL\n");
   fprintf(stderr,"\nUsage: pdb2ms [-s] [-a] [-q] [in.pdb [out.ms]]\n");
   fprintf(stderr,"       -s Write standard data files as well\n");
   fprintf(stderr,"       -a Use alternate atom type radii (as used by \
acall/asurf/access)\n");
   fprintf(stderr,"       -q Operate quietly\n");
   fprintf(stderr,"       -r Use BVal column as radius\n");
   fprintf(stderr,"       -t Use BVal column as atom type\n");
   fprintf(stderr,"\nConverts a PDB file to input for the Connoly MS \
program\n");
   fprintf(stderr,"If standard data files are written they are named \
control.dat and \n");
   fprintf(stderr,"radii.dat. The latter will always be written if the \
-r flag is used.\n");
   fprintf(stderr,"With -r and -t, the default values will be \
substituted if the BVal column\n");
   fprintf(stderr,"contains 0.00. This allows you to set the BVal column \
to 0.00 throughout\n");
   fprintf(stderr,"the structure, but to give explicit radii for certain \
atoms (for example\n");
   fprintf(stderr,"HETATMs).\n\n");
}


/************************************************************************/
/*>BOOL ConvertPDB2MS(FILE *out, PDB *pdb, BOOL Alt, BOOL GotRad, 
                      BOOL GotType)
   --------------------------------------------------------------
*//**

   Does the actual conversion work

-  24.01.96 Original   By: ACRM
-  29.01.96 Added Alt parameter
-  01.02.96 Added GotRad and GotType handling
*/
BOOL ConvertPDB2MS(FILE *out, PDB *pdb, BOOL Alt, BOOL GotRad, 
                   BOOL GotType)
{
   PDB  *p;
   int  type;
   REAL *typerad = NULL,
        radius;

   if(GotRad)
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         radius = p->bval;
         
         if(radius == (REAL)0.0)
         {
            switch(p->atnam[0])
            {
            case 'C':
               if(p->atnam[1] == ' ')
                  radius = 1.76;
               else
                  radius = 1.87;
               break;
            case 'N':
               radius = 1.65;
               break;
            case 'O':
               radius = 1.40;
               break;
            case 'S':
               radius = 1.85;
               break;
            case 'H':
               radius = 1.40;
               break;
            default:
               radius = 1.70;
               break;
            }
         }
         
         if((type=RadiusSeen(typerad, radius))==0)
         {
            if((typerad = StoreRadius(typerad, radius, &type))==NULL)
            {
               fprintf(stderr,"pdb2ms: No memory for radius table\n");
               return(FALSE);
            }
         }

         fprintf(out,"%10.5f%10.5f%10.5f%5d    2    1\n",
                 p->x, p->y, p->z, type);

         WriteRadii(typerad);
      }
   }
   else if(GotType)
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         type = (int)p->bval;
         if(type == 0)
            type = AtomType(p->resnam, p->atnam, Alt);

         fprintf(out,"%10.5f%10.5f%10.5f%5d    2    1\n",
                 p->x, p->y, p->z, type);
      }
   }
   else
   {
      for(p=pdb; p!=NULL; NEXT(p))
      {
         if((type = AtomType(p->resnam, p->atnam, Alt))==(-1))
         {
            fprintf(stderr,"Warning, Atom %s %s not known; \
type set to 1\n",
                    p->resnam, p->atnam);
            type = 1;
         }
         fprintf(out,"%10.5f%10.5f%10.5f%5d    2    1\n",
                 p->x, p->y, p->z, type);
      }
   }

   if(typerad != NULL)
      free(typerad);

   return(TRUE);
}


/************************************************************************/
/*>int AtomType(char *inresnam, char *atnam, BOOL Alt)
   ---------------------------------------------------
*//**

   Converts a residue name / atom name combination into an atom type
   number

-  24.01.96 Original   By: ACRM
-  29.01.96 Added Alt set
-  01.02.96 Returns default types if residue not found
*/
int AtomType(char *inresnam, char *atnam, BOOL Alt)
{
   int              i;
   char             resnam[8];
   static ATOMTYPES atomtypes[] =
   {  
      { "*   ", "OXT ", 8},

      { "ALA ", "C   ", 2, 35},
      { "ALA ", "CA  ", 1, 34},
      { "ALA ", "CB  ", 1, 34},
      { "ALA ", "N   ", 28, 32},
      { "ALA ", "O   ", 9, 36},

      { "ARG ", "C   ", 2, 35},
      { "ARG ", "CA  ", 1, 34},
      { "ARG ", "CB  ", 1, 34},
      { "ARG ", "CD  ", 1, 34},
      { "ARG ", "CG  ", 1, 34},
      { "ARG ", "CZ  ", 2, 34},
      { "ARG ", "N   ", 28, 32},
      { "ARG ", "NE  ", 19, 32},
      { "ARG ", "NH1 ", 19, 32},
      { "ARG ", "NH2 ", 19, 32},
      { "ARG ", "O   ", 9, 36},

      { "ASN ", "C   ", 2, 35},
      { "ASN ", "CA  ", 1, 34},
      { "ASN ", "CB  ", 1, 34},
      { "ASN ", "CG  ", 2, 35},
      { "ASN ", "N   ", 28, 32},
      { "ASN ", "ND2 ", 28, 32},
      { "ASN ", "O   ", 9, 36},
      { "ASN ", "OD1 ", 9, 36},

      { "ASP ", "C   ", 2, 35},
      { "ASP ", "CA  ", 1, 34},
      { "ASP ", "CB  ", 1, 34},
      { "ASP ", "CG  ", 2, 35},
      { "ASP ", "N   ", 28, 32},
      { "ASP ", "O   ", 9, 36},
      { "ASP ", "OD1 ", 9, 36},
      { "ASP ", "OD2 ", 9, 36},

      { "CYS ", "C   ", 2, 35},
      { "CYS ", "CA  ", 1, 34},
      { "CYS ", "CB  ", 1, 34},
      { "CYS ", "N   ", 28, 32},
      { "CYS ", "O   ", 9, 36},
      { "CYS ", "SG  ", 10, 37},

      { "GLU ", "C   ", 2, 35},
      { "GLU ", "CA  ", 1, 34},
      { "GLU ", "CB  ", 1, 34},
      { "GLU ", "CG  ", 1, 34},
      { "GLU ", "CD  ", 2, 35},
      { "GLU ", "N   ", 28, 32},
      { "GLU ", "O   ", 9, 36},
      { "GLU ", "OE1 ", 9, 36},
      { "GLU ", "OE2 ", 9, 36},

      { "GLY ", "C   ", 2, 35},
      { "GLY ", "CA  ", 1, 34},
      { "GLY ", "N   ", 28, 32},
      { "GLY ", "O   ", 9, 36},

      { "ILE ", "C   ", 2, 35},
      { "ILE ", "CA  ", 1, 34},
      { "ILE ", "CB  ", 1, 34},
      { "ILE ", "CD1 ", 1, 34},
      { "ILE ", "CG1 ", 1, 34},
      { "ILE ", "CG2 ", 1, 34},
      { "ILE ", "N   ", 28, 32},
      { "ILE ", "O   ", 9, 36},

      { "LEU ", "C   ", 2, 35},
      { "LEU ", "CA  ", 1, 34},
      { "LEU ", "CB  ", 1, 34},
      { "LEU ", "CD1 ", 1, 34},
      { "LEU ", "CD2 ", 1, 34},
      { "LEU ", "CG  ", 1, 34},
      { "LEU ", "N   ", 28, 32},
      { "LEU ", "O   ", 9, 36},

      { "PHE ", "C   ", 2, 35},
      { "PHE ", "CA  ", 1, 34},
      { "PHE ", "CB  ", 1, 34},
      { "PHE ", "CD1 ", 3, 35},
      { "PHE ", "CD2 ", 3, 35},
      { "PHE ", "CE1 ", 3, 35},
      { "PHE ", "CE2 ", 3, 35},
      { "PHE ", "CG  ", 3, 35},
      { "PHE ", "CZ  ", 3, 35},
      { "PHE ", "N   ", 28, 32},
      { "PHE ", "O   ", 9, 36},

      { "PRO ", "C   ", 2, 35},
      { "PRO ", "CA  ", 1, 34},
      { "PRO ", "CB  ", 1, 34},
      { "PRO ", "CD  ", 1, 34},
      { "PRO ", "CG  ", 1, 34},
      { "PRO ", "N   ", 28, 32},
      { "PRO ", "O   ", 9, 36},

      { "SER ", "C   ", 2, 35},
      { "SER ", "CA  ", 1, 34},
      { "SER ", "CB  ", 1, 34},
      { "SER ", "N   ", 28, 32},
      { "SER ", "O   ", 9, 36},
      { "SER ", "OG  ", 8, 36},

      { "THR ", "C   ", 2, 35},
      { "THR ", "CA  ", 1, 34},
      { "THR ", "CB  ", 1, 34},
      { "THR ", "CG2 ", 1, 34},
      { "THR ", "N   ", 28, 32},
      { "THR ", "O   ", 9, 36},
      { "THR ", "OG1 ", 8, 36},

      { "TYR ", "C   ", 2, 35},
      { "TYR ", "CA  ", 1, 34},
      { "TYR ", "CB  ", 1, 34},
      { "TYR ", "CD1 ", 3, 35},
      { "TYR ", "CD2 ", 3, 35},
      { "TYR ", "CE1 ", 3, 35},
      { "TYR ", "CE2 ", 3, 35},
      { "TYR ", "CG  ", 3, 35},
      { "TYR ", "CZ  ", 3, 35},
      { "TYR ", "N   ", 28, 32},
      { "TYR ", "O   ", 9, 36},
      { "TYR ", "OH  ", 8, 36},

      { "VAL ", "C   ", 2, 35},
      { "VAL ", "CA  ", 1, 34},
      { "VAL ", "CB  ", 1, 34},
      { "VAL ", "CG1 ", 1, 34},
      { "VAL ", "CG2 ", 1, 34},
      { "VAL ", "N   ", 28, 32},
      { "VAL ", "O   ", 9, 36},

      { "HIS ", "N   ", 28, 32},
      { "HIS ", "CA  ", 1, 34},
      { "HIS ", "C   ", 2, 35},
      { "HIS ", "O   ", 9, 36},
      { "HIS ", "CB  ", 1, 34},
      { "HIS ", "CG  ", 1, 35},
      { "HIS ", "ND1 ", 28, 32},
      { "HIS ", "CD2 ", 3, 35},
      { "HIS ", "CE1 ", 3, 35},
      { "HIS ", "NE2 ", 28, 32},

      { "LYS ", "N   ", 28, 32},
      { "LYS ", "CA  ", 1, 34},
      { "LYS ", "C   ", 2, 35},
      { "LYS ", "O   ", 9, 36},
      { "LYS ", "CB  ", 1, 34},
      { "LYS ", "CG  ", 1, 34},
      { "LYS ", "CD  ", 1, 34},
      { "LYS ", "CE  ", 2, 34},
      { "LYS ", "NZ  ", 19, 7},

      { "MET ", "N   ", 28, 32},
      { "MET ", "CA  ", 1, 34},
      { "MET ", "C   ", 2, 35},
      { "MET ", "O   ", 9, 36},
      { "MET ", "CB  ", 1, 34},
      { "MET ", "CG  ", 1, 34},
      { "MET ", "SD  ", 10, 37},
      { "MET ", "CE  ", 1, 34},

      { "GLN ", "N   ", 28, 32},
      { "GLN ", "CA  ", 1, 34},
      { "GLN ", "C   ", 2, 35},
      { "GLN ", "O   ", 9, 36},
      { "GLN ", "CB  ", 1, 34},
      { "GLN ", "CG  ", 1, 34},
      { "GLN ", "CD  ", 1, 35},
      { "GLN ", "OE1 ", 9, 36},
      { "GLN ", "NE2 ", 28, 32},

      { "TRP ", "N   ", 28, 32},
      { "TRP ", "CA  ", 1, 34},
      { "TRP ", "C   ", 2, 35},
      { "TRP ", "O   ", 9, 36},
      { "TRP ", "CB  ", 1, 34},
      { "TRP ", "CG  ", 1, 35},
      { "TRP ", "CD1 ", 3, 35},
      { "TRP ", "CD2 ", 3, 35},
      { "TRP ", "NE1 ", 28, 32},
      { "TRP ", "CE2 ", 3, 35},
      { "TRP ", "CE3 ", 3, 35},
      { "TRP ", "CZ2 ", 3, 35},
      { "TRP ", "CZ3 ", 3, 35},
      { "TRP ", "CH2 ", 3, 35},

      { "PCA ", "C   ", 2, 35},
      { "PCA ", "CA  ", 1, 34},
      { "PCA ", "CB  ", 1, 34},
      { "PCA ", "CD  ", 1, 24},
      { "PCA ", "CG  ", 1, 25},
      { "PCA ", "N   ", 28, 32},
      { "PCA ", "O   ", 9, 36},
      { "PCA ", "OE  ", 8, 35},

      { NULL,   NULL,   0, 0}
   };

   strcpy(resnam, inresnam);

   if(!strncmp(atnam,"OXT ",4))
   {
      strcpy(resnam, "*   ");
   }

   for(i=0; atomtypes[i].resnam != NULL; i++)
   {
      if(!strncmp(atnam,  atomtypes[i].atnam,  4) &&
         !strncmp(resnam, atomtypes[i].resnam, 4))
      {
         if(Alt)
            return(atomtypes[i].attype2);
         else
            return(atomtypes[i].attype1);
      }
   }

   /* We failed to find the required atoms in the special types, so
      give defaults:
   */
   switch(atnam[0])
   {
   case 'C':
      return((Alt)?34:1);
   case 'N':
      return((Alt)?32:28);
   case 'O':
      return((Alt)?36:8);
   case 'S':
      return((Alt)?37:10);
   case 'F':
      if(atnam[1] == 'E')
         return(38);
      /* Drop through                                                   */
   case 'H':
      return(23);
   default:
      return(10);
   }

   return(10);
}

/************************************************************************/
/*>void WriteStdDataFiles(BOOL Quiet, BOOL GotRad)
   -----------------------------------------------
*//**

   Writes a standard control file and radius file

-  24.01.96 Original   By: ACRM
-  29.01.96 Added Quiet flag
-  01.02.96 Added GotRad
            Changed default probe radius from 1.5 to 1.4
*/
void WriteStdDataFiles(BOOL Quiet, BOOL GotRad)
{
   FILE *fp;
   
   if((fp=fopen("control.dat","w"))==NULL)
   {
      fprintf(stderr,"Unable to open control.dat for writing\n");
   }
   else
   {
      fprintf(fp," 5.000     1.40         0    0\n");
      fclose(fp);

      if(!Quiet)
      {
         fprintf(stderr,"\nA standard control file has been written to \
control.dat\n");
         fprintf(stderr,"The fields are: density, probe-size, buried, \
format\n");
         fprintf(stderr,"Buried: 0=normal, 1=only surface buried by \
another molecule, 2=both\n");
         fprintf(stderr,"Format: 0=Long ASCII, 1=Long binary, 2=Short \
ASCII, 3=Short bindary\n\n");
      }
   }

   if(!GotRad)
   {
      if((fp=fopen("radii.dat","w"))==NULL)
      {
         fprintf(stderr,"Unable to open radii.dat for writing\n");
      }
      else
      {
         fprintf(fp,"    0     0.00\n");
         fprintf(fp,"    1     1.53\n");
         fprintf(fp,"    2     1.53\n");
         fprintf(fp,"    3     1.53\n");
         fprintf(fp,"    4     1.54\n");
         fprintf(fp,"    5     1.45\n");
         fprintf(fp,"    6     1.48\n");
         fprintf(fp,"    7     1.50\n");
         fprintf(fp,"    8     1.36\n");
         fprintf(fp,"    9     1.36\n");
         fprintf(fp,"   10     1.70\n");
         fprintf(fp,"   11     1.48\n");
         fprintf(fp,"   12     1.75\n");
         fprintf(fp,"   13     1.08\n");
         fprintf(fp,"   14     1.80\n");
         fprintf(fp,"   15     1.65\n");
         fprintf(fp,"   16     1.30\n");
         fprintf(fp,"   17     2.05\n");
         fprintf(fp,"   18     1.72\n");
         fprintf(fp,"   19     1.50\n");
         fprintf(fp,"   20     0.85\n");
         fprintf(fp,"   21     0.95\n");
         fprintf(fp,"   22     1.33\n");
         fprintf(fp,"   23     0.99\n");
         fprintf(fp,"   24     0.60\n");
         fprintf(fp,"   25     2.05\n");
         fprintf(fp,"   26     0.00\n");
         fprintf(fp,"   27     2.10\n");
         fprintf(fp,"   28     1.45\n");
         fprintf(fp,"   29     1.70\n");
         fprintf(fp,"   30     1.70\n");
         fprintf(fp,"   31     1.45\n");
         fprintf(fp,"   32     1.65\n");
         fprintf(fp,"   33     1.844\n");
         fprintf(fp,"   34     1.87\n");
         fprintf(fp,"   35     1.76\n");
         fprintf(fp,"   36     1.40\n");
         fprintf(fp,"   37     1.85\n");
         fprintf(fp,"   38     1.47\n");
         
         fclose(fp);

         if(!Quiet)
            fprintf(stderr,"A standard atom type radius file has been \
written to radii.dat\n\n");
      }
   }
}

/************************************************************************/
/*>int RadiusSeen(REAL *typerad, REAL radius)
   ------------------------------------------
*//**

   \param[in]      *typerad         Array of seen radii (or NULL)
   \param[in]      radius           Radius for which to search
   \return                     Array offset+1 (0 if not seen)

   Looks to see if given radius has been seen previously (i.e. is in the
   typerad array). Returns the array offset + 1; 0 used to indicate
   not found.

-  01.02.96 Original   By: ACRM
*/
int RadiusSeen(REAL *typerad, REAL radius)
{
   int i;
   
   if(typerad == NULL)
      return(0);

   for(i=0; typerad[i] > (REAL)0.0; i++)
   {
      if(typerad[i] == radius)
         return(i+1);
   }
   
   return(0);
}

/************************************************************************/
/*>REAL *StoreRadius(REAL *typerad, REAL radius, int *type)
   --------------------------------------------------------
*//**

   \param[in]      *typerad     Array of seen radii (or NULL)
   \param[in]      radius       New radius to store in array
   \param[out]     *type        Array offset+1 of stored value
   \return                      New array pointer (when allocations
                                have occurred)

   Allocate storage space (if needed) and store a radius in the typerad
   array.

-  01.02.96 Original   By: ACRM
*/
REAL *StoreRadius(REAL *typerad, REAL radius, int *type)
{
   static int MaxRad = 0,
              NRad   = 0;
   
   /* Allocate space if needed                                          */
   if((NRad + 2) > MaxRad)
   {
      MaxRad += ALLOCQUANTUM;
      if(typerad==NULL)
         typerad = (REAL *)malloc(MaxRad * sizeof(REAL));
      else
         typerad = (REAL *)realloc(typerad, MaxRad * sizeof(REAL));
      if(typerad==NULL)
         return(NULL);
   }
   
   /* Store the value                                                   */
   typerad[NRad++] = radius;
   *type = NRad;                    /* One more than the array position */
   typerad[NRad]   = (REAL)(-1.0);  /* Terminate the array              */

   return(typerad);
}

/************************************************************************/
/*>void WriteRadii(REAL *typerad)
   ------------------------------
*//**

   \param[in]      *typerad     Array of seen radii

   Write the observed radii to a file called radii.dat

-  01.02.96 Original   By: ACRM
*/
void WriteRadii(REAL *typerad)
{
   FILE *fp;
   int  i;
   
   if(typerad==NULL)
      return;
   
   if((fp=fopen("radii.dat","w"))==NULL)
   {
      fprintf(stderr,"Unable to open radii.dat for writing\n");
   }
   else
   {
      for(i=0; typerad[i] > (REAL)0.0; i++)
         fprintf(fp,"%5d%9.2f\n",i+1,typerad[i]);
      fclose(fp);
   }   
}



