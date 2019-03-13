/************************************************************************/
/**

   \file       setpdbnumbering.c
   
   \version    V1.7
   \date       13.03.19
   \brief      Apply standard numbering to a set of PDB files
   
   \copyright  (c) Dr. Andrew C. R. Martin 1996-2019
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
   Takes a PIR sequence alignment file as input. This must have the
   corresponding PDB filenames in the comment line before the - sign
   (or as the sole contents of the line). The PDB files are read and
   rewritten using the numbering scheme of the *first* one in the
   PIR file. Insert codes are added as necessary.

   The files are written out with the added .num extension.

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

-  V1.1  22.07.14 Renamed deprecated functions with bl prefix.
                  Added doxygen annotation. By: CTP
-  V1.2  06.11.14 Renamed from numpdb By: ACRM
-  V1.3  25.11.14 Initialized a variable
-  V1.4  05.03.15 Now calls pdbpatchnumbering rather than assuming the
                  link to patchpdbnum is available. Improved help message
-  V1.5  12.03.15 Changed to allow multi-character chain names
-  V1.6  28.01.18 Increased label buffer sizes
-  V1.7  13.03.19 Increase buffer sizes

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"
#include "bioplib/general.h"
#include "bioplib/pdb.h"
#include "bioplib/seq.h"
#include "bioplib/macros.h"
#include "bioplib/array.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXCHAIN  80
#define MAXBUFF  160

typedef struct _namseq
{
   struct _namseq *next;
   char   name[MAXBUFF],
          *seq;
}  NAMSEQ;


/************************************************************************/
/* Globals
*/
char gPatchFile[MAXBUFF];


/************************************************************************/
/* Prototypes
*/
int main(int argc, char **argv);
void Usage(void);
char *BuildSeqString(char **seqs, int nchain);
BOOL ParseCmdLine(int argc, char **argv, char *infile);
NAMSEQ *ReadSequenceData(FILE *in, int *nres);
BOOL GetNumbering(NAMSEQ *namseq, char **numbering);
BOOL ApplyNumbering(NAMSEQ *namseq, char **Numbering);
void BumpLabel(char *label);


/************************************************************************/
/*>int main(int argc, char **argv)
   -------------------------------
*//**

   Main program for applying numbering to a set of PDB files

-  05.02.96 Original    By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
*/
int main(int argc, char **argv)
{
   FILE    *in  = stdin,
           *out = stdout;
   char    AlnFile[MAXBUFF],
           **Numbering;
   NAMSEQ  *namseq;
   int     nres;
   pid_t   pid;

   /* Save a unique name for temp files                                 */
   pid = getpid();
   sprintf(gPatchFile,"/tmp/Patch.in.%d",pid);

   if(ParseCmdLine(argc, argv, AlnFile))
   {
      if(blOpenStdFiles(AlnFile, NULL, &in, &out))
      {
         if((namseq = ReadSequenceData(in, &nres))!=NULL)
         {
            /* Allocate memory for numbering                            */
            if((Numbering=(char **)blArray2D(sizeof(char),nres,8))==NULL)
            {
               fprintf(stderr,"No memory to store numbering\n");
               return(1);
            }
            
            if(GetNumbering(namseq, Numbering))
            {
               ApplyNumbering(namseq, Numbering);
            }
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

   Prints a usage message

-  05.02.96 Original    By: ACRM
-  22.07.14 V1.1 By: CTP
-  06.11.14 V1.2 By: ACRM
-  25.11.14 V1.3 By: ACRM
-  05.03.15 V1.4
-  12.03.15 V1.5
-  28.01.18 V1.6
-  13.03.19 V1.7
*/
void Usage(void)
{
   fprintf(stderr,"\nsetpdbnumbering V1.7 (c) 1996-2019 Dr. Andrew C.R. \
Martin, UCL\n");

   fprintf(stderr,"\nUsage: setpdbnumbering alnfile\n");

   fprintf(stderr,"\nApplies a standard numbering scheme to a set of PDB \
files. The input \n");
   fprintf(stderr,"'alnfile' is an alignment file in PIR format where \
the comment line for \n");
   fprintf(stderr,"each sequence entry contains the name of the input \
PDB file. The first\n");
   fprintf(stderr,"PDB file will be used to supply the numbering scheme; \
insertion codes\n");
   fprintf(stderr,"will be supplied relative to this file for the other \
structures.\n");

   fprintf(stderr,"\nAll you need as input it a PIR style alignment file \
with the name of\n");
   fprintf(stderr,"each PDB file in the comment line:\n");

   fprintf(stderr,"\n   >P1;1abc\n");
   fprintf(stderr,"   pdb1abc.ent\n");
   fprintf(stderr,"   ACTDFGIDEFGH--LIPNQRST-VLY*\n");
   fprintf(stderr,"   >P2;2def\n");
   fprintf(stderr,"   pdb2def.ent\n");
   fprintf(stderr,"   ACSEYG--EFGRTLLVPQQKSSRVLY*\n");

   fprintf(stderr,"\n2def will then be read and rewritten as pdb2def.num \
using the numbering\n");
   fprintf(stderr,"scheme of 1abc (The file can contain a multiple \
alignment - everything\n");
   fprintf(stderr,"will be written numbered according to 1abc.)\n");

   fprintf(stderr,"\nNote that the program makes use of \
pdbpatchnumbering program which must\n");
   fprintf(stderr,"be in your path.\n\n");
}


/************************************************************************/
/*>char *BuildSeqString(char **seqs, int nchain)
   ---------------------------------------------
*//**

   Takes a multi-chain sequence specification and builds it into a single
   string. Returns a malloc'd character pointer

-  05.02.96 Original    By: ACRM
*/
char *BuildSeqString(char **seqs, int nchain)
{
   int  i, 
        seqlen;
   char *string;
   
   /* Find the total length of the string                               */
   seqlen = 0;
   for(i=0; i<nchain; i++)
      seqlen += strlen(seqs[i]);

   /* Allocate this much space                                          */
   if((string = (char *)malloc(seqlen * sizeof(char)))==NULL)
      return(NULL);
   
   /* Copy the strings into the buffer                                  */
   string[0] = '\0';
   for(i=0;i<nchain;i++)
      strcat(string,seqs[i]);

   /* Return the buffer pointer                                         */
   return(string);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile)
   ------------------------------------------------------
*//**

   \param[in]      argc         Argument count
   \param[in]      **argv       Argument array
   \param[out]     *infile      Input file (or blank string)
   \return                     Success?

   Parse the command line
   
-  05.02.96 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile)
{
   argc--;
   argv++;

   infile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there is only 1 argument left                    */
         if(argc > 1)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}


/************************************************************************/
/*>NAMSEQ *ReadSequenceData(FILE *in, int *nres)
   ---------------------------------------------
*//**

   Reads a PIR file and builds a linked list of names (read from the
   comment line) and sequence strings

-  05.02.96 Original    By: ACRM
-  25.11.14 Initialize ns
*/
NAMSEQ *ReadSequenceData(FILE *in, int *nres)
{         
   char    *seqs[MAXCHAIN],
           *seqstring;
   SEQINFO seqinfo;
   int     nchain,
           NLastStruc = 0;
   BOOL    punct, error;
   NAMSEQ  *namseq = NULL,
           *ns     = NULL;

   while((nchain = blReadPIR(in, TRUE, seqs, MAXCHAIN, 
                             &seqinfo, &punct, &error))!=0)
   {
      seqstring = BuildSeqString(seqs,nchain);
      if(namseq == NULL)
      {
         INIT(namseq, NAMSEQ);
         ns = namseq;
      }
      else
      {
         ALLOCNEXT(ns, NAMSEQ);
      }
      
      if(ns==NULL || seqstring==NULL)
      {
         fprintf(stderr,"No memory for sequence list\n");
         return(NULL);
      }
      
      ns->seq = seqstring;
      strcpy(ns->name, seqinfo.name);

      if(NLastStruc)
      {
         if(NLastStruc != strlen(seqstring))
         {
            fprintf(stderr,"Your alignment file must contain sequences \
of identical length (once\n");
            fprintf(stderr,"the alignments have been made with - \
characters.\n");
         }
      }
      else
      {
         NLastStruc = strlen(seqstring);
      }
   }

   *nres = NLastStruc;
   return(namseq);
}


/************************************************************************/
/*>BOOL GetNumbering(NAMSEQ *namseq, char **numbering)
   ---------------------------------------------------
*//**

   Gets the numbering out of the first PDB file adding insert codes for
   for the sequence from the alignments

-  05.02.96 Original    By: ACRM
-  22.07.14 Renamed deprecated functions with bl prefix. By: CTP
-  12.03.15 Changed to allow multi-character chain names
-  28.01.18 Changed label[8] to label[32]
*/
BOOL GetNumbering(NAMSEQ *namseq, char **numbering)
{
   PDB  *pdb, *p, 
        *prev = NULL;
   FILE *fp;
   char label[32],
        nextlabel[32],
        InsertLabel = ' ',
        RevLabel    = 'Z';
   int  nres, rescount, natoms;
   
   /* Attempt to open the reference PDB file                            */
   if((fp=fopen(namseq->name, "r"))==NULL)
   {
      fprintf(stderr,"Unable to open reference PDB file: %s\n",
              namseq->name);
      return(FALSE);
   }

   /* Read the PDB file                                                 */
   if((pdb=blReadPDB(fp, &natoms))!=NULL)
   {
      nres = strlen(namseq->seq);
      p    = pdb;
      
      /* Step through the sequence array                                */
      for(rescount=0; rescount<nres; rescount++)
      {
         if(namseq->seq[rescount] == '-')
         {
            /* We hit a deletion c.f. the PDB linked list.
               p is the next residue in the PDB; build its label
            */
            if(prev == NULL)
            {
               sprintf(label, "%s.0 ",p->chain);
            }
            else
            {
               sprintf(label,"%s.%d%c",
                       prev->chain,prev->resnum,prev->insert[0]);
            }
            
            if(p!=NULL && (p->insert[0] != ' '))
            {
               /* If this has an insert code we have a problem; 
                  warn the user that insert code from the end of the
                  alphabet will be used
               */
               sprintf(nextlabel,"%s.%d%c",
                       p->chain,p->resnum,p->insert[0]);
               printf("Warning: Insertion occurs before residue %s\n",
                      nextlabel);
               printf("         Will use insertion codes from the end of \
the alphabet\n");
               
               if(prev == NULL)
               {
                  sprintf(label, "%s.0%c",p->chain,RevLabel);
               }
               else
               {
                  sprintf(label,"%s.%d%c",
                          prev->chain,prev->resnum,RevLabel);
               }
               
               RevLabel--;
            }
            else
            {
               /* Simply use the next alphabetical insert code          */
               BumpLabel(&InsertLabel);

               if(prev == NULL)
               {
                  sprintf(label, "%s.0%c",p->chain,InsertLabel);
               }
               else
               {
                  sprintf(label,"%s.%d%c",
                          prev->chain,prev->resnum,InsertLabel);
               }
            }
            
            /* Store the label                                          */
            strcpy(numbering[rescount], label);
         }
         else
         {
            /* Sequence maps directly to PDB file                       */

            /* Build a label                                            */
            sprintf(label,"%s.%d%c",p->chain,p->resnum,p->insert[0]);
            
            /* Store the label                                          */
            strcpy(numbering[rescount], label);
            
            /* Reset the insert labels                                  */
            InsertLabel = p->insert[0];
            RevLabel    = 'Z';

            /* Step on to the next residue                              */
            prev = p;
            p=blFindNextResidue(p);
         }
         
      }  /* End of loop through sequence array                          */

      FREELIST(pdb, PDB);
   }  /* Read PDB file OK                                               */

   return(TRUE);
}


/************************************************************************/
/*>BOOL ApplyNumbering(NAMSEQ *namseq, char **Numbering)
   -----------------------------------------------------
*//**

   Applies the numbering scheme to each PDB file in turn by calling out
   to the pdbpatchnumbering program to do the work

-  05.02.96 Original    By: ACRM
-  28.01.18 Increased buffer size to *3
-  13.03.19 Increased buffer size to *4
*/
BOOL ApplyNumbering(NAMSEQ *namseq, char **Numbering)
{
   NAMSEQ *ns;
   FILE   *fp;
   int    i;
   char   buffer[MAXBUFF*4];
   
   for(ns=namseq; ns!=NULL; NEXT(ns))
   {
      /* Open a patch file with a unique id                             */
      if((fp=fopen(gPatchFile,"w"))==NULL)
      {
         fprintf(stderr,"Unable to open temp file (%s) for writing.\n",
                 gPatchFile);
         return(FALSE);
      }

      /* Write the patch file                                           */
      for(i=0; ns->seq[i]; i++)
      {
         if(ns->seq[i] != '-')
            fprintf(fp,"%s %c\n",Numbering[i], ns->seq[i]);
      }
      fclose(fp);

      /* Run the patch program                                          */
      sprintf(buffer,"pdbpatchnumbering %s %s %s.num",gPatchFile,
              ns->name,ns->name);
      system(buffer);
   }

   /* Remove the patch file                                             */
   unlink(gPatchFile);

   return(TRUE);
}

      
/************************************************************************/
/*>void BumpLabel(char *label)
   ---------------------------
*//**

   Bumps an insert label.

-  05.02.96 Original    By: ACRM
*/
void BumpLabel(char *label)
{
   if(*label == ' ')
      *label = 'A';
   else
      (*label)++;
}

