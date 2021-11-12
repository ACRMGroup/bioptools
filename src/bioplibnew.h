#ifndef _BIOPLIBNEW_H
#define _BIOPLIBNEW_H 1

STRINGLIST *blCreateSEQRES(PDB *pdb);
void blRenumResiduesPDB(PDB *pdb, int offset);
void blReplacePDBHeader(WHOLEPDB *wpdb, char *recordType,
                        STRINGLIST *replacement);
char *blFixSequence(char *seqresSequence, char *atomSequence,
                    char **seqresChains, char **atomChains,
                    char **outchains, BOOL IgnoreSEQRES, int nAtomChains,
                    BOOL upper, BOOL quiet, char *label);
#endif
