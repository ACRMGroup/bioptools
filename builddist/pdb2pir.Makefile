XML    = # -DXML_SUPPORT $(shell xml2-config --cflags)
XMLLIB = # -lxml2
#
GUNZIP = -DGUNZIP_SUPPORT
#
CC     = cc
COPT   = $(XML) $(GUNZIP)
LIBS   = $(XMLLIB) -lm
LFILES = bioplib/ReadPDB.o bioplib/throne.o bioplib/fsscanf.o \
	 bioplib/padterm.o bioplib/FreeStringList.c bioplib/StoreString.c \
	 bioplib/FindNextResidue.o bioplib/chindex.o \
	 bioplib/align.o bioplib/PDB2Seq.o bioplib/countchar.o  \
	 bioplib/hash.o bioplib/prime.o bioplib/GetWord.o \
	 bioplib/array2.o bioplib/stringutil.o bioplib/OpenFile.o \
	 bioplib/WritePDB.o bioplib/PDBHeaderInfo.o \
	 bioplib/strcatalloc.o bioplib/stringcat.o \
	 bioplib/GetPDBChainLabels.o bioplib/WritePIR.o bioplib/ParseRes.o \
	 bioplib/FindResidue.o bioplib/BuildConect.o bioplib/IndexPDB.o

pdb2pir : pdb2pir.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdb2pir.o $(LFILES)

distclean : clean
	\rm -f pdb2pir

