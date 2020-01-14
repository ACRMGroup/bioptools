DEST   = $(HOME)/bin
XML    = # -DXML_SUPPORT $(shell xml2-config --cflags)
XMLLIB = # -lxml2
#
GUNZIP = -DGUNZIP_SUPPORT
#
CC     = cc
COPT   = $(XML) $(GUNZIP)
LIBS   = $(XMLLIB) -lm
LFILES = bioplib/BuildConect.o     \
         bioplib/FreeStringList.o  \
         bioplib/padterm.o         \
         bioplib/chindex.o         \
         bioplib/fsscanf.o         \
         bioplib/ParseRes.o        \
         bioplib/ReadPDB.o         \
         bioplib/WritePDB.o        \
         bioplib/IndexPDB.o        \
         bioplib/OpenStdFiles.o    \
         bioplib/FindResidue.o     \
         bioplib/FindResidueSpec.o \
         bioplib/FindNextResidue.o \
         bioplib/StructurePDB.o    \
         bioplib/CopyPDB.o         \
         bioplib/phi.o             \
         bioplib/StoreString.o     \
         bioplib/OpenFile.o        \
         bioplib/openorpipe.o      \
         bioplib/secstr.o          \
         bioplib/array3.o          \
         bioplib/array2.o 	   \
         bioplib/hash.o 	   \
         bioplib/prime.o 	   \
         bioplib/stringutil.o 	   \
         bioplib/GetPDBChainLabels.o \
         bioplib/PDBHeaderInfo.o   \
         bioplib/throne.o	   \
         bioplib/strcatalloc.o	   \
         bioplib/GetWord.o	   \
         bioplib/stringcat.o	   \



pdbsecstr : pdbsecstr.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdbsecstr.o $(LFILES)

distclean : clean
	\rm -f pdbsecstr

install :
	mkdir -p $(DEST)
	cp pdbsecstr $(DEST)
