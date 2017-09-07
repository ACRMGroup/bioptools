XML    = # -DXML_SUPPORT $(shell xml2-config --cflags)
XMLLIB = # -lxml2
#
GUNZIP = -DGUNZIP_SUPPORT
#
CC     = cc
COPT   = -O3 $(XML) $(GUNZIP)
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
         bioplib/FindNextResidue.o \
         bioplib/StripHPDB.o       \
         bioplib/CopyPDB.o         \
         bioplib/StoreString.o     \
         bioplib/ExtractZonePDB.o  \
         bioplib/DupePDB.o         \
         bioplib/filename.o        \
         bioplib/GetWord.o         \
         bioplib/array2.o          \
         bioplib/phi.o             \
         bioplib/throne.o          \
         bioplib/FindAtomInRes.o   \
         bioplib/SelAtPDB.o        \
         bioplib/StructurePDB.o


pdbtorsions : pdbtorsions.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdbtorsions.o $(LFILES)

distclean : clean
	\rm -f pdbtorsions
