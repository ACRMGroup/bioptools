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
         bioplib/ReadPDB.o         \
         bioplib/WritePDB.o        \
         bioplib/IndexPDB.o        \
         bioplib/OpenStdFiles.o    \
         bioplib/FindResidue.o     \
         bioplib/FindNextResidue.o \
         bioplib/StoreString.o     \
         bioplib/GetWord.o         \
         bioplib/FixCterPDB.o      \
         bioplib/CalcCterCoords.o  \
         bioplib/CopyPDB.o         \
         bioplib/RenumAtomsPDB.o   \
         bioplib/array2.o 


pdbcter : pdbcter.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdbcter.o $(LFILES)

distclean : clean
	\rm -f pdbcter
