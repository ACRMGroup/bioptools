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
         bioplib/IndexPDB.o        \
         bioplib/OpenStdFiles.o    \
         bioplib/StoreString.o     \
         bioplib/FindResidue.o     \
         bioplib/FindNextResidue.o \
         bioplib/WritePDB.o        \
         bioplib/array2.o          \
         bioplib/GetWord.o

pdblistss : pdblistss.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdblistss.o $(LFILES)

distclean : clean
	\rm -f pdblistss
