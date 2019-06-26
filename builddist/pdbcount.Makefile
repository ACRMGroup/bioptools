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
         bioplib/ReadPDB.o         \
         bioplib/WritePDB.o        \
         bioplib/chindex.o         \
	 bioplib/padterm.o         \
         bioplib/OpenStdFiles.o    \
         bioplib/FindResidue.o     \
         bioplib/FindNextResidue.o \
         bioplib/IndexPDB.o        \
         bioplib/StoreString.o     \
         bioplib/fsscanf.o         

pdbcount : pdbcount.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdbcount.o $(LFILES)

distclean : clean
	\rm -f pdbcount
