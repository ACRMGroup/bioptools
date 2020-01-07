XML    = # -DXML_SUPPORT $(shell xml2-config --cflags)
XMLLIB = # -lxml2
#
GUNZIP = -DGUNZIP_SUPPORT
#
CC     = cc
COPT   = -O3 $(XML) $(GUNZIP)
LIBS   = $(XMLLIB) -lm
LFILES = bioplib/ReadPDB.o	   \
         bioplib/fsscanf.o         \
         bioplib/padterm.o         \
         bioplib/FreeStringList.o  \
         bioplib/StoreString.o     \
         bioplib/chindex.o         \
         bioplib/FindNextResidue.o \
         bioplib/OpenStdFiles.o    \
         bioplib/ParseRes.o        \
         bioplib/WritePDB.o        \
         bioplib/FindResidue.o     \
         bioplib/BuildConect.o     \
         bioplib/IndexPDB.o

pdbpatchbval : pdbpatchbval.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdbpatchbval.o $(LFILES)

distclean : clean
	\rm -f pdbpatchbval

