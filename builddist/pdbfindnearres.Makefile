XML    = # -DXML_SUPPORT $(shell xml2-config --cflags)
XMLLIB = # -lxml2
#
GUNZIP = -DGUNZIP_SUPPORT
#
CC     = cc
COPT   = -O3 $(XML) $(GUNZIP)
LIBS   = $(XMLLIB) -lm
LFILES = bioplib/fsscanf.o         \
         bioplib/ReadPDB.o         \
         bioplib/array2.o          \
         bioplib/InPDBZoneSpec.o   \
         bioplib/FindNextResidue.o \
         bioplib/OpenStdFiles.o    \
         bioplib/padterm.o         \
         bioplib/GetWord.o         \
         bioplib/FreeStringList.o  \
         bioplib/ParseRes.o        \
         bioplib/chindex.o         \
         bioplib/IndexPDB.o        \
         bioplib/BuildConect.o     \
         bioplib/FindResidue.o     \
         bioplib/InPDBZone.o       \
         bioplib/StoreString.o     \
         bioplib/WritePDB.o        


pdbfindnearres : pdbfindnearres.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdbfindnearres.o $(LFILES)

distclean : clean
	\rm -f pdbfindnearres
