XML    = # -DXML_SUPPORT $(shell xml2-config --cflags)
XMLLIB = # -lxml2
#
GUNZIP = -DGUNZIP_SUPPORT
#
CC     = cc
COPT   = -O3 $(XML) $(GUNZIP)
LIBS   = $(XMLLIB) -lm
LFILES = bioplib/ReadPDB.o bioplib/FreeStringList.o bioplib/StoreString.o \
	 bioplib/fsscanf.o bioplib/chindex.o bioplib/FindNextResidue.o \
	 bioplib/padterm.o bioplib/WritePDB.o bioplib/FindResidue.o \
	 bioplib/BuildConect.o bioplib/OpenStdFiles.o bioplib/hash.o \
	 bioplib/CalcExtSD.o bioplib/prime.o bioplib/SelAtPDB.o \
	 bioplib/stringutil.o bioplib/CopyPDB.o bioplib/GetWord.o \
	 bioplib/array2.o bioplib/IndexPDB.o


distmat : distmat.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f distmat.o $(LFILES)

distclean : clean
	\rm -f distmat
