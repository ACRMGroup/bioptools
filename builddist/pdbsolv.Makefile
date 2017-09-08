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
         bioplib/FindResidueSpec.o \
         bioplib/FindNextResidue.o \
         bioplib/CopyPDB.o         \
         bioplib/StoreString.o     \
         bioplib/OpenFile.o        \
         bioplib/openorpipe.o      \
         bioplib/StripWatersPDB.o  \
         bioplib/access.o          \
         bioplib/array2.o 


pdbsolv : pdbsolv.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdbsolv.o $(LFILES)

distclean : clean
	\rm -f pdbsolv
