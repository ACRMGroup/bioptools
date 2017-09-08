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
         bioplib/array2.o \
         bioplib/HAddPDB.o \
         bioplib/atomtype.o \
         bioplib/hbond.o \
         bioplib/FindHetatmResidue.o \
         bioplib/FindNextChainPDB.o \
         bioplib/angle.o \
         bioplib/simpleangle.o \
         bioplib/RenumAtomsPDB.o \
	 bioplib/StructurePDB.o 


pdbhbond : pdbhbond.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdbhbond.o $(LFILES)

distclean : clean
	\rm -f pdbhbond
