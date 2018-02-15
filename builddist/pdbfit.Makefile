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
         bioplib/array2.o          \
	 bioplib/SelAtPDB.o        \
         bioplib/fit.o             \
         bioplib/FitPDB.o          \
         bioplib/FitCaPDB.o        \
         bioplib/FitNCaCPDB.o      \
	 bioplib/fit.h             \
         bioplib/GetCGPDB.o        \
         bioplib/OriginPDB.o       \
         bioplib/TranslatePDB.o    \
         bioplib/ApMatPDB.o        \
         bioplib/matrix.h          \
         bioplib/MatMult3_33.o     \
         bioplib/GetPDBCoor.o      \
         bioplib/CalcRMSPDB.o

pdbfit : pdbfit.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f pdbfit.o $(LFILES)

distclean : clean
	\rm -f pdbfit
