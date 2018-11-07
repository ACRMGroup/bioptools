PROFILE  = 
OPTDEBUG = -O3

EXE      = pdbline
CC       = gcc
OFILES   = pdbline.o
LFILES   = bioplib/ReadPDB.o bioplib/regression.o bioplib/eigen.o \
           bioplib/fsscanf.o bioplib/array2.o bioplib/ParseRes.o \
           bioplib/FindNextResidue.o bioplib/OpenStdFiles.o \
           bioplib/WritePDB.o bioplib/SelAtPDB.o bioplib/ExtractZonePDB.o \
           bioplib/chindex.o bioplib/padterm.o bioplib/CopyPDB.o \
           bioplib/DupePDB.o bioplib/BuildConect.o bioplib/FindResidue.o \
	   bioplib/FreeStringList.o bioplib/StoreString.o bioplib/IndexPDB.o
LOPTS    = $(PROFILE) -L $(HOME)/lib
LIBS     = -lm -lxml2 
COPTS    = $(PROFILE) $(OPTDEBUG) -ansi -Wall -pedantic -I $(HOME)/include/ -DNODEPRECATION -Wno-unused-function

$(EXE) : $(OFILES) $(LFILES)
	$(CC) $(LOPTS) -o $@ $(OFILES) $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPTS) -o $@ -c $<

clean :
	\rm -f $(OFILES) $(LFILES)

test :
	./$(EXE) A12 A25 TEST/1AFV_1.pdb TEST/test.pdb
	diff TEST/regressionLine.pdb TEST/test.pdb
