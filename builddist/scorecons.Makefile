XML    = # -DXML_SUPPORT $(shell xml2-config --cflags)
XMLLIB = # -lxml2
#
GUNZIP = -DGUNZIP_SUPPORT
#
CC     = cc
COPT   = -O3 $(XML) $(GUNZIP)
LIBS   = $(XMLLIB) -lm
LFILES = bioplib/OpenStdFiles.o bioplib/ReadPIR.o bioplib/array2.o \
	bioplib/padchar.o bioplib/align.o bioplib/GetWord.o \
	bioplib/OpenFile.o



scorecons : scorecons.o $(LFILES)
	$(CC) $(COPT) -o $@ $< $(LFILES) $(LIBS)

.c.o :
	$(CC) $(COPT) -c -o $@ $<

clean :
	\rm -f scorecons.o $(LFILES)

distclean : clean
	\rm -f scorecons
