# makefile for building support libraries under UNIX:

# Install destination
DEST=${HOME}

# Destination in which to place the library and the include files
LIBDEST       = $(DEST)/lib
INCDEST       = $(DEST)/include
DATADEST      = $(DEST)/data
SHAREDLIBDEST = $(DEST)/lib

# Before stuff moves into bioplib, we need to know where the libraries 
# and includes are
DEVINC=-I $(INCDEST)

# Doxygen and adddoxy directories
DOCDIR      = ../doc
DOXYDIR     = $(DOCDIR)/doxygen
ADDDOXY     = $(DOCDIR)/adddoxy/adddoxy.pl
ADDDOXYCONF = $(DOCDIR)/adddoxy/adddoxy.conf

# Define the C compiler here
CC = gcc
 
# Define any options for the C compiler here
COPT = -fPIC -ansi -Wall -pedantic -g -DGUNZIP_SUPPORT -finline-functions $(DEVINC)

# Generate error messages for deprecated functions.
# Default is to display deprecated error message unless environment 
# variable BIOPLIB_DEPRECATED_QUIET is set.
# Uncomment compile options below to ignore the environment and either 
# always display the error message or silence the error message:
#COPT := $(COPT) -D BIOPLIB_DEPRECATED_CHECK
#COPT := $(COPT) -D BIOPLIB_DEPRECATED_QUIET

# Include libxml2
# Note: xml2-config is installed with libxml2
# When you compile code you need to link to libxml2 with -lxml2
# Comment out this line if you do not require PDBML (XML) support
COPT := $(COPT) -D XML_SUPPORT $(shell xml2-config --cflags)

# Use single letter check for filetype
# Only check first character of file when detecting file type (compressed
# file or pdbml).
# For systems were ungetc() fails after pushing back a single character to
# the input stream.
# Note: This option is required for MS Windows.
#COPT := $(COPT) -D SINGLE_CHAR_FILECHECK

# Define the archive/library command here
AR = ar r

# If ranlib required when creating libraries, set here; 
# otherwise set to echo
RANLIB = ranlib
#RANLIB = echo

###########################################################
#          DON'T CHANGE ANYTHING BELOW THIS POINT         #
###########################################################

# Version numbers for shared libraries - strictly the MINOR numbers here
# are the minor and release numbers
# libgen
GMAJOR = 0
GMINOR = 0.1
# libbiop
BMAJOR = 0
BMINOR = 0.1

# Files for libgenextras.a
OFILESG = 


# Files for libbiopextras.a
OFILESB = secstr.o 


# Static libraries - the default
all : libgenextras.a libbiopextras.a 

libgenextras.a : $(OFILESG)
	$(AR) libgenextras.a $? 
	$(RANLIB) libgenextras.a
libbiopextras.a : $(OFILESB)
	$(AR) libbiopextras.a $? 
	$(RANLIB) libbiopextras.a

# Shared libraries
shared : libgenextrass.so.$(GMAJOR).$(GMINOR) libbiopextrass.so.$(BMAJOR).$(BMINOR)

libgenextrass.so.$(GMAJOR).$(GMINOR) : $(OFILESG)
	$(CC) -shared -fPIC -Wl,-soname,libgenextrass.so.$(GMAJOR) -o libgenextrass.so.$(GMAJOR).$(GMINOR) $? -lc
libbiopextrass.so.$(BMAJOR).$(BMINOR) : $(OFILESB)
	$(CC) -shared -fPIC -Wl,-soname,libbiopextrass.so.$(BMAJOR) -o libbiopextrass.so.$(BMAJOR).$(BMINOR) $? -lc

# C compilation
.c.o : 
	$(CC) $(COPT) -o $@ -c $<

# Cleanup
clean :
	rm -f $(OFILESB) $(OFILESG) libgenextras.a libbiopextras.a libgenextrass.so.$(GMAJOR).$(GMINOR) libbiopextrass.so.$(BMAJOR).$(BMINOR)

# Documentation
doxygen :
	$(ADDDOXY) -conf=$(ADDDOXYCONF) *.c *.h
	mkdir -p $(DOXYDIR)/docsrc
	cp $(DOXYDIR)/docsrcinput/*.dox $(DOXYDIR)/docsrc
	(cd $(DOXYDIR); doxygen)
doxyclean :
	rm -rf $(DOCDIR)/html $(DOCDIR)/latex $(DOXYDIR)/docsrc


# Installation
install :
	mkdir -p $(LIBDEST)
	mkdir -p $(INCDEST)/bioplib
	cp libgenextras.a $(LIBDEST)
	cp libbiopextras.a $(LIBDEST)
	cp *.h $(INCDEST)/bioplib

installdata :
	mkdir -p $(DATADEST)
	(cd ../data; cp -Rcp * $(DATADEST))

installshared : 
	mkdir -p $(SHAREDLIBDEST)
	mkdir -p $(INCDEST)/bioplib
	cp *.h $(INCDEST)/bioplib
	cp libgenextrass.so.$(GMAJOR).$(GMINOR) $(SHAREDLIBDEST)
	cp libbiopextrass.so.$(GMAJOR).$(GMINOR) $(SHAREDLIBDEST)
	(cd $(SHAREDLIBDEST); ln -s libgenextrass.so.$(GMAJOR).$(GMINOR) libgenextrass.so.$(GMAJOR); ln -s libgenextrass.so.$(GMAJOR) libgenextrass.so)
	(cd $(SHAREDLIBDEST); ln -s libbiopextrass.so.$(BMAJOR).$(BMINOR) libbiopextrass.so.$(BMAJOR); ln -s libbiopextrass.so.$(BMAJOR) libbiopextrass.so)

# Help
help :
	@echo "make               : Build the static libraries"
	@echo "make install       : Install the static libraries"
	@echo "make installdata   : Install the data files"
	@echo "make shared        : Build shared libraries"
	@echo "make installshared : Install shared libraries"
	@echo "make clean         : Remove object and library files from compilation directory"
	@echo "make doxygen       : Generate doxygen documentation"
	@echo "make doxyclean     : Remove doxygen documentation"
