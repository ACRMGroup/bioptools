#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    makemake
#   File:       makemake.pl
#   
#   Version:    V1.4
#   Date:       21.11.14
#   Function:   Build the Makefile for BiopTools
#   
#   Copyright:  (c) Dr. Andrew C. R. Martin, UCL, 2014
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Institute of Structural and Molecular Biology
#               Division of Biosciences
#               University College
#               Gower Street
#               London
#               WC1E 6BT
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying LICENCE file
#
#*************************************************************************
#
#   Description:
#   ============
#   Builds the Makefile for building BiopTools
#
#*************************************************************************
#
#   Usage:
#   ======
#   ./makemake.pl [-bioplib] [-prefix=xxx] [-bindir=xxx] [-datadir=xxx]
#                 [-libdir=xxx] [-incdir=xxx]
#   -bioplib     - Build a local version of BiopLib
#   -prefix=xxx  - Change the prefix in front of all directories from
#                  $HOME to xxx
#   -bindir=xxx  - Change the installation binary directory to xxx
#   -datadir=xxx - Change the installation data directory to xxx
#   -libdir=xxx  - Change the location of the BiopLib libraries to xxx
#   -incdir=xxx  - Change the location of the BiopLib include files to xxx
#
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0    07.11.14  Original   By: ACRM
#   V1.1    07.11.14  Adds the -bioplib handling
#   V1.2    07.11.14  Writes the data directory in install
#   V1.3    07.11.14  Grabs the bioplib submodule if this is unpacked
#                     using git clone. If obtained as an archive then
#                     this grabs the bioplib archive as well
#   V1.4    21.11.14  Added -install option
#
#*************************************************************************
$::biopversion = "3.2";
#*************************************************************************
$::biopgit     = "https://github.com/ACRMGroup/bioplib/archive/V";
$::biopext     = ".tar.gz";
#*************************************************************************
# Deal with the command line
UsageDie() if(defined($::h) || defined($::help));
$::bioplib = 0                 if(!defined($::bioplib));
$::prefix  = $ENV{'HOME'}      if(!defined($::prefix));
$::install = $::prefix         if(!defined($::install));
$::bindir  = "$::install/bin"  if(!defined($::bindir));
$::datadir = "$::install/data" if(!defined($::datadir));
if($::bioplib)
{
    $::libdir  = "./bioplib"     if(!defined($::libdir));
    $::incdir  = "./bioplib"     if(!defined($::incdir));
}
else
{
    $::libdir  = "$::prefix/lib"     if(!defined($::libdir));
    $::incdir  = "$::prefix/include" if(!defined($::incdir));
}

# Main program
GetBiopLib()        if($::bioplib);
my @cFiles = GetCFileList('.');
my @exeFiles = StripExtension(@cFiles);
open(my $makefp, ">Makefile") || die "Can't open Makefile for writing";
WriteFlags($makefp, $::libdir, $::incdir, $::bindir, $::datadir);
WriteTargets($makefp, @exeFiles);
WriteDummyRule($makefp, $::bioplib);
WriteInstallRule($makefp, @exeFiles);
WriteCleanRule($makefp, $::bioplib, @exeFiles);
WriteLinksRule($makefp);
foreach my $cFile (@cFiles)
{
    WriteRule($makefp, $cFile);
}

close $makefp;

#*************************************************************************
# Uses git submodule to obtain the BiopLib code if you have downloaded
# this as a git 
#
# 07.11.14 Original   By: ACRM
sub GetBiopLib
{
    if(-d "../.git")
    {
        system("(cd ..; git submodule init)");
        system("(cd ..; git submodule update)");
    }
    else
    {
        if((! -d "libsrc/bioplib/.git") && (! -d "libsrc/bioplib/src"))
        {
            system("\\rm -rf libsrc/bioplib");
            my $url = "$::biopgit$::biopversion$::biopext";
            my $tar = "V$::biopversion$::biopext";
            system("(cd libsrc; wget $url; tar xvf $tar; mv bioplib-$::biopversion bioplib)");
        }
    }
}
#*************************************************************************
# Writes the rule for installing code in $BINDIR
#
# 06.11.14 Original   By: ACRM
sub WriteInstallRule
{
    my($makefp, @exeFiles) = @_;
    print $makefp <<__EOF;

install : 
\tmkdir -p \$(BINDIR)
\tcp \$(TARGETS) \$(BINDIR)
\tmkdir -p \$(DATADIR)
\tif [ ! -f \$(DATADIR)/radii.dat ] ; then cp ../data/radii.dat \$(DATADIR) ; fi
\t\@echo " "
\t\@echo " --- INSTALL COMPLETE --- "
\t\@echo " "
\t\@echo "To use some programs, you will need to set the environment variable"
\t\@echo "DATADIR to \$(DATADIR)"

__EOF
}

#*************************************************************************
# Writes te rule for creating links from new program names to old names
#
# 06.11.14 Original   By: ACRM
sub WriteLinksRule
{
    my($makefp, @exeFiles) = @_;
    print $makefp <<__EOF;

links : 
\t(cd \$(BINDIR); \\rm -f addhet      ; ln -s pdbaddhet          addhet       )
\t(cd \$(BINDIR); \\rm -f as2bval     ; ln -s naccess2bval       as2bval      )
\t(cd \$(BINDIR); \\rm -f atomcount   ; ln -s pdbatomcount       atomcount    )
\t(cd \$(BINDIR); \\rm -f atomsel     ; ln -s pdbatomsel         atomsel      )
\t(cd \$(BINDIR); \\rm -f avbr        ; ln -s pdbavbr            avbr         )
\t(cd \$(BINDIR); \\rm -f centralres  ; ln -s pdbcentralres      centralres   )
\t(cd \$(BINDIR); \\rm -f chainpdb    ; ln -s pdbchain           chainpdb     )
\t(cd \$(BINDIR); \\rm -f checkforres ; ln -s pdbcheckforres     checkforres  )
\t(cd \$(BINDIR); \\rm -f countpdb    ; ln -s pdbcount           countpdb     )
\t(cd \$(BINDIR); \\rm -f dummystrip  ; ln -s pdbdummystrip      dummystrip   )
\t(cd \$(BINDIR); \\rm -f findresrange; ln -s pdbfindresrange    findresrange )
\t(cd \$(BINDIR); \\rm -f flip        ; ln -s pdbflip            flip         )
\t(cd \$(BINDIR); \\rm -f getchain    ; ln -s pdbgetchain        getchain     )
\t(cd \$(BINDIR); \\rm -f getpdb      ; ln -s pdbgetzone         getpdb       )
\t(cd \$(BINDIR); \\rm -f getresidues ; ln -s pdbgetresidues     getresidues  )
\t(cd \$(BINDIR); \\rm -f hetstrip    ; ln -s pdbhetstrip        hetstrip     )
\t(cd \$(BINDIR); \\rm -f hstrip      ; ln -s pdbhstrip          hstrip       )
\t(cd \$(BINDIR); \\rm -f makepatch   ; ln -s pdbmakepatch       makepatch    )
\t(cd \$(BINDIR); \\rm -f nullstrip   ; ln -s pdbdummystrip      nullstrip    )
\t(cd \$(BINDIR); \\rm -f numpdb      ; ln -s setpdbnumbering    numpdb       )
\t(cd \$(BINDIR); \\rm -f patchbval   ; ln -s pdbpatchbval       patchbval    )
\t(cd \$(BINDIR); \\rm -f patchpdbnum ; ln -s pdbpatchnumbering  patchpdbnum  )
\t(cd \$(BINDIR); \\rm -f renumpdb    ; ln -s pdbrenum           renumpdb     )
\t(cd \$(BINDIR); \\rm -f rmspdb      ; ln -s pdbcalcrms         rmspdb       )
\t(cd \$(BINDIR); \\rm -f rotate      ; ln -s pdbrotate          rotate       )
\t(cd \$(BINDIR); \\rm -f solv        ; ln -s pdbsolv            solv         )
\t(cd \$(BINDIR); \\rm -f splitchains ; ln -s pdbsplitchains     splitchains  )
\t(cd \$(BINDIR); \\rm -f sumbval     ; ln -s pdbsumbval         sumbval      )
\t(cd \$(BINDIR); \\rm -f torsions    ; ln -s pdbtorsions        torsions     )
\t(cd \$(BINDIR); \\rm -f transpdb    ; ln -s pdbtranslate       transpdb     )

__EOF
}

#*************************************************************************
# Writes the rule to remove the executables from the source directory
#
# 06.11.14 Original   By: ACRM
sub WriteCleanRule
{
    my($makefp, $bioplib, @exeFiles) = @_;
    if($bioplib)
    {
        print $makefp <<__EOF;

clean : 
\t\\rm -rf bioplib
\t(cd libsrc/bioplib/src; make clean)
\t\\rm -f \$(TARGETS)

__EOF
    }
    else
    {
        print $makefp <<__EOF;

clean : 
\t\\rm -f \$(TARGETS)

__EOF
    }
}

#*************************************************************************
# Writes a rule to build an executable from a C file
#
# 06.11.14 Original   By: ACRM
sub WriteRule
{
    my($makefp, $cFile) = @_;
    my $exeFile = $cFile;
    $exeFile =~ s/\.c$//;
    print $makefp <<__EOF;

$exeFile : $cFile
\t\$(CC) \$(CFLAGS) -o \$@ \$< \$(LFLAGS)
__EOF

}

#*************************************************************************
# Writes the dummy rule for building everything
#
# 06.11.14 Original   By: ACRM
sub WriteDummyRule
{
    my($makefp, $bioplib) = @_;
    if($bioplib)
    {
        print $makefp <<__EOF;

all : bioplib \$(TARGETS)
\t\@echo " "
\t\@echo " --- BUILD COMPLETE --- "
\t\@echo " "
\t\@echo "To use some programs, you will need to set the environment variable"
\t\@echo "DATADIR to \$(DATADIR)"

bioplib :
\t(cd libsrc/bioplib/src; make)
\tmkdir -p bioplib
\tcp libsrc/bioplib/src/*.[ah] bioplib

__EOF
    }
    else
    {
        print $makefp <<__EOF;

all : \$(TARGETS)
\t\@echo " "
\t\@echo " --- BUILD COMPLETE --- "
\t\@echo " "
\t\@echo "To use some programs, you will need to set the environment variable"
\t\@echo "DATADIR to \$(DATADIR)"

__EOF
    }

}

#*************************************************************************
# Write the flags for the compiler and directories
#
# 06.11.14 Original   By: ACRM
sub WriteFlags
{
    my($makefp, $libdir, $incdir, $bindir, $datadir) = @_;
    print $makefp <<__EOF;
CC      = gcc
BINDIR  = $bindir
DATADIR = $datadir
CFLAGS  = -O3 -ansi -Wall -pedantic -I$incdir -L$libdir
LFLAGS  = -lbiop -lgen -lm -lxml2
__EOF
}

#*************************************************************************
# Write the list of targets
#
# 06.11.14 Original   By: ACRM
sub WriteTargets
{
    my ($makefp, @exeFiles) = @_;
    print $makefp "TARGETS = ";
    foreach my $exe (@exeFiles)
    {
        print $makefp "$exe ";
    }
    print $makefp "\n";
}

#*************************************************************************
# Build a list of target excutables by remove the extensions from the
# C source files
#
# 06.11.14 Original   By: ACRM
sub StripExtension
{
    my (@inputs) = @_;
    my @outputs = ();

    foreach my $input (@inputs)
    {
        my $output = $input;
        $output =~ s/\.c$//;
        push @outputs, $output;
    }

    return(@outputs);
}


#*************************************************************************
# Get a list of C source files
#
# 06.11.14 Original   By: ACRM
sub GetCFileList
{
    my($dir) = @_;
    my @files;
    if(opendir(my $dirfp, $dir))
    {
        @files = grep(/\.c$/, readdir $dirfp);
        closedir $dirfp;
    }
    else
    {
        print STDERR "Can't open current directory\n";
        exit 1;
    }
    return(@files);
}


#*************************************************************************
# Print a usage message and exit
#
# 07.11.14 Original   By: ACRM
sub UsageDie
{
    print <<__EOF;

makemake.pl (c) 2014 UCL, Dr. Andrew C.R. Martin

Usage: ./makemake.pl [-bioplib] [-prefix=xxx] [-bindir=xxx] [-datadir=xxx]
                     [-libdir=xxx] [-incdir=xxx]
                     
       -bioplib     - Build a local version of BiopLib

    The following options are shown in reverse order of priority - in
    other words the more specific options (such as -bindir) will take
    priority over global settings such as -prefix and -install
       -prefix=xxx  - Change the prefix in front of all directories from
                      \$HOME to xxx
       -install=xxx - Change the prefix in front of the installation 
                      directories from \$HOME (or whatever is specified
                      with -prefix) to xxx
       -bindir=xxx  - Change the installation binary directory to xxx
       -datadir=xxx - Change the installation data directory to xxx
       -libdir=xxx  - Change the location of the BiopLib libraries to xxx
       -incdir=xxx  - Change the location of the BiopLib include files to xxx

Builds the Makefile for building BiopTools

By default, this script assumes that you you have installed BiopLib in
\$HOME/lib with include files in \$HOME/include and you wish to
install the programs in \$HOME/bin  These defaults can be changed with
the command line switches.

__EOF

    exit 0;
}
