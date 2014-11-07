#!/usr/bin/perl -s
#*************************************************************************
#
#   Program:    makemake
#   File:       makemake.pl
#   
#   Version:    V1.0
#   Date:       07.11.14
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
#   ./makemake.pl [-prefix=xxxx] [-bindir=xxxx] 
#                 [-libdir=xxxx] [-incdir=xxxx]
#   -prefix=xxxx - Change the prefix in front of all directories from
#                  $HOME to xxxx
#   -bindir=xxxx - Change the installation binary directory to xxxx
#   -libdir=xxxx - Change the location of the BiopLib libraries to xxxx
#   -incdir=xxxx - Change the location of the BiopLib include files to xxxx
#
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0    07.11.14  Original   By: ACRM
#
#*************************************************************************
# Add the path of the executable to the library path
use FindBin;
use lib $FindBin::Bin;
# Or if we have a bin directory and a lib directory
#use Cwd qw(abs_path);
#use FindBin;
#use lib abs_path("$FindBin::Bin/../lib");


# Deal with the command line
UsageDie() if(defined($::h) || defined($::help));
$::prefix = $ENV{'HOME'} if(!defined($::prefix));
$::libdir = "$::prefix/lib" if(!defined($::libdir));
$::incdir = "$::prefix/include" if(!defined($::incdir));
$::bindir = "$::prefix/bin" if(!defined($::bindir));


my @cFiles = GetCFileList('.');
my @exeFiles = StripExtension(@cFiles);
open(my $makefp, ">Makefile") || die "Can't open Makefile for writing";
WriteFlags($makefp, $::libdir, $::incdir, $::bindir);
WriteTargets($makefp, @exeFiles);
WriteDummyRule($makefp);
WriteInstallRule($makefp, @exeFiles);
WriteCleanRule($makefp, @exeFiles);
WriteLinksRule($makefp);
foreach my $cFile (@cFiles)
{
    WriteRule($makefp, $cFile);
}

close $makefp;

#*************************************************************************
# Writes the rule for installing code in $BINDIR
#
# 06.11.14 Original   By: ACRM
sub WriteInstallRule
{
    my($makefp, @exeFiles) = @_;
    print $makefp <<__EOF;

install : 
\tcp \$(TARGETS) \$(BINDIR)

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
    my($makefp, @exeFiles) = @_;
    print $makefp <<__EOF;

clean : 
\t\\rm -f \$(TARGETS)

__EOF
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
    my($makefp) = @_;
    print $makefp <<__EOF;

all : \$(TARGETS)

__EOF
}

#*************************************************************************
# Write the flags for the compiler and directories
#
# 06.11.14 Original   By: ACRM
sub WriteFlags
{
    my($makefp, $libdir, $incdir, $bindir) = @_;
    print $makefp <<__EOF;
CC     = gcc
BINDIR = $bindir;
CFLAGS = -O3 -ansi -Wall -pedantic -I$incdir -L$libdir
LFLAGS = -lbiop -lgen -lm -lxml2
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

Usage: ./makemake.pl [-prefix=xxxx] [-bindir=xxxx] 
                     [-libdir=xxxx] [-incdir=xxxx]
       -prefix=xxxx - Change the prefix in front of all directories from
                      \$HOME to xxxx
       -bindir=xxxx - Change the installation binary directory to xxxx
       -libdir=xxxx - Change the location of the BiopLib libraries to xxxx
       -incdir=xxxx - Change the location of the BiopLib include files to xxxx

Builds the Makefile for building BiopTools

By default, this script assumes that you you have installed BiopLib in
\$HOME/lib with include files in \$HOME/include and you wish to
install the programs in \$HOME/bin  These defaults can be changed with
the command line switches.

__EOF

    exit 0;
}
