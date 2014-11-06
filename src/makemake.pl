#!/usr/bin/perl

my @cFiles = GetCFileList('.');
my @exeFiles = StripExtension(@cFiles);
open(my $makefp, ">Makefile") || die "Can't open Makefile for writing";
WriteFlags($makefp);
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

sub WriteInstallRule
{
    my($makefp, @exeFiles) = @_;
    print $makefp <<__EOF;

install : 
\tcp \$(TARGETS) \$(BINDIR)

__EOF
}

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

sub WriteCleanRule
{
    my($makefp, @exeFiles) = @_;
    print $makefp <<__EOF;

clean : 
\t\\rm -f \$(TARGETS)

__EOF
}

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

sub WriteDummyRule
{
    my($makefp) = @_;
    print $makefp <<__EOF;

all : \$(TARGETS)

__EOF
}

sub WriteFlags
{
    my($makefp) = @_;
    print $makefp <<__EOF;
CC     = gcc
BINDIR = \$(HOME)/bin
CFLAGS = -O3 -ansi -Wall -pedantic -I\$(HOME)/include -L\$(HOME)/lib
LFLAGS = -lbiop -lgen -lm -lxml2
__EOF
}

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
