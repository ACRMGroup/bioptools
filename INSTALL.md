BiopTools - Installation
========================

If you have downloaded and installed BiopLib already
----------------------------------------------------

        cd src
        ./makemake.pl
        make
        make install

This assumes that you have installed BiopLib in $HOME/lib with include
files in $HOME/include and you wish to install the programs in
$HOME/bin

You may change these default locations with

-   -prefix=xxxx - Change the prefix in front of all directories from
                   $HOME to xxxx
-   -bindir=xxxx - Change the installation binary directory to xxxx
-   -libdir=xxxx - Change the location of the BiopLib libraries to xxxx
-   -incdir=xxxx - Change the location of the BiopLib include files to xxxx

If you have been a previous user of tools distributed by us, then you
will find that many of the program names have changed. To create
symbolic links in the binary directory to the old names, do:

        make links



If you have NOT downloaded and installed BiopLib already
--------------------------------------------------------

        cd src
        ./makemake.pl -bioplib
        make
        make install

This assumes that you wish to install the programs in $HOME/bin

You may change this default location with

-   -prefix=xxxx - Change the prefix in front of all directories from
                   $HOME to xxxx
-   -bindir=xxxx - Change the installation binary directory to xxxx

If you have been a previous user of tools distributed by us, then you
will find that many of the program names have changed. To create
symbolic links in the binary directory to the old names, do:

        make links



