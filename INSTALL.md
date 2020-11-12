BiopTools - Installation
========================

There are two modes of installation that depend on whether you have
already installed our BiopLib library.

If you have NOT downloaded and installed BiopLib already
--------------------------------------------------------

**Note!!! This will only work if you have downloaded one of the
  *releases* of BiopTools rather than using `git clone`.**

e.g.

```
wget https://github.com/ACRMGroup/bioptools/archive/V1.9.tar.gz
```
(This is because the latest code obtained with `git clone` may use
BiopLib routines not contained in the version of BiopLib that will be
downloaded by makemake.pl.)

### Ensure you have libxml2 installed

This will normally be already installed and available on Linux
systems. If not then it is installed on Fedora/CentOS systems using
(as root):

```
sudo yum -y install libxml2 libxml2-devel
```

or on Debian/Ubuntu systems using:

```
sudo apt-get install libxml2 libxml2-dev
```

On other systems, you will need to install libxml2 manually from 
http://xmlsoft.org/downloads.html


### Ensure you have 'git' installed

```
sudo yum -y install git
```

or on Debian/Ubuntu systems using:

```
sudo apt-get install git
```


### Now compile and install bioptools

```
cd src
./makemake.pl -bioplib
make
make install
```

This assumes that you wish to install the programs in `$HOME/bin`

You may change this default location with the following options for
the makemake.pl script

- `-prefix=xxxx`  - Change the prefix in front of all directories from
                    $HOME to xxxx
- `-install=xxxx` - Change the prefix in front of the installation 
                    directories from $HOME (or whatever was specified
                    with -prefix) to xxxx. Has the same effect as 
                    -prefix, but will override whatever was given 
                    with -prefix
- `-bindir=xxxx`  - Change the installation binary directory to xxxx
- `-datadir=xxxx` - Change the installation data directory to xxxx

If you have been a previous user of tools distributed by us, then you
will find that many of the program names have changed. To create
symbolic links in the binary directory to the old names, do:

```
make links
```

If you have downloaded and installed BiopLib already
----------------------------------------------------

*(Note you will need to have installed libxml2 - see above - to do
anything with BiopLib that involves PDB files.)*

```
cd src
./makemake.pl
make
make install
```

This assumes that you have installed BiopLib in $HOME/lib with include
files in $HOME/include and you wish to install the programs in
$HOME/bin

You may change these default locations with the following options for
the makemake.pl script

- `-prefix=xxxx`  - Change the prefix in front of all directories from
                    $HOME to xxxx
- `-install=xxxx` - Change the prefix in front of the installation 
                    directories from $HOME (or whatever was specified
                    with -prefix) to xxxx
- `-libdir=xxxx`  - Change the location of the BiopLib libraries to xxxx
- `-incdir=xxxx`  - Change the location of the BiopLib include files to xxxx
- `-bindir=xxxx`  - Change the installation binary directory to xxxx
- `-datadir=xxxx` - Change the installation data directory to xxxx

If you have been a previous user of tools distributed by us, then you
will find that many of the program names have changed. To create
symbolic links in the binary directory to the old names, do:

```
make links
```


Configuration to access data files
----------------------------------

You now need to create an environment variable, `DATADIR`, so point to
where the data files are installed.  Assuming you have used the
default install, the programs are installed in `$HOME/bin/` and the
data files are installed in `$HOME/data/`.

**If using the sh/BASH shell**, (the default in Linux), then add the
following line to your `~/.bashrc` file:

```
export DATADIR=$HOME/data
```

**If using the csh/tcsh shell**, then add the following line to your
`~/.cshrc` file:

```
setenv DATADIR $HOME/data
```

Don't forget that this will not become active until you have started a
new shell or typed the appropriate `setenv` or `export` command at the
command line.

If you installed the software elsewhere (with the `-prefix`,
`-install` or `-datadir`, then you need to alter these appropriately.

