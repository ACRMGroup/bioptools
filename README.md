bioptools
=========

Tools making use of BiopLib


See `INSTALL.md` for installation instructions!


chaincontacts
-------------
Calculates contacts between chains at the atom and residue level. You
can get all contacts or contacts between specified chains.

checkpdb
--------
Start of a program for performing various checks on a PDB file

naccess2bval
------------
Rewrites the output from naccess solvent accessibility as a standard PDB
format file with accessibility in the B-val column and radius in the
occupancy column.

pdb2ms
------
Converts a PDB file to input for the Connoly MS program

pdb2pdbml
---------
Converts a PDB file to PDBML format

pdbml2pdb
---------
Converts a PDBML file to PDB format

pdb2pir
-------
Extracts a PIR sequence file from a PDB file.

pdb2xyz
-------
Convert PDB format to GROMOS XYZ. N.B. Does NOT correct atom order.

pdbatomcount 
------------
Counts the number of atoms within the specified radius of each atom in
a PDB structure.

pdbatoms
--------
Extracts only the coordinates (i.e. ATOM and HETATM records) from a
PDB file, discarding the header and footer records

pdbatomsel
----------
Selects specfied atom types from a PDB file. Assumes C-alpha if no atoms
are specified.

pdbavbr
-------
Calculates means and standard deviations for B-values per residue.

pdbcalcrms
----------
Calculates an RMS between 2 PDB files. No fitting is performed.

pdbcentralres
-------------
Identifies the residue closest to the centroid of a protein.

pdbchain
--------
Splits a PDB file into chains using distance criteria.

pdbcheckforres
--------------
Checks whether a specified residue exists in a PDB file.

pdbconect
---------
Builds CONECT records for a PDB file. Deletes existing records.

pdbcount
--------
Counts chains, residues & atoms in a PDB file.

pdbcter
-------
Rename C-terminal oxygens in standard style and generate second one if required.

pdbdummystrip
-------------
Removes atoms from a PDB file which have NULL coordinates (i.e.
x = y = z = 9999.0)

pdbfindnearres
--------------
Takes one or more zones in a PDB file and outputs all residues of a 
specified type within a given distance of the zone(s).

pdbfindresrange
---------------
Takes a PDB file as input and given a key residue and a number of
residues will return the residue identifiers for residues this number
of positions before and width after the key residue.

pdbflip
-------
pdbflip is a rather crude and simple program for correcting the atom
naming of equivalent atoms about freely rotable bonds.

pdbgetchain
-----------
pdbgetchain reads a PDB file and write out only those chains specified
on the command line.

pdbgetresidues
--------------
Takes a list of residue specifications and extracts just those residues 
from a PDB file

pdbgetzone
----------
Extracts a specified zone from a PDB file

pdbhadd
-------
Add hydrogens to a PDB file.

pdbheader
---------
Prints the key header information from a PDB (title; molecule and
species for each chain, etc.)

pdbhetstrip
-----------
Removes het atoms from a PDB file.

pdbhphob
--------
Patches the residue hydrophobicity values into the B-value 
(temperature factor) column to allow colouring by hydrophobicity.

pdbhstrip
---------
Removes hydrogens from a PDB file.

pdbline
-------
Calculates the best fit line through a set of PDB coordinates.

pdblistss
---------
Calculates and lists the disulphides in a PDB file.

pdbmakepatch
------------
Generates a patch around a specified residue.

pdborder
--------
Correct atom order of a PDB file.

pdborigin
---------
Moves a set of PDB coordinates to the origin.

pdbpatchbval
------------
Takes a patch file containing residue specifications and
values one to a line and patches the B-value (or occupancy) for that
residue with the specified values.

pdbpatchnumbering
-----------------
Patches the numbering of a PDB file from a patch file containing
residue numbers.

pdbrenum
--------
Renumber a PDB file

pdbrotate
---------
Rotates a PDB file

pdbsecstr
---------
Calculate secondary structure for a PDB file

pdbselect
---------
Allows the extraction of different occupancies and different models
from a PDB file

pdbsolv
-------
Performs solvent accessibility calculations according to the method of
Lee and Richards. 

pdbsphere
---------
pdbsphere identifies residues within a specified radius of a specified
residue. 

pdbsplitchains
--------------
Split a PDB file into separate files for each chain

pdbsumbval
----------
Sums the b-values over each residue and places the summed values in the
b-value column. Can also calculate averages.

pdbsymm
-------
Applies non-cystollographic symmetry operations specified in REMARK 350
records to a PDB file. Note that this program has certain limitations.
Run pdbsymm -h for details.

pdbtorsions
-----------
Generates a set of backbone torsions from a PDB file.

pdbtranslate
------------
Translates a PDB file

rangecontacts
-------------
Performs a contact analysis at the residue level to find residues whose
sidechains contact any atom of the residues specified in a given range

scorecons
---------
A program similar to (and which predates) Will Valdar's scorecons program
for scoring conservation in a sequence alignment.

setpdbnumbering
---------------
Applies a standard numbering scheme to a set of PDB files.
