pdb=test2.pdb

../../pdbfindnearres -l L50-L56 lys $pdb
# L45
../../pdbfindnearres -l L50 lys $pdb
# L53
../../pdbfindnearres -l L24-L34,L50-L56 lys $pdb
# L45
../../pdbfindnearres -l L24-L34 lys $pdb
# L53
