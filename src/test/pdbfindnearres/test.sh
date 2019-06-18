pdb=test.pdb

../../pdbfindnearres -l L50-L56 lys $pdb
# Expect L45
../../pdbfindnearres -l L50 lys $pdb
# Expect L53
../../pdbfindnearres -l L24-L34,L50-L56 lys $pdb
# Expect L45
../../pdbfindnearres -l L24-L34 lys $pdb
# Expect L53
