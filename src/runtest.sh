make; make install
(cd pdbsecstr; source ./build.sh)
(cd pdbsecstr; ./pdbsecstr test/pdb1yqv.ent > test/1yqv.ss2; diff test/1yqv.ss test/1yqv.ss2)
