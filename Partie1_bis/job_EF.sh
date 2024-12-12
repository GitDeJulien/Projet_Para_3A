#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH -C miriel

# Get input parameter

# Load dependencies
#module load compiler/g++/9.3.0
#module load mpi/openmpi/3.1.4-all
#module load compiler/gfortran/*

# cd in the project folder
#cd $HOME/Documents/Partie1

#Génère le mesh de lecture pour metis
cat build/meshprogc.data metis/proc1/dualformetis.data.dgraph.part.1 > metis/proc1/mesh_for_progc.data
cat build/meshprogc.data metis/proc2/dualformetis.data.dgraph.part.2 > metis/proc2/mesh_for_progc.data
cat build/meshprogc.data metis/proc3/dualformetis.data.dgraph.part.3 > metis/proc3/mesh_for_progc.data
cat build/meshprogc.data metis/proc4/dualformetis.data.dgraph.part.4 > metis/proc4/mesh_for_progc.data
cat build/meshprogc.data metis/proc5/dualformetis.data.dgraph.part.5 > metis/proc5/mesh_for_progc.data
cat build/meshprogc.data metis/proc10/dualformetis.data.dgraph.part.10 > metis/proc10/mesh_for_progc.data
cat build/meshprogc.data metis/proc15/dualformetis.data.dgraph.part.15 > metis/proc15/mesh_for_progc.data
cat build/meshprogc.data metis/proc20/dualformetis.data.dgraph.part.20 > metis/proc20/mesh_for_progc.data
cat build/meshprogc.data metis/proc24/dualformetis.data.dgraph.part.24 > metis/proc24/mesh_for_progc.data

gcc -o fromscotch.exe fromscotch.c

#Génère le mesh de lecture pour scotch
./fromscotch.exe build/meshprogc.data scotch/proc1/dualforscotch.map scotch/proc1/mesh_for_progc.data
./fromscotch.exe build/meshprogc.data scotch/proc2/dualforscotch.map scotch/proc2/mesh_for_progc.data
./fromscotch.exe build/meshprogc.data scotch/proc3/dualforscotch.map scotch/proc3/mesh_for_progc.data
./fromscotch.exe build/meshprogc.data scotch/proc4/dualforscotch.map scotch/proc4/mesh_for_progc.data
./fromscotch.exe build/meshprogc.data scotch/proc5/dualforscotch.map scotch/proc5/mesh_for_progc.data
./fromscotch.exe build/meshprogc.data scotch/proc10/dualforscotch.map scotch/proc10/mesh_for_progc.data
./fromscotch.exe build/meshprogc.data scotch/proc15/dualforscotch.map scotch/proc15/mesh_for_progc.data
./fromscotch.exe build/meshprogc.data scotch/proc20/dualforscotch.map scotch/proc20/mesh_for_progc.data
./fromscotch.exe build/meshprogc.data scotch/proc24/dualforscotch.map scotch/proc24/mesh_for_progc.data

gcc -o Preprocess.exe Preprocess.c

#Génère les fichier Data**.In metis
./Preprocess.exe 1 metis/proc1/mesh_for_progc.data metis/proc1/
./Preprocess.exe 2 metis/proc2/mesh_for_progc.data metis/proc2/
./Preprocess.exe 3 metis/proc3/mesh_for_progc.data metis/proc3/
./Preprocess.exe 4 metis/proc4/mesh_for_progc.data metis/proc4/
./Preprocess.exe 5 metis/proc5/mesh_for_progc.data metis/proc5/
./Preprocess.exe 10 metis/proc10/mesh_for_progc.data metis/proc10/
./Preprocess.exe 15 metis/proc15/mesh_for_progc.data metis/proc15/
./Preprocess.exe 20 metis/proc20/mesh_for_progc.data metis/proc20/
./Preprocess.exe 24 metis/proc24/mesh_for_progc.data metis/proc24/

./Preprocess.exe 1 scotch/proc1/mesh_for_progc.data scotch/proc1/ 
./Preprocess.exe 2 scotch/proc2/mesh_for_progc.data scotch/proc2/ 
./Preprocess.exe 3 scotch/proc3/mesh_for_progc.data scotch/proc3/
./Preprocess.exe 4 scotch/proc4/mesh_for_progc.data scotch/proc4/
./Preprocess.exe 5 scotch/proc5/mesh_for_progc.data scotch/proc5/
./Preprocess.exe 10 scotch/proc10/mesh_for_progc.data scotch/proc10/
./Preprocess.exe 15 scotch/proc15/mesh_for_progc.data scotch/proc15/
./Preprocess.exe 20 scotch/proc20/mesh_for_progc.data scotch/proc20/
./Preprocess.exe 24 scotch/proc24/mesh_for_progc.data scotch/proc24/


mpicc -lm -o fem.exe FemPar.c

# #Run FemPar.c for metis
mpirun -np 1 ./fem.exe metis/proc1/ ./time/metis/
mpirun -np 2 ./fem.exe metis/proc2/ ./time/metis/
mpirun -np 3 ./fem.exe metis/proc3/ ./time/metis/
mpirun -np 4 ./fem.exe metis/proc4/ ./time/metis/
# mpirun -np 5 ./fem.exe metis/proc5/
# mpirun -np 10 ./fem.exe metis/proc10/
# mpirun -np 15 ./fem.exe metis/proc15/
# mpirun -np 20 ./fem.exe metis/proc20/
# mpirun -np 24 ./fem.exe metis/proc24/

# #Run FemPar.c for scotch
mpirun -np 1 ./fem.exe scotch/proc1/ ./time/scotch/
mpirun -np 2 ./fem.exe scotch/proc2/ ./time/scotch/
mpirun -np 3 ./fem.exe scotch/proc3/ ./time/scotch/
mpirun -np 4 ./fem.exe scotch/proc4/ ./time/scotch/
# mpirun -np 5 ./fem.exe shotch/proc5/
# mpirun -np 10 ./fem.exe shotch/proc10/
# mpirun -np 15 ./fem.exe shotch/proc15/
# mpirun -np 20 ./fem.exe shotch/proc20/
# mpirun -np 24 ./fem.exe shotch/proc24/















