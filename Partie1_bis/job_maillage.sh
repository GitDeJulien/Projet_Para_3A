#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH -C miriel

# Get input parameter

# Load dependencies
#module load compiler/gcc/9.3.0
#module load mpi/openmpi/3.1.4-all
#module load compiler/gfortran/*

# cd in the project folder

gfortran data2tec.f90 -o data2tec.exe
./data2tec.exe

metis-4.0.3/mesh2dual build/dualformetis.data

metis-4.0.3/kmetis build/dualformetis.data.dgraph 2
metis-4.0.3/kmetis build/dualformetis.data.dgraph 3
metis-4.0.3/kmetis build/dualformetis.data.dgraph 4
metis-4.0.3/kmetis build/dualformetis.data.dgraph 5
metis-4.0.3/kmetis build/dualformetis.data.dgraph 10
metis-4.0.3/kmetis build/dualformetis.data.dgraph 15
metis-4.0.3/kmetis build/dualformetis.data.dgraph 20
metis-4.0.3/kmetis build/dualformetis.data.dgraph 24

mv build/dualformetis.data.dgraph.part.2 metis/proc2/dualformetis.data.dgraph.part.2
mv build/dualformetis.data.dgraph.part.3 metis/proc3/dualformetis.data.dgraph.part.3
mv build/dualformetis.data.dgraph.part.4 metis/proc4/dualformetis.data.dgraph.part.4
mv build/dualformetis.data.dgraph.part.5 metis/proc5/dualformetis.data.dgraph.part.5
mv build/dualformetis.data.dgraph.part.10 metis/proc10/dualformetis.data.dgraph.part.10
mv build/dualformetis.data.dgraph.part.15 metis/proc15/dualformetis.data.dgraph.part.15
mv build/dualformetis.data.dgraph.part.20 metis/proc20/dualformetis.data.dgraph.part.20
mv build/dualformetis.data.dgraph.part.24 metis/proc24/dualformetis.data.dgraph.part.24


g++ metisdual2scotchdual.cpp -o metisdual2scotchdual

./metisdual2scotchdual build/dualformetis.data.dgraph scotch/proc1/dualforscotch.grf
./metisdual2scotchdual build/dualformetis.data.dgraph scotch/proc2/dualforscotch.grf
./metisdual2scotchdual build/dualformetis.data.dgraph scotch/proc3/dualforscotch.grf
./metisdual2scotchdual build/dualformetis.data.dgraph scotch/proc4/dualforscotch.grf
./metisdual2scotchdual build/dualformetis.data.dgraph scotch/proc5/dualforscotch.grf
./metisdual2scotchdual build/dualformetis.data.dgraph scotch/proc10/dualforscotch.grf
./metisdual2scotchdual build/dualformetis.data.dgraph scotch/proc15/dualforscotch.grf
./metisdual2scotchdual build/dualformetis.data.dgraph scotch/proc20/dualforscotch.grf
./metisdual2scotchdual build/dualformetis.data.dgraph scotch/proc24/dualforscotch.grf

scotch_5.1.11/bin/gtst scotch/proc1/dualforscotch.grf
scotch_5.1.11/bin/gtst scotch/proc2/dualforscotch.grf
scotch_5.1.11/bin/gtst scotch/proc3/dualforscotch.grf
scotch_5.1.11/bin/gtst scotch/proc4/dualforscotch.grf
scotch_5.1.11/bin/gtst scotch/proc5/dualforscotch.grf
scotch_5.1.11/bin/gtst scotch/proc10/dualforscotch.grf
scotch_5.1.11/bin/gtst scotch/proc15/dualforscotch.grf
scotch_5.1.11/bin/gtst scotch/proc20/dualforscotch.grf
scotch_5.1.11/bin/gtst scotch/proc24/dualforscotch.grf

echo cmplt 1 | scotch_5.1.11/bin/gmap scotch/proc1/dualforscotch.grf - scotch/proc1/dualforscotch.map
echo cmplt 2 | scotch_5.1.11/bin/gmap scotch/proc2/dualforscotch.grf - scotch/proc2/dualforscotch.map
echo cmplt 3 | scotch_5.1.11/bin/gmap scotch/proc3/dualforscotch.grf - scotch/proc3/dualforscotch.map
echo cmplt 4 | scotch_5.1.11/bin/gmap scotch/proc4/dualforscotch.grf - scotch/proc4/dualforscotch.map
echo cmplt 5 | scotch_5.1.11/bin/gmap scotch/proc5/dualforscotch.grf - scotch/proc5/dualforscotch.map
echo cmplt 10 | scotch_5.1.11/bin/gmap scotch/proc10/dualforscotch.grf - scotch/proc10/dualforscotch.map
echo cmplt 15 | scotch_5.1.11/bin/gmap scotch/proc15/dualforscotch.grf - scotch/proc15/dualforscotch.map
echo cmplt 20 | scotch_5.1.11/bin/gmap scotch/proc20/dualforscotch.grf - scotch/proc20/dualforscotch.map
echo cmplt 24 | scotch_5.1.11/bin/gmap scotch/proc24/dualforscotch.grf - scotch/proc24/dualforscotch.map

gfortran postprocess_metis.f90 -o post_metis.exe
./post_metis.exe 2
./post_metis.exe 3
./post_metis.exe 4
./post_metis.exe 5
./post_metis.exe 10
./post_metis.exe 15
./post_metis.exe 20
./post_metis.exe 24

gfortran postprocess_scotch.f90 -o post_scotch.exe
./post_scotch.exe 1
./post_scotch.exe 2
./post_scotch.exe 3
./post_scotch.exe 4
./post_scotch.exe 5
./post_scotch.exe 10
./post_scotch.exe 15
./post_scotch.exe 20
./post_scotch.exe 24













