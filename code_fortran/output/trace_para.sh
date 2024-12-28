nit=$1

for ((n=0; n<=nit; n++))
do
    cat dat/sol.me*.tps$n.dat > dat/sol.tps$n.dat
    cat dat/exact/sol.me*.tps$n.dat > dat/exact/sol.tps$n.dat

    rm -f dat/sol.me*.tps$n.dat
    rm -f dat/exact/sol.me*.tps$n.dat
done

#gnuplot trace.gnu
