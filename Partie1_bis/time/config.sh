#!/bin/bash

#Nombre de proc total (24 sur PlaFrim)
nbp_tot=$1

for ((n=1; n<=nbp_tot; n++))
do
    cat ./metis/proc.*.nbproc.$n.plt > ./metis/nbproc.$n.plt
    cat ./scotch/proc.*.nbproc.$n.plt > ./scotch/nbproc.$n.plt

    rm -f ./metis/proc.*.nbproc.$n.plt
    rm -f ./scotch/proc.*.nbproc.$n.plt
done

python3 trace.py $nbp_tot