#!/usr/bin/env bash
#pocet cisel bud zadam nebo 10 :)
if [ $# -lt 1 ];then
    echo "Invalid input"
    exit
else
    nodes=$1;
fi;

STRLENGTH=$(echo -n $nodes | wc -m)


procs="$((2*$STRLENGTH -2))"
echo $nodes >> tree


##preklad cpp zdrojaku
mpic++ -o pro pro.cpp

##spusteni
mpirun -np $procs pro
##uklid
rm -f pro tree
