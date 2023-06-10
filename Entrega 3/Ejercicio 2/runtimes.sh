#!/bin/bash

N=(512 1024 2048 4096)
T=(8 16 32)
P=(1 2 4) #Nodos

join() {
    local v=$(printf "%s$1" "${@:2}")
    echo ${v::-1}
}

run() {
    local out="$1_times.csv"

    if [ ! -f "$1" ]; then
        echo "Compilando..."
        ./build.sh
    fi

    join ';' T ${N[@]} > ${out}
    for t in ${T[@]}; do
        [ ! -z "$2" ] && [ $t != 1 ] && break
        local a=()
        for p in ${P[@]}; do
            ([$p == 1] || [$t == 8]) && ["$1" == "mat_hybrid"] && break
            #if (p == 1 && "$1" == "mat_hybrid"); then
            #    continue;
            #fi
            for n in ${N[@]}; do
                echo -n "$1 T=$t P=$p N=$n: "

                #Cambiar para poder ejecutar en el cluster
                local v=$(mpirun -np $p $1 $n 64 $t)
                ##

                echo "${v}s"
                a+=($v)
            done
            join ';' $t ${a[@]} >> ${out}
        done
    done
}


run mat_mpi
run mat_hybrid

