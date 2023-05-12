#!/bin/bash

N=(512 1024 2048 4096)
T=(1 2 4 8)

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
        for n in ${N[@]}; do
            echo -n "$1 T=$t N=$n: "
            local v=$(./$1 $n 8 $t)
            echo "${v}s"
            a+=($v)
        done
        join ';' $t ${a[@]} >> ${out}
    done
}

run mat 1
run mat_openmp
run mat_pthreads