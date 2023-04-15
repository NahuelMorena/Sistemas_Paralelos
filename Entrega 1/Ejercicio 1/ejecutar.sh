#!/bin/sh

set -e

compile_run() {
    gcc -O$3 -o $2 $1 -lm > /dev/null
    ./$2
}

prueba() {
    f=$1
    t=$2

    g="gen/$f"_"$t"
    sed "s/#define TIMES .*/#define TIMES $t/g" $f.c > $g.c
    for o in $3; do
        echo "Compilando y ejecutando '$f.c' (TIMES = $t, O = $o)..."
        echo "$f;$t;-O$o;$(compile_run $g.c $g"_"O$o $o)" | tr '.' ',' >> resultado.csv
    done
}

mkdir -p gen

echo "ARCHIVO;TIMES;OPT;timed;timef" > resultado.csv

echo "Pruebas con optimizacion..."
for f in quadatric2 quadatric3; do
    for t in $(seq 200 100 1000); do
        prueba $f $t "s 1 2 3"
    done
done
