#!/bin/sh

set -e

run() {
    gcc -O3 -o $2 $1 -lm > /dev/null
    ./$2
}

mkdir -p gen

echo "TIMES,quadatric2,quadatric3" > resultados.csv
for i in quadatric2 quadatric3; do
    for j in $(seq 100 100 1000); do
        echo "Compilando y ejecutando '$i.c' (TIMES = $j)..."
        f=$i-$j
        sed "s/#define TIMES .*/#define TIMES $j/g" $i.c > gen/$f.c
        echo "$j,$(run gen/$f.c gen/$f)" >> resultados.csv
    done
done
