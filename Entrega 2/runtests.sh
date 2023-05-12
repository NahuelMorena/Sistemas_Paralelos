#!/bin/sh

N=256
B=8
T=4

testsln() {
    local out="${TMPDIR}/output_$1.mat"

    if [ ! -f "$2" ]; then
        echo "Compilando..."
        ./build.sh
    fi

    # Ejecuta la solucion
    ./$2 ${N} ${B} ${T} ${IMAT} ${out} > /dev/null

    # Compara los resultados
    echo -n "$1: "
    ./matrixtool.py compare ${N} ${OMAT} ${out} && echo "OK" || echo "FAILED"
}

# Crea un directorio temporal
TMPDIR="$(mktemp -d)"

# Crea archivo con las matrices de entrada
IMAT="${TMPDIR}/input.mat"
./matrixtool.py generate ${N} "${IMAT}" f64 f64 f64 i32

# Crea archivo con el resultado de referencia
OMAT="${TMPDIR}/output.mat"
./matrixtool.py calculate ${N} "${IMAT}" "${OMAT}"

testsln secuencial mat
testsln pthreads mat_pthreads
testsln openmp mat_openmp

# Elimina el directorio temporal
rm -r "${TMPDIR}"