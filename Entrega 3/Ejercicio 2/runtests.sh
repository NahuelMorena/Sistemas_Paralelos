#!/bin/sh

set -e

N=512
B=8
T=4
P=4

testsln() {
    local out="${TMPDIR}/output_$1.mat"

    if [ ! -f "$2" ]; then
        echo "Compilando..."
        ./build.sh
    fi

    # Ejecuta la solucion
    local t=$(mpirun -np ${P} $2 ${N} ${B} ${T} ${IMAT} ${out})

    # Compara los resultados
    echo -n "$1: "
    ./matrixtool.py compare ${N} ${OMAT} ${out} && echo -n "OK" || echo -n "FAILED"
    echo " (time = ${t}s)"
    #./matrixtool.py print ${N} ${OMAT} ${out}
}

# Crea un directorio temporal
TMPDIR="$(mktemp -d)"

# Crea archivo con las matrices de entrada
IMAT="${TMPDIR}/input.mat"
./matrixtool.py generate ${N} "${IMAT}" f64 f64 f64 i32

# Crea archivo con el resultado de referencia
OMAT="${TMPDIR}/output.mat"
./matrixtool.py calculate ${N} "${IMAT}" "${OMAT}"

testsln mpi mat_mpi
testsln mpi_hybrid mat_hybrid

# Elimina el directorio temporal
rm -r "${TMPDIR}"