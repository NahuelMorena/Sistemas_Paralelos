#!/bin/sh

N=256
B=8
T=4

testsln() {
    local exe="${TMPDIR}/$1"
    local out="${TMPDIR}/output_$1.mat"
    # Compila la solucion
    gcc -fopenmp -O3 -o "${exe}" "$1/"*.c > /dev/null

    # Ejecuta la solucion
    "${exe}" ${N} ${B} $2 ${IMAT} ${out}

    # Compara los resultados
    echo "$1:"
    ./matrixtool.py print ${N} ${OMAT} ${out}
    ./matrixtool.py compare ${N} ${OMAT} ${out}
    echo
}

# Crea un directorio temporal
TMPDIR="$(mktemp -d)"

# Crea archivo con las matrices de entrada
IMAT="${TMPDIR}/input.mat"
./matrixtool.py generate ${N} "${IMAT}" f64 f64 f64 i32

# Crea archivo con el resultado de referencia
OMAT="${TMPDIR}/output.mat"
./matrixtool.py calculate ${N} "${IMAT}" "${OMAT}"

testsln Secuencial
#testsln Pthreads ${T}
#testsln OpenMP ${T}

# Elimina el directorio temporal
rm -r "${TMPDIR}"
