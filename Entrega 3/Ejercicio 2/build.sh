#!/bin/sh

build() {
    local n="$1"
    shift
    mpicc "$@" -O3 -o "${n}" "src/${n}.c"
}

build mat_mpi 
build mat_hybrid -fopenmp
