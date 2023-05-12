#!/bin/sh

build() {
    local n="$1"
    shift
    gcc "$@" -O3 -o "${n}" "src/${n}.c"
}

build mat
build mat_openmp -fopenmp
build mat_pthreads -pthread