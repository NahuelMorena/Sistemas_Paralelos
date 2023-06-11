#!/bin/sh

build() {
    local n="$1"
    shift
    mpicc "$@" -O3 -o "${n}" "src/${n}.c"
}

for i in blocking non-blocking non-blocking-nowait blocking-ring non-blocking-ring; do
    build $i
done
