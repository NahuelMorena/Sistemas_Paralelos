#!/usr/bin/env python3

import numpy
import numpy.random
import sys


_PROGNAME = '(null)'


def usage():
    print(f'usage: {_PROGNAME} print size files...')
    print(f'usage: {_PROGNAME} compare size files...')
    print(f'usage: {_PROGNAME} generate size output_file types...')
    print(f'usage: {_PROGNAME} calculate size input_file output_file')
    sys.exit(1)


def getnptype(t):
    nptypes = {
        'i32': numpy.int32,
        'i64': numpy.int64,
        'f32': numpy.float32,
        'f64': numpy.float64
    }

    return nptypes[t]


def loadmat(fp, size, t):
    dtype = getnptype(t)
    
    m = numpy.fromfile(fp, dtype=dtype, count=(size * size))
    m = m.reshape(size, size)

    return m


def calc(ma, mb, mc, md):
    min_a = numpy.min(ma)
    max_a = numpy.max(ma)
    sum_a = numpy.sum(ma)
    avg_a = sum_a / ma.size
    min_b = numpy.min(mb)
    max_b = numpy.max(mb)
    sum_b = numpy.sum(mb)
    avg_b = sum_b / mb.size

    e = (max_a * max_b - min_a * min_b) / (avg_a * avg_b)

    m_ab = numpy.matmul(ma, mb)
    m_cd = numpy.matmul(mc, (md * md))

    return m_ab * e + m_cd


def prt(args):
    if len(args) < 2:
        usage()
    
    n = int(args[0])
    files = args[1:]

    matrices = []
    for filename in files:
        with open(filename, 'rb') as fp:
            m = loadmat(fp, n, 'f64')
        matrices.append((filename, m))
    
    for filename, mat in matrices:
        print(f'{filename}:')
        print(mat)
        print()


def compare(args):
    if len(args) < 3:
        usage()
    
    n = int(args[0])
    files = args[1:]

    matrices = []
    for filename in files:
        with open(filename, 'rb') as fp:
            m = loadmat(fp, n, 'f64')
        matrices.append(m)
    
    equal = True
    base, *matrices = matrices
    for m in matrices:
        if numpy.max(numpy.abs(m - base)) > 1e-10:
            equal = False
            break

    return equal


def generate(args):
    if len(args) < 3:
        usage()
    
    n = int(args[0])
    ofile = args[1]
    types = args[2:]

    matrices = []
    for t in types:
        if t == 'i32':
            m = numpy.random.randint(1, 40, size=(n, n), dtype=numpy.int32)
        elif t == 'i64':
            maxval = (1 << 63) - 1
            m = numpy.random.randint(1, 40, size=(n, n), dtype=numpy.int64)
        elif t == 'f32':
            maxval = (1 << 20) - 1
            m = numpy.random.randint(0, maxval, size=(n, n), dtype=numpy.int32)
            m = m.astype(numpy.float32, copy=False) / maxval
        elif t == 'f64':
            maxval = (1 << 20) - 1
            m = numpy.random.randint(0, maxval, size=(n, n), dtype=numpy.int32)
            m = m.astype(numpy.float64, copy=False) / maxval

        matrices.append(m)
    
    with open(ofile, 'wb') as fp:
        for m in matrices:
            m.tofile(fp)


def calculate(args):
    if len(args) < 3:
        usage()
    
    n = int(args[0])
    ifile = args[1]
    ofile = args[2]

    matrices = []
    with open(ifile, 'rb') as fp:
        for t in ('f64', 'f64', 'f64', 'i32'):
            m = loadmat(fp, n, t)
            matrices.append(m)
    
    mr = calc(*matrices)

    with open(ofile, 'wb') as fp:
        mr.tofile(fp)


_COMMANDS = {
    'print': prt,
    'compare': compare,
    'generate': generate,
    'calculate': calculate,
}


def main(args):
    global _PROGNAME
    _PROGNAME = args[0]

    command = _COMMANDS.get(args[1]) if len(args) >= 2 else None
    if command is None:
        usage()
    
    r = command(args[2:])
    sys.exit(0 if r is not None and r else 1)


if __name__ == '__main__':
    main(sys.argv)