

import numpy as np


# file type (6 bytes), ordering ('F' or 'C', 1 byte), version (1 byte)
#mgn_diag = "BR=RΛF1"
#mgn_ints = "∫udVF1"
#mgn_rate = "<ω̇>F1"
mgn_diag = bytes([0x42, 0x52, 0x3d, 0x52, 0xce, 0x9b, 0x46, 0x01])
mgn_ints = bytes([0xe2, 0x88, 0xab, 0x75, 0x64, 0x56, 0x46, 0x01])
mgn_rate = bytes([0x3c, 0xcf, 0x89, 0xcc, 0x87, 0x3e, 0x46, 0x01])


def symlist_withlen(syms):
    return [
            c for b in 
            zip(
                [len(sym).to_bytes(1, 'little') for sym in syms],
                [bytearray(sym, "ascii") for sym in syms]
                )
            for c in b
            ]

def symlist_fixedlen(syms, l=8):
    return [bytearray(sym.ljust(l), "ascii") for sym in syms]


# note that np.tofile saves in C order regardless of ordering in memory

def write_rate_bin(outfile, spcsyms, temperatures, ratemats, verb=False):

    m, (n, _) = len(temperatures), ratemats[0].shape

    datamat = np.empty((n, n, m), order='F')
    for j in range(m):
        datamat[:,:,j] = ratemats[j]

    with open(outfile, 'wb') as f:
        f.write(mgn_rate)
        f.write(n.to_bytes(4, 'little'))
        f.write(m.to_bytes(4, 'little'))
        for bits in symlist_fixedlen(spcsyms):
            f.write(bits)
        np.array(temperatures).tofile(f)
        datamat.T.tofile(f)


def write_diag_bin(outfile, spcsyms, temperatures, lambdas, Rs, Rinvs, verb=False):

    m, (n, _) = len(temperatures), Rs[0].shape

    datamat = np.empty((2*n+1, n, m), order='F')
    for j in range(m):
        datamat[0,:,j] = lambdas[j]
        datamat[1:(n+1),:,j] = Rs[j]
        datamat[(n+1):(2*n+1),:,j] = Rinvs[j]

    with open(outfile, 'wb') as f:
        f.write(mgn_diag)
        f.write(n.to_bytes(4, 'little'))
        f.write(m.to_bytes(4, 'little'))
        for bits in symlist_fixedlen(spcsyms):
            f.write(bits)
        np.array(temperatures).tofile(f)
        datamat.T.tofile(f)


def write_ints_bin(outfile, lambdas, ints, verb=False):

    with open(outfile, 'wb') as f:
        f.write(mgn_ints)
        f.write(len(lambdas).to_bytes(4, 'little'))
        np.array(lambdas).tofile(f)
        np.array(ints).tofile(f)


