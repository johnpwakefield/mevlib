

from functools import reduce


import numpy as np


# file type (6 bytes), ordering ('F' or 'C', 1 byte), version (1 byte)
mgn_rate = bytes([0x3c, 0xcf, 0x89, 0xcc, 0x87, 0x3e, 0x46, 0x01])
mgn_ints = bytes([0xe2, 0x88, 0xab, 0x75, 0x64, 0x56, 0x46, 0x01])
mgn_diag = bytes([0x42, 0x52, 0x3d, 0x52, 0xce, 0x9b, 0x46, 0x01])
mgn_full = bytes([0x52, 0x5e, 0x7b, 0x2d, 0x31, 0x7d, 0x46, 0x01])


def symlist_fixedlen(syms, s=8):
    return [bytearray(sym.ljust(s), "ascii") for sym in syms]


# note that np.tofile saves in C order regardless of ordering in memory

def w_rate_bin(outfile, diagset, shape, precision, verb=False):

    M, Ng, Ns = len(diagset.Ts), diagset.mech.Ng, diagset.mech.Ns

    spcsyms = [spc.symb for spc in diagset.mech.spcs]
    datamat = np.empty((Ng + Ns, Ng, M), order='F')
    for j, A in enumerate(diagset.get_transforms(shape, precision)):
        datamat[:, :, j] = A

    with open(outfile, 'wb') as f:
        f.write(mgn_rate)
        f.write(Ng.to_bytes(1, 'little'))
        f.write(Ns.to_bytes(1, 'little'))
        f.write(M.to_bytes(2, 'little'))
        for bits in symlist_fixedlen(spcsyms, s=8):
            f.write(bits)
        np.array(diagset.Ts).tofile(f)
        datamat.T.tofile(f)


def w_ints_bin(outfile, lambdas, ints, verb=False):

    with open(outfile, 'wb') as f:
        f.write(mgn_ints)
        f.write(len(lambdas).to_bytes(4, 'little'))
        np.array(lambdas).tofile(f)
        np.array(ints).tofile(f)


def w_diag_bin(outfile, diagset, shape, precision, verb=False):

    M, Ng, Ns = len(diagset.Ts), diagset.mech.Ng, diagset.mech.Ns

    spcsyms = [spc.symb for spc in diagset.mech.spcs]

    with open(outfile, 'wb') as f:
        f.write(mgn_diag)
        f.write(Ng.to_bytes(1, 'little'))
        f.write(Ns.to_bytes(1, 'little'))
        f.write(M.to_bytes(2, 'little'))
        for bits in symlist_fixedlen(spcsyms, s=8):
            f.write(bits)
        np.array(diagset.Ts).tofile(f)
        for D, lams, R, Bs in zip(
            diagset.get_Dis(), diagset.get_evals(),
            diagset.get_evects(), diagset.Bss
        ):
            assert(np.all(
                np.dot(
                    D.reshape((1,-1)),
                    np.vstack((
                        np.dot(
                            np.dot(R, np.diag(lams)),
                            np.linalg.inv(R)
                        ),
                        Bs
                    ))
                )
            ) < 1e-12)
            D.tofile(f)
            lams.tofile(f)
            R.T.tofile(f)
            np.dot(Bs, R).T.tofile(f)


def w_full_bin(outfile, diagset, shape, precision, verb=False):

    M, Ng, Ns = len(diagset.Ts), diagset.mech.Ng, diagset.mech.Ns

    spcsyms = [spc.symb for spc in diagset.mech.spcs]

    with open(outfile, 'wb') as f:
        f.write(mgn_full)
        f.write(Ng.to_bytes(1, 'little'))
        f.write(Ns.to_bytes(1, 'little'))
        f.write(M.to_bytes(2, 'little'))
        for bits in symlist_fixedlen(spcsyms, s=8):
            f.write(bits)
        np.array(diagset.Ts).tofile(f)
        for D, lams, R, Bs, Rinv in zip(
            diagset.get_Dis(), diagset.get_evals(),
            diagset.get_evects(), diagset.Bss,
            diagset.get_evectinvs()
        ):
            D.tofile(f)
            lams.tofile(f)
            R.T.tofile(f)
            np.dot(Bs, R).T.tofile(f)
            Rinv.T.tofile(f)


