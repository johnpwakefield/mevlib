

# TODO it would be better if all the math was migrated to diagonalization
import numpy as np
import scipy.linalg as la

from tqdm import tqdm

from mevlib.mechanisms import Mechanism
from mevlib.diagonalization import computediagonalization, computeintegral
from mevlib.parsing.auto import parse_dynamic

from mevlib.outformats.fortran import write_mat_f03, write_mat_f90
# TODO add mat and mex
from mevlib.outformats.matlab import write_mat_mat
from mevlib.outformats.python import write_mat_pkl
from mevlib.outformats.binary import write_rate_bin
from mevlib.outformats.binary import write_diag_bin
from mevlib.outformats.binary import write_ints_bin


# extension, all, mevs, rates, ints, diag
outfmts = {
    'f03': ('f03', None, write_mat_f03, write_mat_f03, None, None),
    'f90': ('f90', None, write_mat_f90, write_mat_f90, None, None),
    'mat': ('mat', None, write_mat_mat, None, None, None),
    'bin': ('dat', None, None, write_rate_bin, write_ints_bin, write_diag_bin),
    'pkl': ('pickle', None, write_mat_pkl, None, None, None)
    # TODO add python module
}


# there are three things we might want to write here: effectiveness factors,
# diagonalizations, or rate matrices


class OutputFormatError(Exception):
    pass


def parse_file(infile, verb=True):

    (
        precision, shape, temperatures, species, reactions
    ) = parse_dynamic(infile, verb)
    mech = Mechanism(species, reactions)

    if verb:
        print("Successfully read species {} with reactions {}.".format(
            ", ".join(map(str, mech.spcs)), ", ".join(map(str, mech.rxns))
        ))
        undef = len(mech.findundefined(True))
        if undef:
            print((
                "There are{} any undefined species that appear in reactions."
            ).format(undef))
        print("This reaction is{} reversible.".format(
            "n't" if mech.isreversible() else ""
        ))
        mech.findnonreacting(True)
        mech.findproducts(True)

    return precision, shape, temperatures, mech


def compute_diag(temperatures, mech, verb=True):

    if verb:
        print("Computing diagonalizations...")

    lambdas, Rs = zip(*[
        computediagonalization(mech, T)
        for T in (tqdm(temperatures) if verb else temperatures)
        ])

    return lambdas, Rs


def compute_remesh(lambdas):

    print("TODO remeshing not yet implemented; returning original mesh")

    return lambdas


def compute_integrals(lambdas, shpe, precision, verb=True):

    if verb:
        print("Computing integrals...")

    ints = [computeintegral(shpe, lam, precision) for lam in lambdas]

    return ints


def sort_integrals(lambdas, ints, verb=True):

    if verb:
        print("Sorting integrals...")

    stdints = sorted([
        (lamb, intl)
        for lambs, intls in zip(lambdas, ints)
        for lamb, intl in zip(lambs, intls)
        ], key=lambda tpl: tpl[0])

    return zip(*stdints)


def compute_effvmatrices(ints, Rs, verb=True):

    if verb:
        print("Generating MEVs...")

    As = [
            np.dot(R, np.dot(np.diag(intg), np.linalg.inv(R)))
            for intg, R in zip(ints, Rs)
            ]

    return As


def compute_ratematrices(lambdas, ints, Rs, verb=True):

    if verb:
        print("Generating rate matrices...")

    Bs = [
            np.dot(R, np.dot(lams * np.diag(intg), np.linalg.inv(R)))
            for lams, intg, R in zip(lambdas, ints, Rs)
            ]

    return Bs


def compute_augmenteddiag(spcs, Ts, Rs, Rinvs, ndlen=1.0, verb=True):

    if verb:
        print("Computing augmented diagonalization matrices...")

    Es = [np.dot(np.diag([
        spc.effective_diffusion(T) / ndlen**2 for spc in spcs
        ]), R) for T, R in zip(Ts, Rs)]

    return Es, Rinvs


def make_table(
    src, fmt,
    writeall, writeint, writediag, writemev, writeaugd, writerate,
    verb
):

    if fmt in outfmts.keys():
        ext, allf, mevsf, ratesf, intsf, diagf = outfmts[fmt]
    else:
        raise OutputFormatError("Unknown output format '{}'.".format(fmt))

    if writeall:
        if allf is None:
            writeall = False
            (
                writeint, writediag, writemev, writeaugd, writerate
            ) = True, True, True, True

    precision, shape, Ts, mech = parse_file(src, verb=verb)

    # TODO allow a lambda range as inputs for lambda only output
    lambdas, Rs = compute_diag(Ts, mech, verb=verb)
    if writeint or writemev or writerate or writeaugd:
        ints = compute_integrals(lambdas, shape, precision, verb=verb)
    if writemev:
        As = compute_effvmatrices(ints, Rs, verb=verb)
    if writerate:
        Bs = compute_ratematrices(lambdas, ints, Rs, verb=verb)
    if writediag or writemev or writerate or writeaugd:
        Rinvs = [la.inv(R) for R in Rs]
    if writeaugd:
        augd = compute_augmenteddiag(mech.spcs, Ts, Rs, Rinvs, verb=verb)
    syms = [spc.symb for spc in mech.spcs]

    if writeall:
        if allf is not None:
            allf("mevtable_all."+ext, syms, Ts, Bs, verb=verb)
        else:
            print("Method 'all' not available for this format.")

    if writeint:
        if intsf is not None:
            intsf(
                "mevtable_ints."+ext,
                *sort_integrals(lambdas, ints, verb=verb),
                verb=verb
            )
        else:
            print("Method 'ints' not available for this format.")

    if writediag:
        if diagf is not None:
            diagf(
                "mevtable_diag."+ext, syms, Ts, lambdas, Rs, Rinvs,
                verb=verb
            )
        else:
            print("Method 'diag' not available for this format.")

    if writemev:
        if mevsf is not None:
            mevsf("mevtable_mevs."+ext, syms, Ts, As, verb=verb)
        else:
            print("Method 'mevs' not available for this format.")

    if writeaugd:
        if augd is not None:
            diagf("mevtable_augd."+ext, syms, Ts, lambdas, *augd, verb=verb)
        else:
            print("Method 'augd' not available for this format.")

    # TODO why don't we have the same checking here?
    if writerate:
        ratesf("mevtable_rates."+ext, syms, Ts, Bs, verb=verb)


