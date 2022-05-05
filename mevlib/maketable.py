

import sys

from mevlib.mechanisms import Mechanism
from mevlib.diagonalization import DiagonalizationSet
from mevlib.parsing.auto import parse_dynamic

from mevlib.outfmt.fortran import w_rate_f03, w_rate_f90
from mevlib.outfmt.matlab import w_rate_mat, w_ints_mat, w_diag_mat, w_full_mat
from mevlib.outfmt.python import w_rate_pkl
from mevlib.outfmt.binary import w_rate_bin, w_ints_bin, w_diag_bin, w_full_bin


# extension, rates (single matrix), ints (single stage),
#       diag (Dbar/L^2, R, lambdas, B_s R), full (diag and Rinv)
outfmts = {
    'f03': ('f03', w_rate_f03, None,       None,       None      ), #noqa E202
    'f90': ('f90', w_rate_f90, None,       None,       None      ), #noqa E202
    'mat': ('mat', w_rate_mat, w_ints_mat, w_diag_mat, w_full_mat), #noqa E202
    'bin': ('dat', w_rate_bin, w_ints_bin, w_diag_bin, w_full_bin), #noqa E202
    'pkl': ('pkl', w_rate_pkl, None,       None,       None      ), #noqa E202
    # TODO add python module
}
OUTNAME = "mevtable_{}.{}"


# there are three things we might want to write here: effectiveness factors,
# diagonalizations, or rate matrices


class OutputFormatError(Exception):
    pass


def parse_file(infile, verb=True):

    (
        precision, shape, temperatures, tempspacing, species, reactions
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

    return precision, shape, temperatures, tempspacing, mech


def compute_remesh(lambdas):

    print("TODO implement nontrivial remeshing")

    # warn about huge numbers
    for x in lambdas:
        if abs(x) > 1e80:
            print("Warning: lambda = {}.")

    # take out redundant zeros
    lambdas = sorted(list(filter(lambda x: x, lambdas)) + [0.0])

    return lambdas


def compute_integrals(lambdas, shpe, precision):
    return [
        1.0 if lam == 0.0 else shpe.intgtd(lam, precision) for lam in lambdas
    ]

def sort_integrals(lambdas, ints):
    return zip(*sorted([
        (lamb, intl)
        for lambs, intls in zip(lambdas, ints)
        for lamb, intl in zip(lambs, intls)
        ], key=lambda tpl: tpl[0]))


def lookup_outputformat(fmt, verb=False):

    if fmt in outfmts.keys():
        if verb:
            print("Using output format '{}'.".format(fmt))
        return outfmts[fmt]
    else:
        raise OutputFormatError("Unknown output format '{}'.".format(fmt))


def make_table(
    src, fmt,
    writeall, writerate, writeints, writediag, writefull,
    verb
):

    ext, ratef, intsf, diagf, fullf = lookup_outputformat(fmt)

    if not any([
        writeall, writerate, writeints, writediag, writefull
    ]):
        print("No output data specified; exiting...")
        sys.exit(0)

    if writeall:
        (
            writerate, writeints, writediag, writefull
        ) = True, True, True, True

    precision, shape, Ts, Ts_type, mech = parse_file(src, verb=verb)

    diagset = DiagonalizationSet(mech, Ts, spacing_type=Ts_type)

    if writeints:
        if intsf is None:
            print("Method 'ints' not available for this format.")
        else:
            # TODO allow a lambda range as inputs for lambda only output
            lambdas = compute_remesh(sorted([
                lm for lams in diagset.get_evals() for lm in lams
            ]))
            if verb:
                print("Final mesh has {} lambda values.".format(len(lambdas)))
                print("Computing integrals...")
            ints = compute_integrals(lambdas, shape, precision)
            intsf(OUTNAME.format('ints', ext), lambdas, ints, verb=verb)

    for name, longname, reqd, outf in [
        ('rate', 'rates',                writerate, ratef),
        ('diag', 'diagonalization',      writediag, diagf),
        ('full', 'full diagonalization', writefull, fullf)
    ]:
        if reqd:
            if outf is None:
                print(
                    "Method '{}' not available for this format.".format(name)
                )
            else:
                if verb:
                    print("Computing {}...".format(longname))
                outf(
                    "mevtable_{}.{}".format(name, ext),
                    diagset, shape, precision, verb=verb
                )

    sys.exit(0)


