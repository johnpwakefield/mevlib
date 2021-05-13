

from tqdm import tqdm

from mevlib.diagonalization import computetransform
from mevlib.parsing.auto import parse_dynamic

from mevlib.outformats.fortran import f03, f90
#TODO add mat and mex
#from mevlib.outformats.matlab import mat, mex


outfmts = {
    'f03': f03,
    'f90': f90,
#   'mat': mat,
#   'mex': mex
}


class OutputFormatError(Exception):
    pass


def compute_matrices(infile, keep_products, verb=True):

    # parse file
    precision, shape, temperatures, mech = parse_dynamic(
        infile, verb
    )
    if verb:
        print("Successfully read species {} with reactions {}.".format(
            ", ".join(map(str, mech.spcs)), ", ".join(map(str, mech.rxns))
        ))

    # print fun facts
    if verb:
        undef = len(mech.findundefined(True))
        if undef:
            print((
                "There are{} any undefined species that appear in reactions."
            ).format(undef))
        print("This reaction is{} reversible.".format(
            "n't" if mech.reversible() else ""
        ))
        mech.findnonreacting(True)
        mech.findproducts(True)

    # make data table
    if verb:
        print("Computing...")
        if keep_products:
            print("Nonreacting species will be kept.")
        else:
            print("Nonreacting species will be dropped.")
            #TODO
#       print((
#           "The resulting function will expect a {}-vector and return a "
#           "{}-vector."
#       ).format())
    matrices = [
        computetransform(
            shape, mech, T, precision,
            dropnonreacting=(not keep_products)
        )
        for T in (tqdm(temperatures) if verb else temperatures)
    ]

    return matrices, temperatures


def make_table(src, dst, ext, keep_products, verb=True):
    if ext is None and len(dst.name.split('.')) > 1:
        ext = dst.name.split('.')[-1]
    if ext is None:
        print(
            "Could not infer output format; this is typically a result of"
            "specifying an output file without a file extension (e.g. stdout)"
            "and not specifying a format.  Try the previous command with the"
            "additional option `--fmt <format>' where <format> is one of: "
            + ", ".join(outfmts.keys())
        )
        raise OutputFormatError("Output format not specified.")
    if ext in outfmts.keys():
        writer = outfmts[ext]
    else:
        raise OutputFormatError("Unknown output format '{}'.".format(ext))
    writer(dst, *compute_matrices(src, verb, keep_products), verb)


