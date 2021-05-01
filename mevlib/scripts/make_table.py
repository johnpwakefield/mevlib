#!/usr/bin/env python3


import os
import pkgutil
import argparse
from tqdm import tqdm


from mevlib.parsing.auto import parse_dynamic


# input arguments and usage messages

parser = argparse.ArgumentParser(description=(
    """Generates a lookup table given a description of chemical kinetics to """
    """be used for simulation of diffusion-limited reactions in pellets.    """
))
parser.add_argument(
    'fname', metavar='fname', type=str, help=(
        """Name of input file name. Currently only our own file format      """
        """(.sbl files) are supported, but others may be added in the       """
        """future.                                                          """
    )
)
parser.add_argument(
    '-v', '--verbose', action='store_true', help=(
        """Print additional information related to the parsing of generated """
        """files and generation of the lookup table.                        """
    )
)
parser.add_argument(
    '-k', '--keepnonreacting', action='store_true', help=(
        """By default nonreacting species are dropped from the resulting    """
        """transforms to save space.  Use this option to keep them instead. """
    )
)
args = parser.parse_args()


# main

def main():

    # parse file
    if not os.path.isfile(args.fname):
        print("'{}' does not appear to be a file.".format(args.fname))
    precision, shape, temperatures, mech = parse_dynamic(
        args.fname, args.verbose
    )
    if args.verbose:
        print("Successfully read species {} with reactions {}.".format(
            ", ".join(map(str, mech.spcs)), ", ".join(map(str, mech.rxns))
        ))

    # print fun facts
    if args.verbose:
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
    if args.verbose:
        print("Computing...")
        if args.keepnonreacting:
            print("Nonreacting species will be kept.")
        else:
            print("Nonreacting species will be dropped.")
            #TODO
#       print((
#           "The resulting function will expect a {}-vector and return a "
#           "{}-vector."
#       ).format())
    matrices = [
        shape.computetransform(
            mech, T, precision['ntrunc'], precision['ktrunc'],
            dropnonreacting=(not args.keepnonreacting)
        )
        for T in (tqdm(temperatures) if args.verbose else temperatures)
    ]

    setd  = r"this%d  = {}".format(matrices[0].shape[1])
    setad = r"this%ad = {}".format(matrices[0].shape[0])
    allav = r"allocate(this%axisvals({}))".format(len(temperatures))
    alldv = r"allocate(this%datavals(this%ad, this%d, {}))".format(
        len(temperatures)
    )
    setav = r"this%axisvals = (/ {} /)".format(
        ", ".join(map(str, temperatures))
    )
    setdvs = [
        r"this%datavals(:, :, {}) = reshape((/ {} /), (/ {} /))".format(
            i + 1,
            ", ".join(map(str, matrices[i].flatten('F'))),
            ", ".join(map(str, matrices[i].shape))
        )
        for i in range(len(temperatures))
    ]

    def combinelines(limit, prespaces, lines):
        def splitline(line):
            inwords = line.split()
            outwords = []
            while len(inwords):
                nl = []
                for w in inwords:
                    if limit - prespaces - 2 - len(" ".join(nl)) - len(w) > 0:
                        nl.append(w)
                    else:
                        break
                if not len(nl):
                    print(
                        "Cannot break into lines; either lines are too short "
                        "or words are too big."
                    )
                    raise Exception("Cannot break lines.")
                outwords.append(nl)
                inwords = inwords[len(nl):]
            return [
                (" ".join(nl), i + 1 == len(outwords))
                for i, nl in enumerate(outwords)
            ]
        return "\n".join([
            prespaces * " " + nl + (" &" if not last else "")
            for ol in lines for nl, last in splitline(ol)
        ])

    initlines = combinelines(
        79, 8, [setd, setad, allav, alldv, setav] + setdvs
    )

    try:
        tpl = pkgutil.get_data('mevlib', '/templates/mevlookup_template.f03')
        if tpl is None:
            raise FileNotFoundError("template not found")
        program = tpl.decode().replace('!inline_init_here!', initlines)
    except FileNotFoundError:
        print("Could not find template file.")
        exit(1)

    if args.verbose:
        print("Writing output file 'mevlookup.f03'.")
    with open("mevlookup.f03", 'w') as f:
        f.write(program)


main()

