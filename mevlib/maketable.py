

import pkgutil
from tqdm import tqdm

from mevlib.diagonalization import computetransform
from mevlib.parsing.auto import parse_dynamic


# main


#TODO add other output langauge targets and further separate that from this (a
# MVC paradigm seems appropriate here, with the controller located in main.py).


def make_table(infile, outfile, verb=True, keepnonreacting=True):

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
        if keepnonreacting:
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
            dropnonreacting=(not keepnonreacting)
        )
        for T in (tqdm(temperatures) if verb else temperatures)
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

    outfile.write(program)


