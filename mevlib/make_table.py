#!/usr/bin/env python3


from functools import reduce, partial
from abc import ABC, abstractmethod

import os
import argparse
import pyparsing as pp
from tqdm import tqdm

from math import sqrt, exp, inf, floor

import numpy as np
# the scipy linear algebra package has some advantages, but this can but
# substituted for the numpy version if scipy is not available
from scipy.special import i0e, i1e, jn_zeros
from scipy import linalg as la


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


# file parsing

# check if is zero; ideally this will be rewritten to be type-agnostic
def iszero(x):
    return x == 0

class FileTypeException(Exception):
    pass

def parse_chemkin(fn):
    raise NotImplementedError("Chemkin parser not written.")

def parse_ini(fn):
    raise NotImplementedError("Ini parser not written.")

def parse_sensible(fn, verb=False):

    # building blocks
    def keylist(ls):
        return reduce(lambda x, y: x | y, map(pp.CaselessKeyword, ls))
    pp_end = pp.Suppress(pp.CaselessKeyword("end"))
    pp_word = pp.Word(pp.alphas)
    pp_symb = ~pp_end + pp.Word(pp.alphas, pp.alphanums)
    pp_qstr = pp.QuotedString("\"")
    pp_real = pp.pyparsing_common.real | pp.pyparsing_common.integer
    pp_prec = pp.CaselessKeyword("precision")
    pp_shpe = pp.CaselessKeyword("shape")
    pp_temp = pp.CaselessKeyword("temperature")
    pp_diff = pp.CaselessKeyword("diffusion")
    pp_spec = pp.CaselessKeyword("species")
    pp_reac = pp.CaselessKeyword("reactions")
    pp_type = pp.CaselessKeyword("type")
    pp_tble = pp.CaselessKeyword("table")
    pp_comment = pp.oneOf('! # //') + pp.restOfLine
    pp_to = pp.Suppress(pp.Literal('->') | pp.Literal('='))
    pp_value = pp_real | keylist(["default", "specified", "omitted"])
    pp_valpair = pp.Group(pp_word + pp_value)
    def simplesection(key, ls):
        pp_simp_typ = pp.Group(pp_type + keylist(ls))
        pp_simp_blk = pp.Group(
            key + pp.Group(pp.ZeroOrMore(pp_simp_typ | pp_valpair)) + pp_end
        )
        return pp_simp_blk

    # precision section
    pp_prec_blk = simplesection(pp_prec, ["trunc", "tol"])

    # shape section
    pp_shpe_blk = simplesection(pp_shpe, ["cylinder", "box", "sphere"])

    # temperature section
    pp_temp_blk = simplesection(pp_temp, ["uniform", "explicit"])

    # diffusion section
    pp_diff_typ = pp.Group(pp_type + keylist(["knudsen"]))
    pp_diff_blk = pp.Group(
        pp_diff + pp.Group(pp.ZeroOrMore(pp_diff_typ | pp_valpair)) + pp_end
    )

    # species section
    pp_spec_row = pp.Group(
        pp_symb + pp.Optional(pp_qstr) + pp.Optional(pp_real)
    )
    pp_spec_tbl = pp.Group(
        pp_tble + pp.Group(pp.ZeroOrMore(pp_spec_row)) + pp_end
    )
    pp_spec_blk = pp.Group(
        pp_spec + pp.Group(pp.ZeroOrMore(pp_spec_tbl | pp_valpair)) + pp_end
    )

    # reaction section
    pp_reac_typ = pp.Group(pp_type + pp_word)
    pp_reac_row = pp.Group(pp_symb + pp_to + pp_symb + pp.OneOrMore(pp_real))
    pp_reac_tbl = pp.Group(
        pp_tble + pp.Group(pp.ZeroOrMore(pp_reac_row)) + pp_end
    )
    pp_reac_blk = pp.Group(
        pp_reac
        + pp.Group(pp.ZeroOrMore(pp_reac_typ | pp_reac_tbl | pp_valpair))
        + pp_end
    )

    parser = pp.ZeroOrMore(
        pp_prec_blk | pp_shpe_blk | pp_temp_blk
        | pp_diff_blk | pp_spec_blk | pp_reac_blk
    ).ignore(pp_comment)

    with open(fn, "r") as f:
        d = {
            h: {k.lower(): v for k, v in u}
            for h, u in parser.parseFile(f).asList()
        }

    # make sure we have the key sections
    req_sections = [
        'diffusion', 'species', 'reactions',
        'precision', 'shape', 'temperature'
    ]
    for sec in req_sections:
        if sec not in d.keys():
            raise FileTypeException(
                "Could not parse required section '{}'.".format(sec)
            )

    # read series truncation or error tolerance
    precision = d['precision']
    if 'type' not in d['precision']:
        print("Precision type not specified.")
        raise FileTypeException("Precision type not specified.")
    if d['precision']['type'] == 'trunc':
        if not ('ntrunc' in d['precision'] and 'ktrunc' in d['precision']):
            if 'trunc' in d['precision']:
                if verb:
                    print("Single truncation given, using this for both sums.")
                precision['ntrunc'] = d['precision']['trunc']
                precision['ktrunc'] = d['precision']['trunc']
            else:
                print("Missing truncation information.")
                raise FileTypeException("Missing truncation information.")
    elif d['precision']['type'] == 'tol':
        raise NotImplementedError("Convergence tolerance not implemented yet.")

    # read shape
    if 'type' not in d['shape']:
        print("Shape not specified.")
        raise FileTypeException("Shape not specified.")
    if d['shape']['type'] == "cylinder":
        if 'radius' not in d['shape'] and 'diameter' not in d['shape']:
            print("Radius missing from shape dimensions.")
            raise FileTypeException("Missing shape dimensions.")
        if 'radius' in d['shape'] and 'diameter' in d['shape']:
            print("Radius and diameter cannot both be specified.")
            raise FileTypeException("Shape dimensions overspecified.")
        if 'height' not in d['shape']:
            print("Height missing from shape dimensions.")
            raise FileTypeException("Missing shape dimensions.")
        if 'radius' in d['shape']:
            shape = Cylinder(d['shape']['radius'], d['shape']['height'])
        else:
            shape = Cylinder(d['shape']['diameter'] / 2, d['shape']['height'])
    elif d['shape']['type'] == "box":
        raise NotImplementedError("Shape 'box' not implemented.")
    elif d['shape']['type'] == "sphere":
        raise NotImplementedError("Shape 'sphere' not implemented.")
    else:
        print("Shape '{}' not recognized.".format(d['shape']['type']))
        raise FileTypeException(
            "Shape '{}' not recognized.".format(d['shape']['type'])
        )

    # read temperature ranges
    if 'type' not in d['temperature']:
        print("Type of temperature specification not present.")
        raise FileTypeException("Missing type for temperature range.")
    if d['temperature']['type'] == "uniform":
        if 'start' not in d['temperature'] or 'stop' not in d['temperature']:
            print("Start and stop of temperature range must be given.")
            raise FileTypeException("Range parameters unspecified.")
        if 'num' not in d['temperature'] and 'step' not in d['temperature']:
            print("Either number of points or spacing must be specified.")
            raise FileTypeException("Range parameters unspecified.")
        if 'num' in d['temperature'] and 'step' in d['temperature']:
            print("Only one of 'num' or 'step' may be specified.")
            raise FileTypeException("Range parameters overspecified.")
        l = d['temperature']['stop'] - d['temperature']['start']
        if 'num' in d['temperature']:
            num = d['temperature']['num']
        else:
            num = floor(l / d['temperature']['step'])
        step = l / (num - 1)
        temperatures = [
            d['temperature']['start'] + i * step for i in range(num)
        ]
    elif d['temperature']['type'] == "explicit":
        raise NotImplementedError("Explicit temperatures not implemented.")
    else:
        print((
            "Temperature type '{}' not recognized."
        ).format(d['temperature']['type']))
        raise FileTypeException("Temperature type not recognized.")

    # read diffusion parameters
    if d['diffusion']['type'] == 'knudsen':
        knudsenparams = ['porediameter', 'voidage', 'tortuosity']
        for param in knudsenparams:
            if param not in d['diffusion'].keys():
                print("Diffusion section missing '{}'.".format(param))
                raise FileTypeException("Missing diffusion parameters.")
        speciesconstructor = partial(KnudsenSpecies, *[
            d['diffusion'][param] for param in knudsenparams
        ])
    else:
        print((
            "Did not recognize diffusion type '{}'."
        ).format(d['diffusion']['type']))
        raise FileTypeException("Unrecognized diffusion type.")

    # make species list
    if 'table' not in d['species']:
        raise FileTypeException("Missing species table.")
    name = d['species'].get('longname', 'default')
    if name == 'default':
        name == 'omitted'
    molweight = d['species'].get('molweight', 'default')
    if molweight == 'default':
        molweight == 'specified'
    def checkrow(row):
        if (
            len(row)
            != (1 + (name == 'specified') + (molweight == 'specified'))
        ):
            print("Unexpected length of table row.")
            print("Expected {} columns: {}.".format(", ".join(
                ['symbol']
                + (['name'] if name == 'specified' else [])
                + (['molweight'] if molweight == 'specified' else [])
            )))
            raise FileTypeException("Unexpected length of table row.")
    def getname(row):
        if name == 'omitted':
            return row[0]
        elif name == 'specified':
            return row[1]
        else:
            if verb:
                print((
                    "Using chemical name '{}' for all species; this was "
                    "likely unintended. Try using 'DEFAULT', 'OMITTED', or "
                    "'SPECIFIED'. "
                ).format(name))
            return name
    def getmolweight(row):
        if molweight == 'omitted':
            raise FileTypeException("Cannot ommit molecular weight.")
        elif molweight == 'specified':
            return row[2]
    def makespecies(row):
        checkrow(row)
        return speciesconstructor(row[0], getname(row), getmolweight(row))
    species = [makespecies(row) for row in d['species']['table']]

    # make reaction list
    reactiontype = d['reactions'].get('type', 'arrcoeffs').lower()
    if reactiontype == 'arrcoeffs':
        prefactor = d['reactions'].get('prefactor', 'default')
        exponent = d['reactions'].get('exponent', 'default')
        actenergy = d['reactions'].get('actenergy', 'default')
        reftemp = d['reactions'].get('reftemp', 'default')
        if prefactor == 'default':
            prefactor = 'specified'
        if prefactor == 'omitted':
            print("Prefactor may not be omitted.")
            raise FileTypeException("Prefactor omitted.")
        if exponent == 'default' or exponent == 'omitted':
            exponent = None
        if actenergy == 'default':
            actenergy = 'specified'
        if actenergy == 'omitted':
            print("Activation energy may not be omitted.")
            raise FileTypeException("Activation energy omitted.")
        if reftemp == 'default' or reftemp == 'omitted':
            reftemp = None
        expected_row_length = (
            2 + (prefactor == 'specified') + (exponent == 'specified')
            + (actenergy == 'specified') + (reftemp == 'specified')
        )
        def makereaction(row):
            if len(row) != expected_row_length:
                print("Unexpected row length in reaction table.")
                raise FileTypeException("Unexpected row length.")
            args = row[:2] + [prefactor, exponent, actenergy, reftemp]
            i = 2
            for j in range(2, 6):
                if args[j] == 'specified':
                    args[j] = row[i]
                    i += 1
                if args[j] == 'omitted':
                    args[j] = None
            return ArrheniusReaction(*args)
    else:
        print("Reaction type '{}' not recognized.".format(reactiontype))
        raise FileTypeException("Unrecognized reaction type.")
    reactions = [makereaction(row) for row in d['reactions']['table']]

    return precision, shape, temperatures, Mechanism(species, reactions)


file_types = {
#   '.inp'  : ("Chemkin-II", parse_chemkin),
#   '.ini'  : ("ini", parse_ini),
    '.sbl'  : ("sensible", parse_sensible)
}

def parse_attempt(f, ext, verb):
    name, parser = file_types[ext]
    try:
        if verb:
            print("Attempting to parse input as a {} file.".format(name))
        params = parser(f, verb=verb)
        if verb:
            print("Input successfully parsed as a {} file.".format(name))
        return params
    except (FileTypeException,) as err:
        print("File could not be parsed as a {} file.".format(name))
        print(err)
        return None

def parse_dynamic(f, verb):
    ext = f[:-4] if f[:-4] in file_types.keys() else None
    if ext is not None:
        if verb:
            print("Inferring file type from file extension.")
        params = parse_attempt(f, ext, verb)
        if params is not None:
            return params
    for k in file_types.keys():
        if k != ext:
            params = parse_attempt(f, k, verb)
            if params is not None:
                return params
    print("Unable to parse input file.")
    raise Exception("Unable to parse input file.")

class Species(ABC):
    name = None                 # string containing species name

    def __init__(self, symb, name):
        self.symb, self.name = symb, name

    def __str__(self):
        return "{} ({})".format(self.name, self.symb)

    @abstractmethod
    def effective_diffusion(self, T):
        pass

class KnudsenSpecies(Species):

    # put the species-specific parameters last so we can curry the constructor
    def __init__(self, Dpore, epspore, tau, symb, name, Mi):
        super().__init__(symb, name)
        self.Dpore, self.epspore, self.tau, self.Mi = Dpore, epspore, tau, Mi

    def __repr__(self):
        return (
            "KnudsenSpecies(Dpore = {}, epspore = {}, tau = {}, symb = {}, "
            "name = {}, Mi = {})"
        ).format(
            self.Dpore, self.epspore, self.tau, self.symb, self.name, self.Mi
        )

    def effective_diffusion(self, T):
        return 48.5 * self.Dpore * sqrt(T / self.Mi) * self.epspore / self.tau

class Reaction(ABC):
    src, dst = None, None       # names of species

    def __init__(self, src, dst):
        self.src, self.dst = src, dst

    def __str__(self):
        return "{} -> {}".format(self.src, self.dst)

    @abstractmethod
    def kij(self, T):
        pass

class ArrheniusReaction(Reaction):

    GASCONST = 8.31446261815324e-3  # kJ per mol Kelvin

    def __init__(self, src, dst, A, b, Ea, T0):
        super().__init__(src, dst)
        self.A, self.b, self.Ea, self.T0 = A, b, Ea, T0

    def __repr__(self):
        return (
            "ArrheniusReaction(src = {}, dst = {}, A = {}, b = {}, Ea = {}, "
            "T0 = {})"
        ).format(
            self.src, self.dst, self.A, self.b, self.Ea, self.T0
        )

    def kij(self, T):
        A = self.A
        if self.b is not None and self.b != 0.0:
            A *= T**self.b
        if self.T0 is None or self.T0 == inf or self.T0 == -inf:
            T0inv = 0.0
        else:
            T0inv = self.T0**(-1)
        return A * exp(-self.Ea / self.GASCONST * (1 / T - T0inv))

class Mechanism(object):

    spcs = None             # list of species objects
    rxns = None             # list of reaction objects
    numactive = None      # ignore products

    # providing an explicit constructor prevents the common error of
    # overwriting default values
    def __init__(self, spcs, rxns):
        self.spcs, self.rxns = spcs, rxns

    def reversible(self):
        edges = [(rxn.src, rxn.dst) for rxn in self.rxns]
        def findroot(es):
            ns = list(set([n for e in es for n in e]))
            for n in ns:
                ss = [s for s, d in es if n == d]
                if len(ss) == 0:
                    return n
            return None
        root = findroot(edges)
        while root is not None:
            edges = [e for e in edges if root not in e]
            root = findroot(edges)
        return len(edges) == 0

    def findnonreacting(self, verb):
        nonreacting = [
            spc.symb for spc in self.spcs
            if not any([
                spc.symb == rxn.src or spc.symb == rxn.dst
                for rxn in self.rxns
            ])
        ]
        if verb:
            if len(nonreacting):
                print("Some species are inert, they are: {}".format(
                    ", ".join(nonreacting)
                ))
            else:
                print("All species are involved in at least one reaction.")
        return nonreacting

    def findproducts(self, verb):
        products = [
            spc.symb for spc in self.spcs
            if not any([spc.symb == rxn.src for rxn in self.rxns])
        ]
        self.spcs.sort(key=lambda spc: spc.symb in products)
        self.numactive = len(self.spcs) - len(products)
        if verb:
            print((
                "Found {} product species:\n\t{}"
            ).format(len(products), ", ".join(products)))
        return products

    def getnumactive(self):
        if self.numactive is None:
            self.findproducts(False)
        return self.numactive

    def findundefined(self, verb):
        spcnames = [spc.symb for spc in self.spcs]
        undef = []
        for rxn in self.rxns:
            if rxn.src not in spcnames:
                undef.append(rxn.src)
                if verb:
                    print("Source '{}' not in species list.".format(rxn.src))
            if rxn.dst not in spcnames:
                undef.append(rxn.dst)
                if verb:
                    print(
                        "Destination '{}' not in species list.".format(rxn.dst)
                    )
        return undef

    def getmatrix(self, charlen, T):
        n = len(self.spcs)
        symbols = [spc.symb for spc in self.spcs]
        kijs, Dis = np.zeros((n, n)), np.empty((n,))
        for i, spcs in enumerate(self.spcs):
            Dis[i] = spcs.effective_diffusion(T)
        for rxn in self.rxns:
            i, j = symbols.index(rxn.src), symbols.index(rxn.dst)
            kijs[i, j] = rxn.kij(T)
        phiij2s = (
            kijs * charlen**2 / np.tile(Dis.reshape((-1,1)), len(self.spcs))
        )
        B = np.diag(np.sum(phiij2s, axis=0)) - phiij2s
        return B


# shapes

class Shape(ABC):
    pass

class Cylinder(Shape):

    def __init__(self, R, H, charlen=None):
        if charlen is None:
            charlen = 0.5 * (R * H / (R + H))
        self.R, self.H, self.charlen, self.zero_cache = R, H, charlen, None

    def gn(self, a2, n):
        return np.sqrt(a2 + (np.pi * (2 * n + 1) / self.H)**2)

    def ln(self, a2, alphak):
        return np.sqrt(a2 + (alphak / self.R)**2)

    def s1(self, a2, trunc, bratcutoff=None):
        ns = np.array([n for n in range(trunc-1, -1, -1)])
        gns = self.gn(a2, ns)
        # the ratio I1/I0 tends to 1 (inf/inf) very quickly
        brat = i1e(gns * self.R) / i0e(gns * self.R)
        if bratcutoff is not None:
            brat = np.where(gns * self.R < bratcutoff, brat, 1.0)
        return sum(16.0 / (np.pi**2 * self.R * gns * (2 * ns + 1)**2) * brat)

    def s2(self, a2, trunc):
        if self.zero_cache is None or len(self.zero_cache) < trunc:
            self.zero_cache = jn_zeros(0, trunc)
        lnks = self.ln(a2, self.zero_cache[::-1])
        return 8.0 * sum(
            np.tanh(lnks * self.H / 2)
            / (self.H * lnks * self.zero_cache[::-1]**2)
        )

    def s(self, a2, ntrunc, ktrunc):
        s1val = self.s1(a2, ntrunc)
        s2val = self.s2(a2, ktrunc)
        if np.any(np.logical_not(np.isfinite(s1val))):
            print("Nonfinite values in s1.")
            print(s1val)
        if np.any(np.logical_not(np.isfinite(s2val))):
            print("Nonfinite values in s2.")
            print(s2val)
        return s1val + s2val

    def computetransform(self, mech, T, ntrunc, ktrunc, dropnonreacting=False):
        #TODO interrogate the Mechanism object to choose an ideal method
        lams, R = la.eig(mech.getmatrix(self.charlen, T))
        if np.max(np.abs(lams.imag)) > 1e-8:
            print("Matrix has imaginary eigenvalues (irreversible).")
        lams, R = np.real_if_close(lams), np.real_if_close(R)
        mult = np.array([
            1.0 if lam == 0.0 else self.s(lam, ntrunc, ktrunc) for lam in lams
        ])
        A = np.dot(R, np.dot(np.diag(mult), la.inv(R)))
        if dropnonreacting:
            return A[:mech.getnumactive(), :]
        else:
            return A


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
        with open("template.f03", 'r') as f:
            program = str(f.read()).replace('!inline_init_here!', initlines)
    except FileNotFoundError:
        print("Could not find template file.")
        exit(1)

    if args.verbose:
        print("Writing output file 'mevlookup.f03'.")
    with open("mevlookup.f03", 'w') as f:
        f.write(program)


main()


