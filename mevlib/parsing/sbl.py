

from functools import reduce, partial
from math import floor

import pyparsing as pp

from mevlib.shapes import Cylinder
from mevlib.mechanisms import KnudsenSpecies, SolidSpecies, ArrheniusReaction


class SBLParsingException(Exception):
    pass


def parse_sensible(f, verb=False):
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
    pp_phase = keylist(["s", "g"])
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
        pp_symb + pp.Optional(pp_qstr) + pp.Optional(pp_phase)
        + pp.Optional(pp_real)
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

    # actually parse the file
    d = {
        h: {k.lower(): v for k, v in u}
        for h, u in parser.parseFile(f).asList()
    }

    # read series truncation or error tolerance
    if 'precision' in d:
        if verb:
            print("Found section '{}'.".format('precision'))
        precision = d['precision']
        if 'type' not in d['precision']:
            print("Precision type not specified.")
            raise SBLParsingException("Precision type not specified.")
        if d['precision']['type'] == 'trunc':
            if not ('ntrunc' in d['precision'] and 'ktrunc' in d['precision']):
                if 'trunc' in d['precision']:
                    if verb:
                        print(
                            "Single truncation given, using this for both"
                            "sums."
                        )
                    precision['ntrunc'] = d['precision']['trunc']
                    precision['ktrunc'] = d['precision']['trunc']
                else:
                    print("Missing truncation information.")
                    raise SBLParsingException(
                        "Missing truncation information."
                    )
        elif d['precision']['type'] == 'tol':
            raise NotImplementedError(
                "Convergence tolerance not implemented yet."
            )
    else:
        precision = None

    # read shape
    if 'shape' in d:
        if verb:
            print("Found section '{}'.".format('shape'))
        if 'type' not in d['shape']:
            print("Shape type not specified.")
            raise SBLParsingException("Shape type not specified.")
        if d['shape']['type'] == "cylinder":
            if 'radius' not in d['shape'] and 'diameter' not in d['shape']:
                print("Radius missing from shape dimensions.")
                raise SBLParsingException("Missing shape dimensions.")
            if 'radius' in d['shape'] and 'diameter' in d['shape']:
                print("Radius and diameter cannot both be specified.")
                raise SBLParsingException("Shape dimensions overspecified.")
            if 'height' not in d['shape']:
                print("Height missing from shape dimensions.")
                raise SBLParsingException("Missing shape dimensions.")
            if 'radius' in d['shape']:
                shape = Cylinder(d['shape']['radius'], d['shape']['height'])
            else:
                shape = Cylinder(
                    d['shape']['diameter'] / 2, d['shape']['height']
                )
        elif d['shape']['type'] == "box":
            raise NotImplementedError("Shape 'box' not implemented.")
        elif d['shape']['type'] == "sphere":
            raise NotImplementedError("Shape 'sphere' not implemented.")
        else:
            print("Shape '{}' not recognized.".format(d['shape']['type']))
            raise SBLParsingException(
                "Shape '{}' not recognized.".format(d['shape']['type'])
            )
    else:
        shape = None

    # read temperature ranges
    if 'temperature' in d:
        if verb:
            print("Found section '{}'.".format('temperature'))
        if 'type' not in d['temperature']:
            print("Type of temperature specification not present.")
            raise SBLParsingException("Missing type for temperature range.")
        if d['temperature']['type'] == "uniform":
            if (
                'start' not in d['temperature']
                or
                'stop' not in d['temperature']
            ):
                print("Start and stop of temperature range must be given.")
                raise SBLParsingException("Range parameters unspecified.")
            if (
                'num' not in d['temperature']
                and
                'step' not in d['temperature']
            ):
                print("Either number of points or spacing must be specified.")
                raise SBLParsingException("Range parameters unspecified.")
            if 'num' in d['temperature'] and 'step' in d['temperature']:
                print("Only one of 'num' or 'step' may be specified.")
                raise SBLParsingException("Range parameters overspecified.")
            s = d['temperature']['stop'] - d['temperature']['start']
            if 'num' in d['temperature']:
                num = d['temperature']['num']
            else:
                num = floor(s / d['temperature']['step'])
            step = s / (num - 1)
            temperatures = [
                d['temperature']['start'] + i * step for i in range(num)
            ]
        elif d['temperature']['type'] == "explicit":
            raise NotImplementedError("Explicit temperatures not implemented.")
        else:
            print((
                "Temperature type '{}' not recognized."
            ).format(d['temperature']['type']))
            raise SBLParsingException("Temperature type not recognized.")
    else:
        temperatures = None

    # read diffusion parameters
    if 'diffusion' in d:
        if verb:
            print("Found section '{}'.".format('diffusion'))
        if d['diffusion']['type'] == 'knudsen':
            knudsenparams = ['porediameter', 'voidage', 'tortuosity']
            for param in knudsenparams:
                if param not in d['diffusion'].keys():
                    print("Diffusion section missing '{}'.".format(param))
                    raise SBLParsingException("Missing diffusion parameters.")
            gasconstructor = partial(KnudsenSpecies, *[
                d['diffusion'][param] for param in knudsenparams
            ])
        else:
            print((
                "Did not recognize diffusion type '{}'."
            ).format(d['diffusion']['type']))
            raise SBLParsingException("Unrecognized diffusion type.")
    solidconstructor = SolidSpecies

    # make species list
    if 'species' in d:
        if verb:
            print("Found section '{}'.".format('species'))
        if 'table' not in d['species']:
            raise SBLParsingException("Missing species table.")
        name = d['species'].get('longname', 'default')
        if name == 'default':
            name == 'omitted'
        phase = d['species'].get('phase', 'default')
        if phase == 'default':
            phase == 'specified'
        molweight = d['species'].get('molweight', 'default')
        if molweight == 'default':
            molweight == 'specified'
        def checkrow(row):
            if (
                len(row)
                != (1 + (name == 'specified') + (phase == 'specified') + (molweight == 'specified'))
            ):
                print("Unexpected length of table row.")
                expected_syms = [
                    ['symbol']
                    + (['name'] if name == 'specified' else [])
                    + (['phase'] if phase == 'specified' else [])
                    + (['molweight'] if molweight == 'specified' else [])
                ]
                print("Expected {} columns: {}.".format(
                    len(expected_syms), ", ".join(expected_syms)
                ))
                raise SBLParsingException("Unexpected length of table row.")
        def getname(row):
            if name == 'omitted':
                return row[0]
            elif name == 'specified':
                return row[1]
            else:
                if verb:
                    print((
                        "Using chemical name '{}' for all species; this was "
                        "likely unintended. Try using 'DEFAULT', 'OMITTED', "
                        "or 'SPECIFIED'. "
                    ).format(name))
                return name
        def getphase(row):
            if phase == 'omitted':
                raise SBLParsingException("Cannot omit phase.")
            else:
                return row[2]
        def getmolweight(row):
            if molweight == 'omitted':
                raise SBLParsingException("Cannot ommit molecular weight.")
            elif molweight == 'specified':
                return row[3]
        def makespecies(row):
            checkrow(row)
            if getphase(row).strip().lower() == 'g':
                return gasconstructor(row[0], getname(row), getmolweight(row))
            elif getphase(row).strip().lower() == 's':
                return solidconstructor(row[0], getname(row), getmolweight(row))
            else:
                raise SBLParsingException("Phase must be 'g' or 's'.")
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
            raise SBLParsingException("Prefactor omitted.")
        if exponent == 'default' or exponent == 'omitted':
            exponent = None
        if actenergy == 'default':
            actenergy = 'specified'
        if actenergy == 'omitted':
            print("Activation energy may not be omitted.")
            raise SBLParsingException("Activation energy omitted.")
        if reftemp == 'default' or reftemp == 'omitted':
            reftemp = None
        expected_row_length = (
            2 + (prefactor == 'specified') + (exponent == 'specified')
            + (actenergy == 'specified') + (reftemp == 'specified')
        )
        def makereaction(row):
            if len(row) != expected_row_length:
                print("Unexpected row length in reaction table.")
                raise SBLParsingException("Unexpected row length.")
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
        raise SBLParsingException("Unrecognized reaction type.")
    reactions = [makereaction(row) for row in d['reactions']['table']]

    return precision, shape, temperatures, species, reactions


