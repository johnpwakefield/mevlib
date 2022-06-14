

import sys
from os.path import splitext

import click
from mevlib.maketable import make_table


COPYSTATEMENT = """
MEVLIB  Copyright (C) 2021  J Wakefield, A Lattanzi, J Capecelatro
This program comes with ABSOLUTELY NO WARRANTY.  This is free software, and
you are welcome to redistribute it under the GNU General Public License.
"""


def determine_outname(infile, outfile):
    if outfile is None:
        base, ext = splitext(infile)
        return base
    else:
        return outfile



@click.command()
@click.option('-a', '--all', 'wall', is_flag=True, help=(
    "Write all available tables for FMT."
    ))
@click.option('-r', '--rate', 'wrate', is_flag=True, help=(
    "Write lookup table of reaction rate transforms."
    ))
@click.option('-i', '--ints', 'wints', is_flag=True, help=(
    "Write lookup table of single-step effectiveness factors as a function "
    "of Thiele modulus."
    ))
@click.option('-d', '--diag', 'wdiag', is_flag=True, help=(
    "Write lookup table of diagonalizations as a function of temperature."
    ))
@click.option('-l', '--full', 'wfull', is_flag=True, help=(
    "Write lookup table of Multistep Effectiveness Vector (MEV) transforms."
    ))
@click.option('-o', '--outname', type=str, default=None, help=(
    "Output file name.  By default, output file is the input file name with "
    "the appropriate suffix."
))
@click.option('-v/-q', '--verb/--quiet', default=True, help=(
    "Verbosity."
    ))
@click.argument('src', type=click.Path(exists=True))
@click.argument('fmt', type=str)
def maketable(
        src, fmt, verb=None,
        wall=None, wrate=None, wints=None, wdiag=None, wfull=None, outname=None
        ):
    """\b
    Computes lookup tables used for estimating diffusion limitations in porous
    heterogeneous catalyst particles in CFD-DEM simulations.  Input parameters
    describing catalyst geometry and material properties are specified in SRC
    and output tables are written in a format specified in FMT.  Options for
    FMT include:
        - 'f03' (OO FORTRAN)
        - 'f90' (FORTRAN)
        - 'mat' (Matlab, see modules/mevrates.m)
        - 'bin' (binary)
        - 'pkl' (pickle, see modules/mevrates.py)

    A manuscript is in progress describing this work; for citation information
    contact jwake(at)umich.edu.
    """
    click.echo(COPYSTATEMENT)
    outfile = determine_outname(src, outname)
    make_table(src, outfile, fmt, wall, wrate, wints, wdiag, wfull, verb)

    sys.exit(0)


@click.command()
def otherexec():
    click.echo("""
The main executable for this module is `mevtable'.
    """)
    return 1


