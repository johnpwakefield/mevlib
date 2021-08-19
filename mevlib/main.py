

import click
from mevlib.maketable import make_table


def copystatement():
    click.echo("""
MEVLIB  Copyright (C) 2021  J Wakefield, A Lattanzi, J Capecelatro
This program comes with ABSOLUTELY NO WARRANTY.  This is free software, and
you are welcome to redistribute it under the GNU General Public License.
    """)


@click.command()
@click.option('-a', '--all', 'wall', is_flag=True)
@click.option('-i', '--ints', 'wints', is_flag=True)
@click.option('-d', '--diag', 'wdiag', is_flag=True)
@click.option('-m', '--mevs', 'wmevs', is_flag=True)
@click.option('-g', '--augd', 'waugd', is_flag=True)
@click.option('-r', '--rate', 'wrate', is_flag=True)
@click.option('-v/-q', '--verb/--quiet', default=True)
@click.argument('src', type=click.Path(exists=True))
@click.argument('fmt', type=str)
def maketable(
        src, fmt, verb=None,
        wall=None, wints=None, wdiag=None, wmevs=None, waugd=None, wrate=None
        ):
    copystatement()
    make_table(src, fmt, wall, wints, wdiag, wmevs, waugd, wrate, verb)


@click.command()
def otherexec():
    click.echo("""
The main executable for this module is `mevtable'.
    """)
    return 1


