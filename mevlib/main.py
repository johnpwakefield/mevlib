#!/usr/bin/env python3


import click
from mevlib.maketable import make_table


def copystatement():
    click.echo("""
MEVLIB  Copyright (C) 2021  J Wakefield, A Lattanzi, J Capecelatro
This program comes with ABSOLUTELY NO WARRANTY.  This is free software, and
you are welcome to redistribute it under the GNU General Public License.
    """)


@click.command()
@click.option('--fmt', nargs=1, type=str, default=None)
@click.option('-d', '--discard-products/--keep-products', default=False)
@click.argument('src', type=click.Path(exists=True))
@click.argument('dst', type=click.Path(exists=False))
def maketable(src, dst, fmt=None, discard_products=False):
    copystatement()
    make_table(src, dst, fmt, not discard_products, True)


