#!/usr/bin/env python3


import setuptools


setuptools.setup(
    name="MEVLIB",
    version="0.1",
    author="John P Wakefield et al",
    author_email="jwake@umich.edu",
    description="""
Generates pre-computed maps from free stream concentration to Multistep
Effectiveness Vector (MEV) for use in CFD-DEM models for chemical reactions
(e.g. packed-bed catalysis/pyrolysis).
    """,
    packages=["mevlib"]
)


