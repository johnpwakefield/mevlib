#!/usr/bin/env bash


mevtable ../examples/example_file.sbl mevlookup.f90
gfortran -fPIC -c mevlookup.f90
f2py3 -c --fcompiler=gfortran -I. mevlookup.o -m getvals circfort_wrapper.f90


