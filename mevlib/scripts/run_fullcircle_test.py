#!/usr/bin/env python3


import os
import numpy.f2py as f2py


from mevlib.maketable import make_table



builddir = "build_circular_test"
fn = "../../examples/example_file.sbl"   #TODO change this to not be relative
outfile = builddir + '/' + "mevlookup.f03"


# run make_table
if not os.path.isdir(builddir):
    os.mkdir(builddir)
make_table(
    fn, outfile, verb=False, keepnonreacting=True
)

# compile mevlookup.mod
#TODO f2py does not support fortran 03; we will need to figure something else
# out here
with open(outfile, "r") as f:
    f2py.compile(
        f.read(),
        modulename=(builddir + '/' + "mevlookup.mod"),
        source_fn=(builddir + '/' + "mevlookup.f90")
    )


# compile fullcircletest (and link mevlookup)



# test mesh values and compute interpolation error




#TODO delete this when draft of test is complete
assert(False)



