
Multistep Effectiveness Factor Lookup Library
==========================================


<add authors, funding acknowledgement, etc>



<todooverview>





Installation
------------------------------------------

Because the intended audience of this package is small, it is unlikely we will
add it to a repository like PyPi.  This package may either be used in a virtual
environment (recommended) or installed.  For the former, a script has been
included to do this for you (assuming standard packages like venv and pip3 are
already installed); simply run `source activate_venv.sh`.  The latter can be
done for a user-specific (PEP370) install `pip3 install path/to/mevlib` or a
system-wide install `sudo -H pip3 install path/to/mevlib`.













Detailed Implementation Instructions
------------------------------------------










TODO
------------------------------------------

  - migrate verification cases to new format
  - add other shapes to volrenders
  - write better documentation / tutorial
  - verify license decision with collaborators
  - add authorship page
  - implement full test
        series -> fortran lookup tables -> fortran compiled
        -> call from python -> compare to other method
  - it is possible several of the files in this directory should be in the
    package directory; determine if this is the case and fix accordingly
  - several of the scripts in "scripts" are actually tests, make these into
    tests and move them into the tests directory
  - scripts should have a function as an entry point
  - many scripts should maybe be "examples"
  - the scripts.py file referenced by scripts is kinda wonky
  - move unit tests out of package directory and use unittests module for these



