
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











Notes on Programming Style
------------------------------------------

Docstrings are a work in progress; reshuffling some code among modules may be
required to make module groupings make more sense.

All files in this project should adhere to PEP8 except:

  - W391 (avoids confusion with W292 on different text editors)
  - E306 (avoids awkward spacing and is better than violating E731)
  - E302 and E305 (two lines are used to denote different logical groupings of
    functions and classes)




TODO
------------------------------------------

  - migrate verification cases to new format
  - add other shapes to volrenders
  - write better documentation / tutorial
  - verify license decision with collaborators
  - add authorship page
  - add tests (pytest/tox, use pytest-cov)
  - add and sign versions
  - docstrings
  - add travis CI for pull requests
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
  - make sure this file obeys rst syntax, add a picture



