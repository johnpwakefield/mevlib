
Multistep Effectiveness Factor Lookup Library
==========================================


<add authors, funding acknowledgement, etc>



<todooverview>





Getting Started
------------------------------------------

There are five steps required to integrate <libname> into your existing code:

 1. Generate a data file containing precomputed values throughout the
    temperature range using `maketable.py`.

 2. Compile `efactors.f90` through your own build process (or build
    `efactors.mod` using the included makefile).

 3. Initialize the `efactors` object.  (Because object orientation in FORTRAN
    is patchy, the constructor `efactors%build` will need to be evoked
    explicitly.

 4. Call `efactors%getmev` to obtain a multistage effectiveness vector at the
    local temperature and free stream concentration.

 5. Call `efactors%destroy` to deallocate the memory claimed by
    `efactors%build`.


If you are using <list of some libs>, then all you have to do is turn it on.

<list of libs and info on configurations>







Developing with MEVLib
------------------------------------------


from the project directory:
python3 -m venv virtualenv
source ./virtualenv/bin/activate
python3 setup.py develop
# do stuff
deactivate






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



