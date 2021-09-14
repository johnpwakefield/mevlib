

.. image:: logo/mevlogo.png
    :alt: MEVlib Logo
    :align: right


Multistep Effectiveness Factor Lookup Library
==============================================================================

This package is certainly in its early stages.  We hope that others will find
it useful, but it unlikely this package will do everything you want it to out
of the box (though this is of course the goal).  If you intend to use it, feel
free to send an email to the package maintainers (jwake@umich.edu) with any
questions or concerns.




Installation
------------------------------------------

Because the intended audience of this package is small, it is unlikely we will
add it to a repository like PyPi.  This package may either be used in a virtual
environment (recommended) or installed.  For the former, a script has been
included to do this for you (assuming standard packages like venv and pip3 are
already installed); simply run ``source activate_venv.sh``.  The latter can be
done for a user-specific (PEP370) install ``pip3 install path/to/mevlib`` or a
system-wide install ``sudo -H pip3 install path/to/mevlib``.


Usage Example
------------------------------------------

To generate a Matlab lookup tables::

    ~ $ source activate_venv.sh

    ~ $ cd examples/

    ~ $ mevtable --rates example_file.sbl mat

After doing so one will find ``mevtable_rates.mat`` in the current working
directory containing lookup table data.


Integration into CFD Packages
------------------------------------------

This package can produce four different lookuptables:

``rate``
    A table mapping temperature to :math:`N_g \times N` rate matrices (assumes
    infinite Biot number)

``biotrate``
    A table mapping temperature to :math:`N_g \times N` rate matrices (assumes
    infinite Biot number) TODO not yet implemented

``ints``
    A table mapping Thiele moduli to single stage effectiveness factors
    (:math:`\lambda \mapsto \eta`)

``biotints``
    A table mapping a Thiele modulus and a Biot number to a single stage
    effectiveness factor TODO not yet implemented

``diag``
    A table mapping temperatures to diagonalizations without left eigenvectors
    (temperature to :math:`N + (N + 1) N_g` scalars)

``fulldiag``
    A table mapping temperatures to diagonalizations with left eigenvectors
    (temperature to :math:`N + (N + 1) N_g` scalars)



An implementation in a CFD code would consist of either: a ``rate`` table
(where we must have :math:`\psi = 1`) or

  - One of
      - ``fulldiag``
      - ``diag`` and a means of computing :math:`R^{-1}` online
      - a means of computing the entire diagonalization online
  - and one of
      - ``ints`` or ``biotints``
      - a means of computing single stage effectiveness factors online

For some irreversible, cascading reactions :math:`R^{-1}` may be very easy to
compute, making its inclusion in a lookup table suboptimal.  For simple shapes
(e.g.~spheres) maintaining a table of single stage effectiveness factors is
entirely unnecessary.  For small numbers of species (e.g.~2 gaseous and 1
solid) in complex geometries it is likely optimal to do diagonalization online
but use a table for single stage effectiveness factors.

We hope that including the easiest to use option (a single ``rate`` table) in
addition to the ability to mix and match for a particular application will
increase the usefulness of this package.


Input File Formats
------------------------------------------

At present only one input file format is supported.  An example is provided in
``examples/example_file.sbl``.





Output File Formats
------------------------------------------






.mat and Pickle Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^





Binary Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All numbers are little-endian.  The format descriptions below impose the limit
of 8  ascii characters for a species name, 128 gaseous and solid species, and a
maximum lookup table size of 32768.  (The single-stage effectiveness factor has
a 4 byte table length so we mostly stay aligned to 32 bit words.)

For ``rate`` files:

+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
|                         | Element                           | Size (bytes)      | Data Type          |                      |
+=========================+===================================+===================+====================+======================+
|                         | Magic Number                      | 6                 | literal            | ``0x3ccf89cc873e``   |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Array Ordering                    | 1                 | 'F' or 'C' (ascii) | only 'F' implemented |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Version                           | 1                 | integer            | only ``0x01`` so far |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Number of Gases :math:`N_g`       | 1                 | integer            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Number of Solids :math:`N_s`      | 1                 | integer            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Size of Lookup Table :math:`M`    | 2                 | integer            |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
| Repeat                  | Species Identifier                | 8                 | ascii              |                      |
| :math:`N = N_g + N_s`   |                                   |                   |                    |                      |
| times                   |                                   |                   |                    |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
| Repeat :math:`M` times  | Temperature                       | 8                 | double             |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
| Repeat :math:`M` times  | Rate Matrix                       | :math:`8 N N_g`   | doubles            |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+

For ``ints`` files:

+-------------------------+-----------------------------------+--------------+--------------------+----------------------+
|                         | Element                           | Size (bytes) | Data Type          |                      |
+=========================+===================================+==============+====================+======================+
|                         | Magic Number                      | 6            | literal            | ``0xe288ab756456``   |
|                         +-----------------------------------+--------------+--------------------+----------------------+
|                         | Array Ordering                    | 1            | 'F' or 'C' (ascii) | only 'F' implemented |
|                         +-----------------------------------+--------------+--------------------+----------------------+
|                         | Version                           | 1            | integer            | only ``0x01`` so far |
|                         +-----------------------------------+--------------+--------------------+----------------------+
|                         | Number of Theile Moduli :math:`M` | 4            | integer            |                      |
+-------------------------+-----------------------------------+--------------+--------------------+----------------------+
| Repeat :math:`M` times  | Theile Modulus                    | 8            | double             |                      |
+-------------------------+-----------------------------------+--------------+--------------------+----------------------+
| Repeat :math:`M` times  | Single-Step Effectiveness Factor  | 8            | double             |                      |
+-------------------------+-----------------------------------+--------------+--------------------+----------------------+

For ``diag`` files:

+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
|                         | Element                           | Size (bytes)      | Data Type          |                      |
+=========================+===================================+===================+====================+======================+
|                         | Magic Number                      | 6                 | literal            | ``0x42523d52ce9b``   |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Array Ordering                    | 1                 | 'F' or 'C' (ascii) | only 'F' implemented |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Version                           | 1                 | integer            | only ``0x01`` so far |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Number of Gases :math:`N_g`       | 1                 | integer            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Number of Solids :math:`N_s`      | 1                 | integer            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Size of Lookup Table :math:`M`    | 2                 | integer            |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
| Repeat                  | Species Identifier                | 8                 | ascii              |                      |
| :math:`N = N_g + N_s`   |                                   |                   |                    |                      |
| times                   |                                   |                   |                    |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
| Repeat :math:`M` times  | Temperature                       | 8                 | double             |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
| Repeat :math:`M` times  | :math:`\bar{\mathbf{D}} / L^2`    | :math:`8 N`       | doubles            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Eigenvalues                       | :math:`8 N_g`     | doubles            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | :math:`R`                         | :math:`8 N_g^2`   | doubles            | fortran order        |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | :math:`B_s R`                     | :math:`8 N_s N_g` | doubles            | fortran order        |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+

For ``fulldiag`` files:

+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
|                         | Element                           | Size (bytes)      | Data Type          |                      |
+=========================+===================================+===================+====================+======================+
|                         | Magic Number                      | 6                 | literal            | ``0x525e7b2d317d``   |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Array Ordering                    | 1                 | 'F' or 'C' (ascii) | only 'F' implemented |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Version                           | 1                 | integer            | only ``0x01`` so far |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Number of Gases :math:`N_g`       | 1                 | integer            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Number of Solids :math:`N_s`      | 1                 | integer            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Size of Lookup Table :math:`M`    | 2                 | integer            |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
| Repeat                  | Species Identifier                | 8                 | ascii              |                      |
| :math:`N = N_g + N_s`   |                                   |                   |                    |                      |
| times                   |                                   |                   |                    |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
| Repeat :math:`M` times  | Temperature                       | 8                 | double             |                      |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+
| Repeat :math:`M` times  | :math:`\bar{\mathbf{D}} / L^2`    | :math:`8 N`       | doubles            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | Eigenvalues                       | :math:`8 N_g`     | doubles            |                      |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | :math:`R`                         | :math:`8 N_g^2`   | doubles            | fortran order        |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | :math:`B_s R`                     | :math:`8 N_s N_g` | doubles            | fortran order        |
|                         +-----------------------------------+-------------------+--------------------+----------------------+
|                         | :math:`R^{-1}`                    | :math:`8 N_g^2`   | doubles            | fortran order        |
+-------------------------+-----------------------------------+-------------------+--------------------+----------------------+






Fortran Modules
------------------------------------------

To make this package easier to use, some output formats generate code that can
be called from a containing simulation without any complex dependencies or data
files.  These are generated by using `f03` and `f90` as output formats.  `f90`
generates a data type `MEVData` containing the lookup table and a collection of
related functions whereas `f03` generates an class.  These modules contain
`mevdata_getmev`, `mevdata_init`, and `mevdata_destroy`.  The initialization
and destruction functions have no imput arguments (other than the data
structure).  `mevdata_getmev` requires a vector of free stream concentrations
and a temperature.






Note on Naming of the MEVlib Package
------------------------------------------

When referred to as a Python package `mevlib` is written in all lowercase to be
consistent with Python conventions.  As a project it is referred to as MEVlib.
For example, 'Effectiveness factors were computed with MEVlib.' and 'The mevlib
package is not available on PyPi.' are both correct.



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
  - write better documentation / tutorial
  - add tests (pytest/tox, use pytest-cov)
  - add and sign versions
  - ensure docstring coverage
  - several of the scripts in "scripts" are actually tests, make these into
    tests and move them into the tests directory
  - many scripts should maybe be "examples"
  - the options.py file referenced by scripts is kinda wonky
  - make sure this file obeys rst syntax


Attribution
------------------------------------------

This library/tool was written by John Wakefield (jwake@umich.edu) in
collaboration with Aaron Lattanzi, Brennan Pecha, Peter Ciesielski, and Jesse
Capacelatro.

For imformation on citing this paper contact jwake@umich.edu.

This software package was developed based upon funding from the Alliance for
Sustainable Energy, LLC, Managing and Operating Contractor for the National
Renewable Energy Laboratory for the U.S.  Department of Energy.

