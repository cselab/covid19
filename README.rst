Bayesian Inference of Epidemic Models
---------------------------------------

coming soon...


Project structure
=================

- ``epidemics`` - The Python module and files.
- ``epidemics/data`` - Python files for fetch and preprocessing data; data files.
- ``epidemics/data/files`` - Data files.
- ``epidemics/data/downloads`` - Raw files downloaded by the Python scripts.
- ``epidemics/data/cache`` - Files processed by the Python scripts.
- ``src/epidemics`` - The C++ code.
- ``src/epidemics/model/country`` - C++ country-level model implementations.
- ``src/epidemics/model/cantons`` - C++ canton-level model implementations.
- ``test/py`` - Python unit tests.
- ``test/cpp`` - C++ unit tests.


Compilation
===========

From the repository root folder do:

.. code-block:: bash

    git submodule update --init --recursive
    mkdir -p
    cd build
    cmake ..
    make

Tests
=====

To run the tests, run the following command (from the repository root):

.. code-block:: bash

    cd tests
    ./run_all_tests.sh

To run only Python tests, run ``cd tests/py && ./run.sh``.
To run only C++ tests, run ``cd build && ./libepidemics_unittests``.


Troubleshooting
===============

If the changes in the code are not reflected in the results, try erasing ``epidemics/data/cache`` and ``epidemica/data/downloads`` folders.
The cache decorators attempt to detect changes in the code, but may not always be successful.


Adding a new country-level model (C++)
======================================

Follow these steps to create a new C++ country-level model. The steps are shown on an example of creating a XYZ model from the existing SIR model.

1. Make a copy of ``src/epidemics/models/country/sir.h`` and name it ``xyz.h``.

2. Change the ``sir`` namespace to ``xyz``.

2. Update the ``Parameters`` struct. Update the number of derivatives in ``Element``.

3. Update the ``State`` struct: change the number of states in ``StateBase`` parent class, and customize named getters.

4. Make a copy of ``src/epidemics/models/country/sir.cpp`` and name it ``xyz.cpp``.

5. Change the ``sir`` namespace to ``xyz``.

6. Update the model: update the ``make_ad`` at the beginning to initialize AD variables, update formulae and the output.

7. Go to ``src/epidemics/country.cpp`` and make a function ``exportXYZ`` similar to ``exportSIR``.

8. Call ``exportXYZ`` from ``exportCountryModels``.

9. Add the new .cpp file to ``CMakeLists.txt``.

10. Create a ``test/py/test_country_xyz.py`` analoguous to ``test_country_sir.h`` and test your code. You may skip testing the derivatives, since AD should already be tested.

In the case AD does not support some operation, add it in ``src/epidemics/utils/autodiff.h``.
Create a test in ``test/cpp/test_autodiff.cpp``!
