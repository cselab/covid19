#!/bin/bash

set -e

# First try to import libepidemics, then run tests. `unittest` ends up
# importing libepidemics multiple times. If the import fails, multiple errors
# appear. Apart from the first one, the others are super confusing. So, we
# separately try to import libepidemics only once.
PYTHONPATH=../..:../../build:$PYTHONPATH python3 -c "import libepidemics"
PYTHONPATH=../..:../../build:$PYTHONPATH python3 -m unittest "$@"
