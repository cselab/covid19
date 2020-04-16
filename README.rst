Bayesian Inference of Epidemic Models
---------------------------------------

coming soon...


Project structure
=================

- `epidemics` - The Python module and files.
- `epidemics/data` - Python files for fetch and preprocessing data; data files.
- `epidemics/data/files` - Data files.
- `epidemics/data/downloads` - Raw files downloaded by the Python scripts.
- `epidemics/data/cache` - Files processed by the Python scripts.

Troubleshooting
===============

If the changes in the code are not reflected in the results, try erasing `epidemics/data/cache` and `epidemica/data/downloads` folders.
The cache decorators attempt to detect changes in the code, but may not always be successful.
