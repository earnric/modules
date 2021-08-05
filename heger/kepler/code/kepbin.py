"""
Impurt kepler library.

There could be different versions for JMZ, NBUNZ, ...
"""

from ._build import _BuildKepler

# below we have to compile and import different versions
# maybe have a ..setup module.

_BuildKepler().run()

from . import _kepler as _kepler
