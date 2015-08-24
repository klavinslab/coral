# pymbt
[![Build Status](https://travis-ci.org/klavinslab/pymbt.svg?branch=master)](https://travis-ci.org/klavinslab/pymbt)
[![Documentation Status](https://readthedocs.org/projects/pymbt/badge/?version=latest)](https://readthedocs.org/projects/pymbt/?badge=latest)

PyMBT: Python molecular biology tools. Read the documentation at http://pymbt.readthedocs.org. Works with PyPy so long as a PyPy-compatible numpy is installed.

## Installation:

Most users:
```
pip install pymbt
```

To get the latest on git:

```
git clone https://github.com/klavinslab/pymbt.git
cd pymbt
pip install .
```

## Requirements:

###python (pip-compatible):

```
numpy
biopython
intermine
```

optional:
| `matplotlib` | plotting sequencing analysis |
| `cython`     | 300 times faster sequencing alignment |

###system:

| `NuPack` | Structural analysis |
| `ViennaRNA` | Structural analysis |
