# coral
[![Build Status](https://travis-ci.org/klavinslab/coral.svg?branch=master)](https://travis-ci.org/klavinslab/coral)
[![Documentation Status](https://readthedocs.org/projects/coral/badge/?version=latest)](https://readthedocs.org/projects/coral/?badge=latest)

Coral: A synthetic DNA design library. Read the documentation at http://coral.readthedocs.org. Coral works with PyPy so long as a PyPy-compatible numpy is installed.

## Installation:

Most users:
```
pip install coral
```

To get the latest on git:

```
git clone https://github.com/klavinslab/coral.git
cd coral
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

| Package | Added functionality |
| --- | --- |
| `matplotlib` | plotting sequencing analysis |
| `cython`     | 300 times faster sequencing alignment |

###system:

| Package | Added functionality |
| --- | --- |
| `NuPack` | Structural analysis |
| `ViennaRNA` | Structural analysis |
