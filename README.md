# coral
<img align="right" height="256" src="docs/coral_256.png">
[![Build Status](https://travis-ci.org/klavinslab/coral.svg?branch=master)](https://travis-ci.org/klavinslab/coral)
[![Documentation Status](https://readthedocs.org/projects/coral/badge/?version=latest)](https://readthedocs.org/projects/coral/?badge=latest)


Coral: Core tools for synthetic DNA design. Read the documentation at http://coral.readthedocs.org.

Coral encodes synthetic DNA design rules into its core sequence data types (`DNA`, `RNA`, and `Peptide`), enabling concise, dependable methods for automated DNA design.

Coral works with PyPy so long as a PyPy-compatible numpy is installed.

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

## License

MIT
