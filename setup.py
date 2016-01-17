# Check python versions
import sys
if sys.version_info.major > 2:
    print('Coral is currently compatible only with Python 2.')
    sys.exit(1)

try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

import numpy

config = {
    'description': 'coral',
    'author': 'Nick Bolten',
    'url': 'https://github.com/klavinslab/coral',
    'download_url': 'https://github.com/klavinslab/coral.git',
    'author_email': 'nbolten _at_ gmail',
    'version': '0.3.2',
    'install_requires': ['numpy', 'biopython'],
    'extras_require': {'plotting': ['matplotlib'],
                       'yeastdatabases': ['intermine', 'requests'],
                       'documentation': ['sphinx']},
    'packages': ['coral',
                 'coral.analysis',
                 'coral.analysis._sequence',
                 'coral.analysis._sequencing',
                 'coral.analysis._structure',
                 'coral.constants',
                 'coral.database',
                 'coral.design',
                 'coral.design._oligo_synthesis',
                 'coral.design._sequence_generation',
                 'coral.seqio',
                 'coral.reaction',
                 'coral.sequence'],
    'package_data': {'coral': ['coral/analysis/_sequencing/data/*',
                               'coral/sequence/d3-plasmid.js']},
    'include_package_data': True,
    'scripts': [],
    'name': 'coral',
    'license': 'MIT',
    'classifiers': ['Programming Language :: Python',
                    'Programming Language :: Python 2.7',
                    'Programming Language :: Python :: 2 :: Only']
}

seq_extension = Extension('coral.analysis._sequencing.calign',
                          ['coral/analysis/_sequencing/calign.c'],
                          include_dirs=[numpy.get_include()])
EXTENSIONS = [seq_extension]

setup(ext_modules=EXTENSIONS,
      test_suite='nose.collector',
      **config)
