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
    'version': '0.0.2',
    'install_requires': ['nose', 'numpy', 'biopython', 'intermine',
                         'requests'],
    'extras_require': {'plotting': ['matplotlib'],
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
    'package_data': {'coral': ['coral/analysis/_sequencing/data/*']},
    'include_package_data': True,
    'scripts': [],
    'name': 'coral',
    'license': 'Copyright University of Washington'
}

setup(ext_modules=[Extension('coral.analysis._sequencing.calign',
                             ['coral/analysis/_sequencing/calign.c'],
                             include_dirs=[numpy.get_include()])],
      test_suite='nose.collector',
      **config)
