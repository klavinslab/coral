try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension
import numpy

config = {
    'description': 'pymbt',
    'author': 'Nick Bolten',
    'url': 'https://github.com/klavinslab/pymbt',
    'download_url': 'https://github.com/klavinslab/pymbt.git',
    'author_email': 'nbolten _at_ gmail',
    'version': '0.0.2',
    'install_requires': ['nose', 'numpy', 'biopython', 'intermine',
                         'requests'],
    'extras_require': {'plotting': ['matplotlib'],
                       'documentation': ['sphinx']},
    'packages': ['pymbt',
                 'pymbt.analysis',
                 'pymbt.analysis._sequence',
                 'pymbt.analysis._sequencing',
                 'pymbt.analysis._structure',
                 'pymbt.constants',
                 'pymbt.database',
                 'pymbt.design',
                 'pymbt.design._oligo_synthesis',
                 'pymbt.design._sequence_generation',
                 'pymbt.seqio',
                 'pymbt.reaction',
                 'pymbt.sequence'],
    'package_data': {'pymbt': ['pymbt/analysis/_sequencing/data/*']},
    'include_package_data': True,
    'scripts': [],
    'name': 'pymbt',
    'license': 'Copyright University of Washington'
}

setup(ext_modules=[Extension('pymbt.analysis._sequencing.calign',
                             ['pymbt/analysis/_sequencing/calign.c'],
                             include_dirs=[numpy.get_include()])],
      test_suite='nose.collector',
      **config)
