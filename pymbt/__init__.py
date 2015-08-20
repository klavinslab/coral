'''python molecular biology tool, a cloning-oriented DNA design module.'''
try:
    import matplotlib
except ImportError:
    print "Failed to import matplotlib. Plotting sequencing won't work."

from . import analysis
from . import constants
from . import database
from . import design
from . import reaction
from . import seqio
from .sequence import DNA
from .sequence import RNA
from .sequence import Peptide
from .sequence import Primer
from .sequence import RestrictionSite
from .sequence import Feature
