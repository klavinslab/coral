'''Read and write DNA sequences.'''
import coral
import coral.constants.genbank
import re


def parse_genbank(handle):
    '''Parse DNA from a genbank file string. Uses BioPython and coerces to
    coral format.

    :param lines: Lines of the genbank file (e.g. f.read().split('\n')).
    :type lines: list
    :returns: DNA sequence.
    :rtype: coral.DNA

    '''

    read = handle.read()

    lines = read.split('\n')
    # Remove empty last line
    lines.pop()

    # Read in LOCUS information
    locus = lines.pop(0).split()
    name = locus[1]
    bp = locus[2]
    stranded = locus[4]
    topology = locus[5]
    date = locus[6]
    # Read in DEFINITION
    definition = lines.pop(0)
    description = re.search('DEFINITION  (.*)').group(1)
    # Read in FEATURES
    # Read in sequence after ORIGIN
