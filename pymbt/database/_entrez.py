import os
from tempfile import mkdtemp
from Bio import Entrez
import pymbt.sequence


# FIXME: If a Genome is a data structure, it should be in the DNA sequence
# module. Then rename the genome acquirer after Entrez / NCBI, etc.
# FIXME: If a Genome is not a new special data structure, should be a function
# not a class unless there are more operations than fetch.
# TODO: docstring
# TODO: Figure out why reading in the DNA is so slow and if it can be sped up
# - MG1655 takes 30-60 seconds to process into memory and pbt.DNA.
# MG1655 id is 'U00096.3'
def fetch_genome(genome_id):
    '''Acquire a genome from Entrez

    '''
    # TODO: Can strandedness by found in fetched genome attributes?
    # TODO: skip read/write step?
    # Using a dummy email for now - does this violate NCBI guidelines?
    email = 'loremipsum@gmail.com'
    Entrez.email = email

    print 'Downloading Genome...'
    handle = Entrez.efetch(db='nucleotide', id=str(genome_id), rettype='gb',
                           retmode='text')
    print 'Genome Downloaded...'
    tmpfile = os.path.join(mkdtemp(), 'tmp.gb')
    with open(tmpfile, 'w') as f:
        f.write(handle.read())
    genome = pymbt.seqio.read_dna(tmpfile)

    return genome
