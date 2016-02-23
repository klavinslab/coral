'''Read and write DNA sequences.'''
import csv
import os
from . import parsers
from . import writers


class PrimerAnnotationError(ValueError):
    pass


class FeatureNameError(ValueError):
    pass


def read_dna(path):
    '''Read a single DNA sequence from file.

    :param path: Full path to input file.
    :type path: str
    :returns: DNA sequence.
    :rtype: coral.DNA

    '''
    filename, ext = os.path.splitext(os.path.split(path)[-1])

    genbank_exts = ['.gb', '.ape']
    fasta_exts = ['.fasta', '.fa', '.fsa', '.seq']

    with open(path) as f:
        if any([ext == extension for extension in genbank_exts]):
            parser = parsers.parse_genbank
        elif any([ext == extension for extension in fasta_exts]):
            parser = parsers.parse_fasta
        else:
            raise ValueError('File format not recognized.')

        return parser(f)


def read_sequencing(directory):
    '''Read .seq and .abi/.ab1 results files from a dir.

    :param directory: Path to directory containing sequencing files.
    :type directory: str
    :returns: A list of DNA sequences.
    :rtype: coral.DNA list

    '''
    dirfiles = os.listdir(directory)
    seq_exts = ['.seq', '.abi', '.ab1']
    # Exclude files that aren't sequencing results
    seq_paths = [x for x in dirfiles if os.path.splitext(x)[1] in seq_exts]
    paths = [os.path.join(directory, x) for x in seq_paths]
    sequences = [read_dna(x) for x in paths]

    return sequences


def write_dna(dna, path):
    '''Write DNA to a file (genbank or fasta).

    :param dna: DNA sequence to write to file
    :type dna: coral.DNA
    :param path: file path to write. Has to be genbank or fasta file.
    :type path: str

    '''
    # Check if path filetype is valid, remember for later
    ext = os.path.splitext(path)[1]
    with open(path, 'w') as f:
        if ext == '.gb' or ext == '.ape':
            writer = writers.write_genbank
        elif ext == '.fa' or ext == '.fasta':
            writer = writers.write_fasta
        else:
            raise ValueError('Only genbank or fasta files are supported.')

        writer(dna, f)


def write_primers(primer_list, path, names=None, notes=None):
    '''Write a list of primers out to a csv file. The first three columns are
    compatible with the current IDT order form (name, sequence, notes). By
    default there are no notes, which is an optional parameter.

    :param primer_list: A list of primers.
    :type primer_list: coral.Primer list
    :param path: A path to the csv you want to write.
    :type path: str
    :param names: A list of strings to name each oligo. Must be the same length
                  as the primer_list.
    :type names: str list
    :param notes: A list of strings to provide a note for each oligo. Must be
                  the same length as the primer_list.
    :type notes: str list

    '''
    # Check for notes and names having the right length, apply them to primers
    if names is not None:
        if len(names) != len(primer_list):
            names_msg = 'Mismatch in number of notes and primers.'
            raise PrimerAnnotationError(names_msg)
        for i, name in enumerate(names):
            primer_list[i].name = name
    if notes is not None:
        if len(notes) != len(primer_list):
            notes_msg = 'Mismatch in number of notes and primers.'
            raise PrimerAnnotationError(notes_msg)
        for i, note in enumerate(notes):
            primer_list[i].note = note

    # Write to csv
    with open(path, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['name', 'sequence', 'notes'])
        for primer in primer_list:
            string_rep = str(primer.overhang).lower() + str(primer.anneal)
            writer.writerow([primer.name, string_rep, primer.note])
