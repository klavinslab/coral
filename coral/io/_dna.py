'''Read and write DNA sequences.'''
import coral as cr
import csv
import os
from . import parsers
from . import writers
from .exceptions import UnsupportedFileError


class PrimerAnnotationError(ValueError):
    '''Raise if primer has insufficient annotation to be written to file.'''


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
    abi_exts = ['.abi', '.ab1']

    with open(path) as f:
        if ext in genbank_exts:
            return parsers.parse_genbank(f)
        elif ext in fasta_exts:
            return parsers.parse_fasta(f)
        elif ext in abi_exts:
            return cr.DNA(parsers.ABI(f).seq)
        else:
            raise UnsupportedFileError('File format not recognized.')


def read_dnas(directory):
    '''Read all DNA sequences files in a directory.

    :param directory: Path to directory containing sequencing files.
    :type directory: str
    :returns: A list of DNA sequences.
    :rtype: coral.DNA list

    '''
    dirfiles = os.listdir(directory)
    sequences = []
    for dirfile in dirfiles:
        path = os.path.join(directory, dirfile)
        try:
            sequences.append(read_dna(path))
        except UnsupportedFileError:
            pass

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
