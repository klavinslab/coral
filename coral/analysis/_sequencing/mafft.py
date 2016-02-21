'''A Coral wrapper for the MAFFT command line multiple sequence aligner.'''
import coral
import os
import shutil
import subprocess
import tempfile


def MAFFT(sequences, gap_open=1.53, gap_extension=0.0, retree=2):
    '''A Coral wrapper for the MAFFT command line multiple sequence aligner.

    :param sequences: A list of sequences to align.
    :type sequences: List of homogeneous sequences (all DNA, or all RNA,
                     etc.)
    :param gap_open: --op (gap open) penalty in MAFFT cli.
    :type gap_open: float
    :param gap_extension: --ep (gap extension) penalty in MAFFT cli.
    :type gap_extension: float
    :param retree: Number of times to build the guide tree.
    :type retree: int

    '''
    arguments = ['mafft']
    arguments += ['--op', str(gap_open)]
    arguments += ['--ep', str(gap_extension)]
    arguments += ['--retree', str(retree)]
    arguments.append('input.fasta')
    tempdir = tempfile.mkdtemp()
    try:
        with open(os.path.join(tempdir, 'input.fasta'), 'w') as f:
            for i, sequence in enumerate(sequences):
                if hasattr(sequence, 'name'):
                    name = sequence.name
                else:
                    name = 'sequence{}'.format(i)
                f.write('>{}\n'.format(name))
                f.write(str(sequence) + '\n')
        process = subprocess.Popen(arguments, stdout=subprocess.PIPE,
                                   stderr=open(os.devnull, 'w'), cwd=tempdir)
        stdout = process.communicate()[0]
    finally:
        shutil.rmtree(tempdir)

    # Process stdout into something downstream process can use

    records = stdout.split('>')
    # First line is now blank
    records.pop(0)
    aligned_list = []
    for record in records:
        lines = record.split('\n')
        name = lines.pop(0)
        aligned_list.append(coral.DNA(''.join(lines)))

    return aligned_list
