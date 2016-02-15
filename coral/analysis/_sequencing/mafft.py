'''A Coral wrapper for the MAFFT command line multiple sequence aligner.'''
import coral
import os
import shutil
import subprocess
import tempfile


def MAFFT(sequences):
    '''A Coral wrapper for the MAFFT command line multiple sequence aligner.

    :param sequences: A list of sequences to align.
    :type sequences: List of homogeneous sequences (all DNA, or all RNA,
                     etc.)

    '''
    arguments = ['mafft', 'input.fasta']
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
