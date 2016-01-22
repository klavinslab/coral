# -*- coding: utf-8
'''Wrapper for NUPACK 3.0.'''
import multiprocessing
import numpy as np
import os
import re
import subprocess
import time
from coral.utils import tempdirs


class NUPACK(object):
    '''Run NUPACK functions on sequences.'''
    def __init__(self, nupack_home=None):
        '''
        :param nupack_home: NUPACK home dir. If the NUPACK commands aren't in
                            your path and the NUPACKHOME environment variable
                            isn't set, you can manually specify the NUPACK
                            directory here (the directory that contains bin/).
        :type nupack_home: str

        '''
        # Figure out where the NUPACK executables are
        if nupack_home is not None:
            self._nupack_home = nupack_home
        else:
            try:
                self._nupack_home = os.environ['NUPACKHOME']
            except:
                # FIXME: Test to see if it's in PATH
                pass

        # Initialize empty temp dir location
        self._tempdir = ''

    @tempdirs.tempdir
    def pfunc(self, strand, temp=37.0, pseudo=False, material=None,
              dangles='some'):
        '''Compute the partition function for an ordered complex of strands.
        Runs the \'pfunc\' command.

        :param strand: Strand on which to run pfunc. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: list
        :param temp: Temperature setting for the computation. Negative values
                     are not allowed.
        :type temp: float
        :param pseudo: Enable pseudoknots.
        :type pseudo: bool
        :param material: The material setting to use in the computation. If set
                         to None (the default), the material type is inferred
                         from the strands. Other settings available: 'dna' for
                         DNA parameters, 'rna' for RNA (1995) parameters, and
                         'rna1999' for the RNA 1999 parameters.
        :type material: str
        :param dangles: How to treat dangles in the computation. From the
                        user guide: For \'none\': Dangle energies are ignored.
                        For \'some\': \'A dangle energy is incorporated for
                        each unpaired base flanking a duplex\'. For 'all': all
                        dangle energy is considered.
        :type dangles: str
        :returns: Dictionary with the following key:value pairs: 'free_energy':
                  free energy, 'pfunc': partition function.
        :rtype: dict

        '''
        # Check the inputs
        if material is None:
            material = strand.material

        # Set up command flags
        cmd_args = []
        cmd_args += ['-T', temp]
        cmd_args += ['-dangles', dangles]
        cmd_args += ['-material', material]
        if pseudo:
            cmd_args.append('-pseudo')

        # Set up the input file
        with open(os.path.join(self._tempdir, 'pfunc.in'), 'w') as f:
            f.write(self._pfunc_file(strand)[0])

        # Run the command
        pfunc_stdout = self._run('pfunc', 'pfunc', cmd_args).split('\n')

        output = {}
        output['free_energy'] = float(pfunc_stdout[-3])
        output['pfunc'] = float(pfunc_stdout[-2])

        return output

    @tempdirs.tempdir
    def pfunc_multi(self, strands, permutation=None, temp=37.0, pseudo=False,
                    material=None, dangles='some'):
        '''Compute the partition function for an ordered complex of strands.
        Runs the \'pfunc\' command.

        :param strands: List of strands to use as inputs to pfunc -multi.
        :type strands: list
        :param permutation: The circular permutation of strands to test in
                            complex. e.g. to test in the order that was input
                            for 4 strands, the permutation would be [1,2,3,4].
                            If set to None, defaults to the order of the
                            input strands.
        :type permutation: list
        :param temp: Temperature setting for the computation. Negative values
                     are not allowed.
        :type temp: float
        :param pseudo: Enable pseudoknots.
        :type pseudo: bool
        :param material: The material setting to use in the computation. If set
                         to None (the default), the material type is inferred
                         from the strands. Other settings available: 'dna' for
                         DNA parameters, 'rna' for RNA (1995) parameters, and
                         'rna1999' for the RNA 1999 parameters.
        :type material: str
        :param dangles: How to treat dangles in the computation. From the
                        user guide: For \'none\': Dangle energies are ignored.
                        For \'some\': \'A dangle energy is incorporated for
                        each unpaired base flanking a duplex\'. For 'all': all
                        dangle energy is considered.
        :type dangles: str
        :returns: Dictionary with the following key:value pairs: 'energy':
                  free energy, 'pfunc': partition function.
        :rtype: dict

        '''
        # Check the inputs
        self._same_material(strands)

        if material is None:
            material = strands[0].material

        if permutation is None:
            permutation = range(1, len(strands) + 1)

        # Set up command flags
        cmd_args = []
        cmd_args.append('-multi')
        cmd_args += ['-T', temp]
        cmd_args += ['-dangles', dangles]
        cmd_args += ['-material', material]
        if pseudo:
            cmd_args.append('-pseudo')

        # Set up the input file
        with open(os.path.join(self._tempdir, 'pfunc.in'), 'w') as f:
            f.write('\n'.join(self._pfunc_file_multi(strands, permutation)))

        # Run the command
        pfunc_output = self._run('pfunc', 'pfunc', cmd_args).split('\n')

        output = {}
        output['free_energy'] = float(pfunc_output[-3])
        output['pfunc'] = float(pfunc_output[-2])

        return output

    @tempdirs.tempdir
    def pairs(self, strand, temp=37.0, pseudo=False, material=None,
              dangles='some', cutoff=0.001):
        '''Compute the pair probabilities for an ordered complex of strands.
        Runs the \'pairs\' command.

        :param strand: Strand on which to run pairs. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: list
        :param temp: Temperature setting for the computation. Negative values
                     are not allowed.
        :type temp: float
        :param pseudo: Enable pseudoknots.
        :type pseudo: bool
        :param material: The material setting to use in the computation. If set
                         to None (the default), the material type is inferred
                         from the strands. Other settings available: 'dna' for
                         DNA parameters, 'rna' for RNA (1995) parameters, and
                         'rna1999' for the RNA 1999 parameters.
        :type material: str
        :param dangles: How to treat dangles in the computation. From the
                        user guide: For \'none\': Dangle energies are ignored.
                        For \'some\': \'A dangle energy is incorporated for
                        each unpaired base flanking a duplex\'. For 'all': all
                        dangle energy is considered.
        :type dangles: str
        :param cutoff: Only probabilities above this cutoff appear in the
                       output.
        :type cutoff: float
        :returns: The probability matrix, where the (i, j)th entry
                  is the probability that base i is bound to base j. The matrix
                  is augmented (it's N+1 by N+1, where N is the number of bases
                  in the sequence) with an (N+1)th column containing the
                  probability that each base is unpaired.
        :rtype: numpy.array

        '''
        # Check the inputs
        if material is None:
            material = strand.material

        # Set up command flags
        cmd_args = []
        cmd_args += ['-T', temp]
        cmd_args += ['-dangles', dangles]
        cmd_args += ['-material', material]
        if pseudo:
            cmd_args.append('-pseudo')

        # Set up the input file
        with open(os.path.join(self._tempdir, 'pairs.in'), 'w') as f:
            f.write(self._pfunc_file(strand)[0])

        # Run the command. There's no STDOUT.
        self._run('pairs', 'pairs', cmd_args).split('\n')

        # Read the output from file
        with open(os.path.join(self._tempdir, 'pairs.ppairs')) as g:
            ppairs = g.read()
        data = re.search('\n\n\d*\n(.*)', ppairs, flags=re.DOTALL).group(1)
        N = len(strand)
        prob_matrix = np.zeros((N, N + 1))
        data_lines = data.split('\n')
        # Remove the last line (empty)
        data_lines.pop()
        # Convert into probability matrix
        for line in data_lines:
            i, j, prob = line.split('\t')
            prob_matrix[int(i) - 1, int(j) - 1] = float(prob)

        return prob_matrix

    @tempdirs.tempdir
    def pairs_multi(self, strands, permutation=None, temp=37.0, pseudo=False,
                    material=None, dangles='some', cutoff=0.001):
        '''Compute the pair probabilities for an ordered complex of strands.
        Runs the \'pfunc\' command.

        :param strands: List of strands to use as inputs to pfunc -multi.
        :type strands: list
        :param permutation: The circular permutation of strands to test in
                            complex. e.g. to test in the order that was input
                            for 4 strands, the permutation would be [1,2,3,4].
                            If set to None, defaults to the order of the
                            input strands.
        :type permutation: list

        :param temp: Temperature setting for the computation. Negative values
                     are not allowed.
        :type temp: float
        :param pseudo: Enable pseudoknots.
        :type pseudo: bool
        :param material: The material setting to use in the computation. If set
                         to None (the default), the material type is inferred
                         from the strands. Other settings available: 'dna' for
                         DNA parameters, 'rna' for RNA (1995) parameters, and
                         'rna1999' for the RNA 1999 parameters.
        :type material: str
        :param dangles: How to treat dangles in the computation. From the
                        user guide: For \'none\': Dangle energies are ignored.
                        For \'some\': \'A dangle energy is incorporated for
                        each unpaired base flanking a duplex\'. For 'all': all
                        dangle energy is considered.
        :type dangles: str
        :param cutoff: Only probabilities above this cutoff appear in the
                       output.
        :type cutoff: float
        :returns: Two probability matrices: The probability matrix as in the
                  pairs method (but with a dimension equal to the sum of the
                  lengths of the sequences in the permutation), and a similar
                  probability matrix where strands of the same species are
                  considered to be indistinguishable.
        :rtype: list

        '''
        # Check the inputs
        self._same_material(strands)

        if material is None:
            material = strands[0].material

        if permutation is None:
            permutation = range(1, len(strands) + 1)

        # Set up command flags
        cmd_args = []
        cmd_args += ['-T', temp]
        cmd_args += ['-dangles', dangles]
        cmd_args += ['-material', material]
        cmd_args.append('-multi')
        if pseudo:
            cmd_args.append('-pseudo')

        # Set up the input file
        with open(os.path.join(self._tempdir, 'pairs.in'), 'w') as f:
            f.write('\n'.join(self._pfunc_file_multi(strands, permutation)))

        # Run the command
        self._run('pairs', 'pairs', cmd_args).split('\n')

        # Read the output from file
        N = sum([len(s) for s in strands])
        matrices = []
        for mat_type in ['ppairs', 'epairs']:
            with open(os.path.join(self._tempdir, 'pairs.' + mat_type)) as g:
                data = g.read()

            probs = re.search('\n\n\d*\n(.*)', data, flags=re.DOTALL).group(1)
            lines = probs.split('\n')
            # Remove the last line (empty)
            lines.pop()
            prob_matrix = np.zeros((N, N + 1))
            # Convert into probability matrix
            for line in lines:
                split = line.split('\t')
                i = int(split[0]) - 1
                j = int(split[1]) - 1
                prob = float(split[2])
                prob_matrix[i, j] = prob
            matrices.append(prob_matrix)

        return matrices

    @tempdirs.tempdir
    def mfe(self, strand, temp=37.0, pseudo=False, material=None,
            dangles='some', degenerate=False):
        '''Compute the MFE for an ordered complex of strands. Runs the \'mfe\'
        command.

        :param strand: Strand on which to run mfe. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: coral.DNA or coral.RNA
        :param temp: Temperature setting for the computation. Negative values
                     are not allowed.
        :type temp: float
        :param pseudo: Enable pseudoknots.
        :type pseudo: bool
        :param material: The material setting to use in the computation. If set
                         to None (the default), the material type is inferred
                         from the strands. Other settings available: 'dna' for
                         DNA parameters, 'rna' for RNA (1995) parameters, and
                         'rna1999' for the RNA 1999 parameters.
        :type material: str
        :param dangles: How to treat dangles in the computation. From the
                        user guide: For \'none\': Dangle energies are ignored.
                        For \'some\': \'A dangle energy is incorporated for
                        each unpaired base flanking a duplex\'. For 'all': all
                        dangle energy is considered.
        :type dangles: str
        :param degenerate: Setting to True will result in returning a list of
                           dictionaries associated with structures having the
                           same, minimal MFE value.
        :type degenerate: bool
        :returns: A dictionary with keys for 'mfe' (a float), 'dotbracket'
                  (dot-parens notation of the MFE structure), and 'pairlist'
                  (a pair list notation of the MFE structure). Note that the
                  pair list will be an empty list if the MFE is unstructured.
        :rtype: dict

        '''
        # Check the inputs
        if material is None:
            material = strand.material

        # Set up command flags
        cmd_args = []
        cmd_args += ['-T', temp]
        cmd_args += ['-dangles', dangles]
        cmd_args += ['-material', material]
        if pseudo:
            cmd_args.append('-pseudo')
        if degenerate:
            cmd_args.append('-degenerate')

        # Set up the input file
        with open(os.path.join(self._tempdir, 'mfe.in'), 'w') as f:
            f.write(self._pfunc_file(strand)[0])

        # Run the command. There's no STDOUT.
        self._run('mfe', 'mfe', cmd_args).split('\n')

        # Read the output from file
        with open(os.path.join(self._tempdir, 'mfe.mfe')) as g:
            data = g.read()

        structures = self._process_mfe(data)

        if degenerate:
            return structures
        else:
            return structures[0]

    @tempdirs.tempdir
    def mfe_multi(self, strands, permutation=None, temp=37.0, pseudo=False,
                  material=None, dangles='some', degenerate=False):
        '''Compute the MFE for an ordered complex of strands. Runs the \'mfe\'
        command.

        :param strands: Strands on which to run mfe. Strands must be either
                       coral.DNA or coral.RNA).
        :type strands: list
        :param temp: Temperature setting for the computation. Negative values
                     are not allowed.
        :type temp: float
        :param pseudo: Enable pseudoknots.
        :type pseudo: bool
        :param material: The material setting to use in the computation. If set
                         to None (the default), the material type is inferred
                         from the strands. Other settings available: 'dna' for
                         DNA parameters, 'rna' for RNA (1995) parameters, and
                         'rna1999' for the RNA 1999 parameters.
        :type material: str
        :param dangles: How to treat dangles in the computation. From the
                        user guide: For \'none\': Dangle energies are ignored.
                        For \'some\': \'A dangle energy is incorporated for
                        each unpaired base flanking a duplex\'. For 'all': all
                        dangle energy is considered.
        :type dangles: str
        :param degenerate: Setting to True will result in returning a list of
                           dictionaries associated with structures having the
                           same, minimal MFE value.
        :type degenerate: bool
        :returns: A dictionary with keys for 'mfe' (a float), 'dotbracket'
                  (dot-parens notation of the MFE structure), and 'pairlist'
                  (a pair list notation of the MFE structure). Note that the
                  pair list will be an empty list if the MFE is unstructured.
        :rtype: dict

        '''
        # Check the inputs
        self._same_material(strands)
        if material is None:
            material = strands[0].material

        if permutation is None:
            permutation = range(1, len(strands) + 1)

        # Set up command flags
        cmd_args = []
        cmd_args += ['-T', temp]
        cmd_args += ['-dangles', dangles]
        cmd_args += ['-material', material]
        if pseudo:
            cmd_args.append('-pseudo')
        if degenerate:
            cmd_args.append('-degenerate')
        cmd_args.append('-multi')

        # Set up the input file
        with open(os.path.join(self._tempdir, 'mfe.in'), 'w') as f:
            f.write('\n'.join(self._pfunc_file_multi(strands, permutation)))

        # Run the command. There's no STDOUT.
        self._run('mfe', 'mfe', cmd_args).split('\n')

        # Read the output from file
        with open(os.path.join(self._tempdir, 'mfe.mfe')) as g:
            data = g.read()

        structures = self._process_mfe(data)
        if degenerate:
            return structures
        else:
            return structures[0]

    # Helper methods for preparing command input files
    def _pfunc_file(self, strand):
        '''Prepares lines to write to file for pfunc command input.

        :param strand: Strand input (cr.DNA or cr.RNA).
        :type strand: cr.DNA or cr.DNA
        '''
        return [str(strand)]

    def _pfunc_file_multi(self, strands, permutation):
        '''Prepares lines to write to file for pfunc command input.

        :param strand: Strand input (cr.DNA or cr.RNA).
        :type strand: cr.DNA or cr.DNA
        :param permutation: Permutation (e.g. [1, 2, 3, 4]) of the type used
                            by pfunc_multi.
        :type permutation: list
        '''
        lines = []
        # Write the total number of distinct strands
        lines.append(str(len(strands)))
        # Write the distinct strands
        lines += [str(strand) for strand in strands]
        # Write the permutation
        lines.append(' '.join(str(p) for p in permutation))

        return lines

    # Helper methods for processing output files
    def _process_mfe(self, data):
        # Parse the output data
        # Find the text in between the large comment lines
        commentline = '\n% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n'
        # Find everything between two commentlines
        groups = data.split(commentline)
        # Everything before the comment line is notes about the command,
        # the last part is just a newline
        groups.pop(0)
        groups.pop()
        # The remainder is data
        output = []
        # Skip every other one (every 2nd match is empty lines)
        for group in groups[::2]:
            lines = group.split('\n')
            # Line 1 is the strand number (ignored)
            lines.pop(0)
            # Line 2 is the MFE
            mfe = float(lines.pop(0))
            # Line 3 is the dot-bracket structure
            dotbracket = lines.pop(0)
            # If there are any more lines, they are a pair list format
            # structure
            pairlist = []
            for line in lines:
                pair = line.split('\t')
                pairlist.append([int(pair[0]) - 1, int(pair[1]) - 1])
            output.append({'mfe': mfe, 'dotbracket': dotbracket,
                           'pairlist': pairlist})

        return output

    # Helper methods for repetitive tasks
    def _same_material(self, strands):
        if len(set([s.material for s in strands])) > 1:
            raise ValueError('Inputs must all be coral.DNA or all coral.RNA')

    def _run(self, command, prefix, cmd_args):
        arguments = [os.path.join(self._nupack_home, 'bin', command)]
        arguments += cmd_args
        arguments.append(prefix)

        arguments = [str(x) for x in arguments]
        process = subprocess.Popen(arguments, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, cwd=self._tempdir)
        output = process.communicate()[0]
        return output


def nupack_multi(seqs, material, cmd, arguments, report=True):
    '''Split Nupack commands over processors.

    :param inputs: List of sequences, same format as for coral.analysis.Nupack.
    :type inpus: list
    :param material: Input material: 'dna' or 'rna'.
    :type material: str
    :param cmd: Command: 'mfe', 'pairs', 'complexes', or 'concentrations'.
    :type cmd: str
    :param arguments: Arguments for the command.
    :type arguments: str
    :returns: A list of the same return value you would get from `cmd`.
    :rtype: list

    '''
    nupack_pool = multiprocessing.Pool()
    try:
        args = [{'seq': seq,
                 'cmd': cmd,
                 'material': material,
                 'arguments': arguments} for seq in seqs]
        nupack_iterator = nupack_pool.imap(run_nupack, args)
        total = len(seqs)
        msg = ' calculations complete.'
        passed = 4
        while report:
            completed = nupack_iterator._index
            if (completed == total):
                break
            else:
                if passed >= 4:
                    print '({0}/{1}) '.format(completed, total) + msg
                    passed = 0
                passed += 1
                time.sleep(1)
        multi_output = [x for x in nupack_iterator]
        nupack_pool.close()
        nupack_pool.join()
    except KeyboardInterrupt:
        nupack_pool.terminate()
        nupack_pool.close()
        raise KeyboardInterrupt

    return multi_output


def run_nupack(kwargs):
    '''Run picklable Nupack command.

    :param kwargs: keyword arguments to pass to Nupack as well as 'cmd'.
    :returns: Variable - whatever `cmd` returns.

    '''
    run = NUPACK(kwargs['seq'])
    output = getattr(run, kwargs['cmd'])(**kwargs['arguments'])
    return output
