# -*- coding: utf-8
'''Wrapper for ViennaRNA functions.'''
# Due to inconsistency of command inputs, each command is a function
import os
import re
from subprocess import Popen, PIPE, STDOUT
import numpy as np
from coral.utils import tempdirs


# TODO: Generic structure object to return from ViennaRNA, NUPACK classes
class ViennaRNA(object):

    def __init__(self):
        self._tempdir = ''

    @tempdirs.tempdir
    def cofold(self, strand1, strand2, temp=37.0, dangles=2, nolp=False,
               nogu=False, noclosinggu=False, constraints=None,
               canonicalbponly=False, partition=-1, pfscale=None, gquad=False):
        '''Run the RNAcofold command and retrieve the result in a dictionary.

        :param strand1: Strand 1 for running RNAcofold.
        :type strand1: coral.DNA or coral.RNA
        :param strand1: Strand 2 for running RNAcofold.
        :type strand2: coral.DNA or coral.RNA
        :param temp: Temperature at which to run the calculations.
        :type temp: float
        :param dangles: How to treat dangling end energies. Set to 0 to ignore
                        dangling ends. Set to 1 to limit unpaired bases to
                        at most one dangling end (default for MFE calc). Set to
                        2 (the default) to remove the limit in 1. Set to 3 to
                        allow coaxial stacking of adjacent helices in
                        .multi-loops
        :type dangles: int
        :param nolp: Produce structures without lonely pairs (isolated single
                     base pairs).
        :type nolp: bool
        :param nogu: Do not allow GU pairs.
        :type nogu: bool
        :param noclosinggu: Do not allow GU pairs at the end of helices.
        :type noclosinggu: bool
        :param constraints: Any structural constraints to use. Format is
                            defined at
                            http://www.tbi.univie.ac.at/RNA/RNAfold.1.html
        :type constraints: str
        :param canonicalbponly: Remove non-canonical base pairs from the
                                structure constraint (if applicable).
        :type canonicalbponly: bool
        :param partition: Calculates the partition function for the sequence.
        :type partition: int
        :param pfscale: Scaling factor for the partition function.
        :type pfScale: float
        :param gquad: Incorporate G-Quadruplex formation into the structure
                      prediction.
        :type gquad: bool

        :returns: Dictionary of calculated values, defaulting to values of MFE
                  ('mfe': float) and dotbracket structure ('dotbracket': str).
                  More keys are added depending on keyword arguments.
        :rtype: dict

        '''
        cmd_args = []
        cmd_kwargs = {'--temp=': str(temp)}
        cmd_kwargs['--dangles='] = dangles
        if nolp:
            cmd_args.append('--noLP')
        if nogu:
            cmd_args.append('--noGU')
        if noclosinggu:
            cmd_args.append('--noClosingGU')
        if constraints is not None:
            cmd_args.append('--constraint')
            if canonicalbponly:
                cmd_args.append('--canonicalBPonly')
        if partition:
            cmd_args.append('--partfunc')
        if pfscale is not None:
            cmd_kwargs['pfScale'] = float(pfscale)
        if gquad:
            cmd_args.append('--gquad')

        inputs = ['>strands\n{}&{}'.format(str(strand1), str(strand2))]
        if constraints is not None:
            inputs.append(constraints)

        rnafold_output = self._run('RNAcofold', inputs, cmd_args, cmd_kwargs)

        # Process the output
        output = {}
        lines = rnafold_output.splitlines()
        # Line 1 is the name of the sequence input, line 2 is the sequence
        lines.pop(0)
        lines.pop(0)
        # Line 3 is the dotbracket + mfe for strand1
        line3 = lines.pop(0)
        output['dotbracket'] = self._lparse(line3, '^(.*) \(')
        output['mfe'] = float(self._lparse(line3, ' \((.*)\)$'))
        # Optional outputs
        if partition:
            # Line 4 is 'a coarse representation of the pair probabilities' and
            # the ensemble free energy
            line4 = lines.pop(0)
            output['coarse'] = self._lparse(line4, '^(.*) \[')
            output['ensemble'] = float(self._lparse(line4, ' \[(.*)\]$'))
            # Line 5 is the centroid structure, its free energy, and distance
            # to the ensemble
            line5 = lines.pop(0)
            'ensemble (.*),'
            output['frequency'] = float(self._lparse(line5, 'ensemble (.*),'))
            output['deltaG'] = float(self._lparse(line5, 'binding=(.*)$'))
            # Parse the postscript file (the only place the probability matrix
            # is)
            with open(os.path.join(self._tempdir, 'strands_dp.ps')) as f:
                pattern = 'start of base pair probability data\n(.*)\nshowpage'
                dotplot_file = f.read()
                dotplot_data = re.search(pattern, dotplot_file,
                                         flags=re.DOTALL).group(1).split('\n')
                # Dimension of the dotplot - compares seq1, seq2 to self and
                # to each other (concatenation of seq1 and seq2 = axis)
                dim = len(strand1) + len(strand2)
                ensemble_probs = np.zeros((dim, dim))
                optimal_probs = np.zeros((dim, dim))

                for point in dotplot_data:
                    point_split = point.split(' ')
                    # Use zero indexing
                    i = int(point_split[0]) - 1
                    j = int(point_split[1]) - 1
                    sqprob = float(point_split[2])
                    probtype = point_split[3]
                    if probtype == 'ubox':
                        ensemble_probs[i][j] = sqprob**2
                    else:
                        optimal_probs[i][j] = sqprob**2
                output['ensemble_matrix'] = ensemble_probs
                output['optimal_matrix'] = optimal_probs

        return output

    @tempdirs.tempdir
    def fold(self, strand, temp=37.0, dangles=2, nolp=False, nogu=False,
             noclosinggu=False, constraints=None, canonicalbponly=False,
             partition=False, pfscale=None, imfeelinglucky=False, gquad=False):
        '''Run the RNAfold command and retrieve the result in a dictionary.

        :param strand: The DNA or RNA sequence on which to run RNAfold.
        :type strand: coral.DNA or coral.RNA
        :param temp: Temperature at which to run the calculations.
        :type temp: float
        :param dangles: How to treat dangling end energies. Set to 0 to ignore
                        dangling ends. Set to 1 to limit unpaired bases to
                        at most one dangling end (default for MFE calc). Set to
                        2 (the default) to remove the limit in 1. Set to 3 to
                        allow coaxial stacking of adjacent helices in
                        .multi-loops
        :type dangles: int
        :param nolp: Produce structures without lonely pairs (isolated single
                     base pairs).
        :type nolp: bool
        :param nogu: Do not allow GU pairs.
        :type nogu: bool
        :param noclosinggu: Do not allow GU pairs at the end of helices.
        :type noclosinggu: bool
        :param constraints: Any structural constraints to use. Format is
                            defined at
                            http://www.tbi.univie.ac.at/RNA/RNAfold.1.html
        :type constraints: str
        :param canonicalbponly: Remove non-canonical base pairs from the
                                structure constraint (if applicable).
        :type canonicalbponly: bool
        :param partition: Generates the partition function, generating a coarse
                          grain structure ('coarse') in the format described at
                          http://www.itc.univie.ac.at/~ivo/RNA/RNAlib/PF-Fold.h
                          tml, the ensemble free energy ('ensemble'), the
                          centroid structure ('centroid'), the free energy of
                          the centroid structure ('centroid_fe'), and its
                          distance from the ensemble ('centroid_d').
        :type partition: int
        :param pfscale: Scaling factor for the partition function.
        :type pfScale: float
        :param imfeelinglucky: Returns the one secondary structure from the
                               Boltzmann equilibrium according to its
                               probability in the ensemble.
        :type imfeelinglucky: bool
        :param gquad: Incorporate G-Quadruplex formation into the structure
                      prediction.
        :type gquad: bool

        :returns: Dictionary of calculated values, defaulting to values of MFE
                  ('mfe': float) and dotbracket structure ('dotbracket': str).
                  More keys are added depending on keyword arguments.
        :rtype: dict

        '''
        cmd_args = []
        cmd_kwargs = {'--temp=': str(temp)}
        cmd_kwargs['--dangles='] = dangles
        if nolp:
            cmd_args.append('--noLP')
        if nogu:
            cmd_args.append('--noGU')
        if noclosinggu:
            cmd_args.append('--noClosingGU')
        if constraints is not None:
            cmd_args.append('--constraint')
            if canonicalbponly:
                cmd_args.append('--canonicalBPonly')
        if partition:
            cmd_args.append('--partfunc')
        if pfscale is not None:
            cmd_kwargs['pfScale'] = float(pfscale)
        if gquad:
            cmd_args.append('--gquad')

        inputs = [str(strand)]
        if constraints is not None:
            inputs.append(constraints)

        if strand.circular:
            cmd_args.append('--circ')
        rnafold_output = self._run('RNAfold', inputs, cmd_args, cmd_kwargs)

        # Process the output
        output = {}
        lines = rnafold_output.splitlines()
        # Line 1 is the sequence as RNA
        lines.pop(0)
        # Line 2 is the dotbracket + mfe
        line2 = lines.pop(0)
        output['dotbracket'] = self._lparse(line2, '^(.*) \(')
        output['mfe'] = float(self._lparse(line2, ' \((.*)\)$'))
        # Optional outputs
        if partition:
            # Line 3 is 'a coarse representation of the pair probabilities' and
            # the ensemble free energy
            line3 = lines.pop(0)
            output['coarse'] = self._lparse(line3, '^(.*) \[')
            output['ensemble'] = float(self._lparse(line3, ' \[(.*)\]$'))
            # Line 4 is the centroid structure, its free energy, and distance
            # to the ensemble
            line4 = lines.pop(0)
            output['centroid'] = self._lparse(line4, '^(.*) \{')
            output['centroid_fe'] = float(self._lparse(line4, '^.*{(.*) d'))
            output['centroid_d'] = float(self._lparse(line4, 'd=(.*)}$'))

        return output

    def _lparse(self, line, pattern):
        '''Parses a line from STDOUT using a regex pattern.'''
        return re.search(pattern, line).group(1)

    def _run(self, command, inputs, cmd_args, cmd_kwargs):
        arguments = [command]
        arguments += cmd_args
        for flag, value in cmd_kwargs.iteritems():
            arguments.append('{} {}'.format(flag, value))

        process = Popen(arguments, stdin=PIPE, stdout=PIPE, stderr=STDOUT,
                        cwd=self._tempdir)
        command_input = '\n'.join(inputs)
        output = process.communicate(input=command_input)[0]
        return output
