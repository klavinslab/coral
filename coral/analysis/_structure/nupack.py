# -*- coding: utf-8
'''Wrapper for NUPACK 3.0.'''
import multiprocessing
import numpy as np
import os
import re
import subprocess
import time
from coral.utils import tempdirs


class LambdaError(Exception):
    '''Raise if maximum states is exceeded (for \'distributions\' command).'''


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
            except KeyError:
                pfunc_paths = []
                for path in os.environ['PATH'].split(os.pathsep):
                    test = os.path.join(path, 'pfunc')
                    if os.path.isfile(test) and os.access(test, os.X_OK):
                        pfunc_paths.append(test)
                if not pfunc_paths:
                    raise IOError('NUPACK commands not found - see '
                                  'documentation')

        # Initialize empty temp dir location
        self._tempdir = ''

    @tempdirs.tempdir
    def pfunc(self, strand, temp=37.0, pseudo=False, material=None,
              dangles='some', sodium=1.0, magnesium=0.0):
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: A 2-tuple of the free energy of the ordered complex
                  (float) and the partition function (float).
        :rtype: tuple

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strand, material)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=False)

        # Set up the input file and run the command
        stdout = self._run('pfunc', cmd_args, [str(strand)]).split('\n')

        return (float(stdout[-3]), float(stdout[-2]))

    @tempdirs.tempdir
    def pfunc_multi(self, strands, permutation=None, temp=37.0, pseudo=False,
                    material=None, dangles='some', sodium=1.0, magnesium=0.0):
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: A 2-tuple of the free energy of the ordered complex
                  (float) and the partition function (float).
        :rtype: tuple

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strands, material, multi=True)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=True)
        # Set up the input file and run the command
        if permutation is None:
            permutation = range(1, len(strands) + 1)
        lines = self._multi_lines(strands, permutation)
        stdout = self._run('pfunc', cmd_args, lines).split('\n')

        return (float(stdout[-3]), float(stdout[-2]))

    @tempdirs.tempdir
    def pairs(self, strand, cutoff=0.001, temp=37.0, pseudo=False,
              material=None, dangles='some', sodium=1.0, magnesium=0.0):
        '''Compute the pair probabilities for an ordered complex of strands.
        Runs the \'pairs\' command.

        :param strand: Strand on which to run pairs. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: list
        :param cutoff: Only probabilities above this cutoff appear in the
                       output.
        :type cutoff: float
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: The probability matrix, where the (i, j)th entry
                  is the probability that base i is bound to base j. The matrix
                  is augmented (it's N+1 by N+1, where N is the number of bases
                  in the sequence) with an (N+1)th column containing the
                  probability that each base is unpaired.
        :rtype: numpy.array

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strand, material)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=False)

        # Set up the input file and run the command. Note: no STDOUT
        lines = [str(strand)]
        self._run('pairs', cmd_args, lines)

        # Read the output from file
        ppairs = self._read_tempfile('pairs.ppairs')
        data = re.search('\n\n\d*\n(.*)', ppairs, flags=re.DOTALL).group(1)
        N = len(strand)
        data_lines = [line.split('\t') for line in data.split('\n') if line]
        prob_matrix = self._pairs_to_np(data_lines, N)

        return prob_matrix

    @tempdirs.tempdir
    def pairs_multi(self, strands, cutoff=0.001, permutation=None, temp=37.0,
                    pseudo=False, material=None, dangles='some', sodium=1.0,
                    magnesium=0.0):
        '''Compute the pair probabilities for an ordered complex of strands.
        Runs the \'pairs\' command.

        :param strands: List of strands to use as inputs to pairs -multi.
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :param cutoff: Only probabilities above this cutoff appear in the
                       output.
        :type cutoff: float
        :returns: Two probability matrices: The probability matrix as in the
                  pairs method (but with a dimension equal to the sum of the
                  lengths of the sequences in the permutation), and a similar
                  probability matrix where multiple strands of the same species
                  are considered to be indistinguishable.
        :rtype: list

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strands, material, multi=True)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=True)

        # Set up the input file and run the command. Note: no STDOUT
        if permutation is None:
            permutation = range(1, len(strands) + 1)
        lines = self._multi_lines(strands, permutation)
        self._run('pairs', cmd_args, lines)

        # Read the output from file
        N = sum([len(s) for s in strands])
        matrices = []
        for mat_type in ['ppairs', 'epairs']:
            data = self._read_tempfile('pairs.' + mat_type)
            probs = re.search('\n\n\d*\n(.*)', data, flags=re.DOTALL).group(1)
            lines = probs.split('\n')
            # Remove the last line (empty)
            lines.pop()
            pairlist = [line.split('\t') for line in lines]
            prob_matrix = self._pairs_to_np(pairlist, N)
            matrices.append(prob_matrix)

        return matrices

    @tempdirs.tempdir
    def mfe(self, strand, degenerate=False, temp=37.0, pseudo=False,
            material=None, dangles='some', sodium=1.0, magnesium=0.0):
        '''Compute the MFE for an ordered complex of strands. Runs the \'mfe\'
        command.

        :param strand: Strand on which to run mfe. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: coral.DNA or coral.RNA
        :param degenerate: Setting to True will result in returning a list of
                           dictionaries associated with structures having the
                           same, minimal MFE value.
        :type degenerate: bool
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: A dictionary with keys for 'mfe' (a float), 'dotparens'
                  (dot-parens notation of the MFE structure), and 'pairlist'
                  (a pair list notation of the MFE structure). Note that the
                  pair list will be an empty list if the MFE is unstructured.
        :rtype: dict

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strand, material)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=False)
        if degenerate:
            cmd_args.append('-degenerate')

        # Set up the input file and run the command. Note: no STDOUT
        lines = [str(strand)]
        self._run('mfe', cmd_args, lines)

        # Read the output from file
        structures = self._process_mfe(self._read_tempfile('mfe.mfe'))

        if degenerate:
            return structures
        else:
            return structures[0]

    @tempdirs.tempdir
    def mfe_multi(self, strands, permutation=None, degenerate=False, temp=37.0,
                  pseudo=False, material=None, dangles='some', sodium=1.0,
                  magnesium=0.0):
        '''Compute the MFE for an ordered complex of strands. Runs the \'mfe\'
        command.

        :param strands: Strands on which to run mfe. Strands must be either
                       coral.DNA or coral.RNA).
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :param degenerate: Setting to True will result in returning a list of
                           dictionaries associated with structures having the
                           same, minimal MFE value.
        :type degenerate: bool
        :returns: A dictionary with keys for 'mfe' (a float), 'dotparens'
                  (dot-parens notation of the MFE structure), and 'pairlist'
                  (a pair list notation of the MFE structure). Note that the
                  pair list will be an empty list if the MFE is unstructured.
        :rtype: dict

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strands, material, multi=True)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=True)
        if degenerate:
            cmd_args.append('-degenerate')

        # Set up the input file and run the command. Note: no STDOUT
        if permutation is None:
            permutation = range(1, len(strands) + 1)
        lines = self._multi_lines(strands, permutation)
        self._run('mfe', cmd_args, lines)

        # Read the output from file
        structures = self._process_mfe(self._read_tempfile('mfe.mfe'))

        if degenerate:
            return structures
        else:
            return structures[0]

    @tempdirs.tempdir
    def subopt(self, strand, gap, temp=37.0, pseudo=False, material=None,
               dangles='some', sodium=1.0, magnesium=0.0):
        '''Compute the suboptimal structures within a defined energy gap of the
        MFE. Runs the \'subopt\' command.

        :param strand: Strand on which to run subopt. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: coral.DNA or coral.RNA
        :param gap: Energy gap within to restrict results, e.g. 0.1.
        :type gap: float
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: A list of dictionaries of the type returned by .mfe().
        :rtype: list

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strand, material)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=False)

        # Set up the input file and run the command. Note: no STDOUT
        lines = [str(strand), str(gap)]
        self._run('subopt', cmd_args, lines)

        # Read the output from file
        structures = self._process_mfe(self._read_tempfile('subopt.subopt'))

        return structures

    @tempdirs.tempdir
    def subopt_multi(self, strands, gap, permutation=None, temp=37.0,
                     pseudo=False, material=None, dangles='some', sodium=1.0,
                     magnesium=0.0):
        '''Compute the suboptimal structures within a defined energy gap of the
        MFE for an ordered permutation of strands. Runs the \'subopt\' command.

        :param strands: Strands on which to run subopt. Strands must be either
                       coral.DNA or coral.RNA).
        :type strands: list
        :param gap: Energy gap within to restrict results, e.g. 0.1.
        :type gap: float
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: A list of dictionaries of the type returned by .mfe().
        :rtype: list

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strands, material, multi=True)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=True)

        # Set up the input file and run the command. Note: no STDOUT
        if permutation is None:
            permutation = range(1, len(strands) + 1)
        lines = self._multi_lines(strands, permutation)
        lines.append(str(gap))
        self._run('subopt', cmd_args, lines)

        # Read the output from file
        structures = self._process_mfe(self._read_tempfile('subopt.subopt'))

        return structures

    @tempdirs.tempdir
    def count(self, strand, pseudo=False):
        '''Enumerates the total number of secondary structures over the
        structural ensemble Ω(π). Runs the \'count\' command.

        :param strand: Strand on which to run count. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: list
        :param pseudo: Enable pseudoknots.
        :type pseudo: bool
        :returns: The count of the number of structures in the structural
                  ensemble.
        :rtype: int

        '''
        # Set up command flags
        if pseudo:
            cmd_args = ['-pseudo']
        else:
            cmd_args = []

        # Set up the input file and run the command
        stdout = self._run('count', cmd_args, [str(strand)]).split('\n')

        # Return the count
        return int(float(stdout[-2]))

    @tempdirs.tempdir
    def count_multi(self, strands, permutation=None, pseudo=False):
        '''Enumerates the total number of secondary structures over the
        structural ensemble Ω(π) with an ordered permutation of strands. Runs
        the \'count\' command.

        :param strands: List of strands to use as inputs to count -multi.
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
        :returns: Dictionary with the following key:value pairs: 'energy':
                  free energy, 'pfunc': partition function.
        :rtype: dict

        '''
        # Set up command flags
        cmd_args = ['-multi']
        if pseudo:
            cmd_args.append('-pseudo')

        # Set up the input file and run the command
        if permutation is None:
            permutation = range(1, len(strands) + 1)
        lines = self._multi_lines(strands, permutation)
        stdout = self._run('count', cmd_args, lines).split('\n')

        return int(float(stdout[-2]))

    @tempdirs.tempdir
    def energy(self, strand, dotparens, temp=37.0, pseudo=False, material=None,
               dangles='some', sodium=1.0, magnesium=0.0):
        '''Calculate the free energy of a given sequence structure. Runs the
        \'energy\' command.

        :param strand: Strand on which to run energy. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: coral.DNA or coral.RNA
        :param dotparens: The structure in dotparens notation.
        :type dotparens: str
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: The free energy of the sequence with the specified secondary
                  structure.
        :rtype: float

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strand, material)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=False)

        # Set up the input file and run the command. Note: no STDOUT
        lines = [str(strand), dotparens]
        stdout = self._run('energy', cmd_args, lines).split('\n')

        # Return the energy
        return float(stdout[-2])

    @tempdirs.tempdir
    def energy_multi(self, strands, dotparens, permutation=None, temp=37.0,
                     pseudo=False, material=None, dangles='some', sodium=1.0,
                     magnesium=0.0):
        '''Calculate the free energy of a given sequence structure. Runs the
        \'energy\' command.

        :param strands: Strands on which to run energy. Strands must be either
                       coral.DNA or coral.RNA).
        :type strands: list
        :param dotparens: The structure in dotparens notation.
        :type dotparens: str
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: The free energy of the sequence with the specified secondary
                  structure.
        :rtype: float

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strands, material, multi=True)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=True)
        # Set up the input file and run the command
        if permutation is None:
            permutation = range(1, len(strands) + 1)
        lines = self._multi_lines(strands, permutation)
        lines.append(dotparens)
        stdout = self._run('energy', cmd_args, lines).split('\n')

        return float(stdout[-2])

    @tempdirs.tempdir
    def prob(self, strand, dotparens, temp=37.0, pseudo=False, material=None,
             dangles='some', sodium=1.0, magnesium=0.0):
        '''Calculate the equilibrium probability of a given secondary
        structure. Runs the \'prob\' command.

        :param strand: Strand on which to run prob. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: coral.DNA or coral.RNA
        :param dotparens: The structure in dotparens notation.
        :type dotparens: str
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: The equilibrium probability of the sequence with the
                  specified secondary structure.
        :rtype: float

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strand, material)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=False)

        # Set up the input file and run the command.
        lines = [str(strand), dotparens]
        stdout = self._run('prob', cmd_args, lines).split('\n')

        # Return the probabilities
        return float(stdout[-2])

    @tempdirs.tempdir
    def prob_multi(self, strands, dotparens, permutation=None, temp=37.0,
                   pseudo=False, material=None, dangles='some', sodium=1.0,
                   magnesium=0.0):
        '''Calculate the equilibrium probability of a given secondary
        structure of a complex of sequences in a given circular permutation.
        Runs the \'prob\' command.

        :param strands: Strands on which to run prob. Strands must be either
                       coral.DNA or coral.RNA).
        :type strands: list
        :param dotparens: The structure in dotparens notation.
        :type dotparens: str
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: The free energy of the sequence with the specified secondary
                  structure.
        :rtype: float

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strands, material, multi=True)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=True)
        # Set up the input file and run the command
        if permutation is None:
            permutation = range(1, len(strands) + 1)
        lines = self._multi_lines(strands, permutation)
        lines.append(dotparens)
        stdout = self._run('prob', cmd_args, lines).split('\n')

        return float(stdout[-2])

    @tempdirs.tempdir
    def defect(self, strand, dotparens, mfe=False, temp=37.0, pseudo=False,
               material=None, dangles='some', sodium=1.0, magnesium=0.0):
        '''Calculate the ensemble defect for a given sequence and secondary
        structure. From the documentation, the ensemble defect is defined as
        \'the number of incorrectly paired nucleotides at equilibrium evaluated
        over the ensemble of the ordered complex.\' Runs the \'defect\'
        command.

        :param strand: Strand on which to run defect. Strands must be either
                       coral.DNA or coral.RNA).
        :type strand: coral.DNA or coral.RNA
        :param dotparens: The structure in dotparens notation.
        :type dotparens: str
        :param mfe: Return the MFE defect (Zadeh et al., 2010) instead of the
                    ensemble defect.
        :type mfe: bool
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: A 2-tuple of the ensemble defect (float) and normalized
                  ensemble defect (float).
        :rtype: tuple

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strand, material)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=False)
        if mfe:
            cmd_args.append('-mfe')

        # Set up the input file and run the command.
        lines = [str(strand), dotparens]
        stdout = self._run('defect', cmd_args, lines).split('\n')

        # Return the defect [ensemble defect, ensemble defect]
        return (float(stdout[-3]), float(stdout[-2]))

    @tempdirs.tempdir
    def defect_multi(self, strands, dotparens, permutation=None, mfe=False,
                     temp=37.0, pseudo=False, material=None, dangles='some',
                     sodium=1.0, magnesium=0.0):
        '''Calculate the free energy of a given sequence structure. Runs the
        \'energy\' command.

        :param strands: Strands on which to run energy. Strands must be either
                       coral.DNA or coral.RNA).
        :type strands: list
        :param dotparens: The structure in dotparens notation.
        :type dotparens: str
        :param permutation: The circular permutation of strands to test in
                            complex. e.g. to test in the order that was input
                            for 4 strands, the permutation would be [1,2,3,4].
                            If set to None, defaults to the order of the
                            input strands.
        :type permutation: list
        :param mfe: Return the MFE defect (Zadeh et al., 2010) instead of the
                    ensemble defect.
        :type mfe: bool
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: A 2-tuple of the ensemble defect (float) and normalized
                  ensemble defect (float).
        :rtype: tuple

        '''
        # Set the material (will be used to set command material flag)
        material = self._set_material(strands, material, multi=True)

        # Set up command flags
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=True)
        if mfe:
            cmd_args.append('-mfe')

        # Set up the input file and run the command
        if permutation is None:
            permutation = range(1, len(strands) + 1)
        lines = self._multi_lines(strands, permutation)
        lines.append(dotparens)
        stdout = self._run('defect', cmd_args, lines).split('\n')

        # Return the defect [ensemble defect, ensemble defect]
        return (float(stdout[-3]), float(stdout[-2]))

    @tempdirs.tempdir
    def complexes(self, strands, max_size, ordered=False, pairs=False,
                  mfe=False, cutoff=0.001, degenerate=False, temp=37.0,
                  pseudo=False, material=None, dangles='some', sodium=1.0,
                  magnesium=0.0):
        '''
        :param strands: Strands on which to run energy. Strands must be either
                       coral.DNA or coral.RNA).
        :type strands: list of coral.DNA or coral.RNA
        :param max_size: Maximum complex size to consider (maximum number of
                         strand species in complex).
        :type max_size: int
        :param ordered: Consider distinct ordered complexes - all distinct
                        circular permutations of each complex.
        :type ordered: bool
        :param pairs: Calculate base-pairing observables as with .pairs().
        :type pairs: bool
        :param cutoff: A setting when pairs is set to True - only probabilities
                       above this threshold will be returned.
        :type cutoff: float
        :param degenerate: Applies only when \'mfe\' is set to True. If
                           set to True, the 'mfe' value associated with each
                           complex will be a list of degenerate MFEs (as in
                           the case of .mfe()).
        :type degenerate: bool
        :param temp: Temperature.
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
        :param sodium: Sodium concentration in solution (molar), only applies
                       to DNA.
        :type sodium: float
        :param magnesium: Magnesium concentration in solution (molar), only
                          applies to DNA>
        :type magnesium: float
        :returns: A list of dictionaries containing at least 'energy',
                  'complex', and 'strands' keys. If 'ordered' is True, the
                  different possible ordered permutations of complexes are
                  considered. In addition, with 'ordered' set to True, an
                  additional 'order' key describing the exact order of strands
                  and a 'permutation' index (integer) are added. If 'pairs' is
                  True, there is an additional 'epairs' key containing the
                  base-pairing expectation values. If 'mfe' is selected, 'mfe,
                  'dotparens', and 'pairlist' keys in the same as .mfe(). In
                  addition, 'mfe' sets the -ordered flag, so the same keys as
                  when 'ordered' is set to True are added.
        :rtype: float

        '''
        # TODO: Consider returning a pandas dataframe in this (and other)
        # situations to make sorting/selection between results easier.
        material = self._set_material(strands, material, multi=True)
        cmd_args = self._prep_cmd_args(temp, dangles, material, pseudo, sodium,
                                       magnesium, multi=False)
        cmd_args.append('-quiet')
        if mfe:
            cmd_args.append('-mfe')
            ordered = True
            if degenerate:
                cmd_args.append('-degenerate')
        if ordered:
            cmd_args.append('-ordered')
        if pairs:
            cmd_args.append('-pairs')
            cmd_args.append('-cutoff')
            cmd_args.append(cutoff)

        dim = sum([len(s) for s in strands])
        nstrands = len(strands)
        # Set up the input file and run the command
        lines = self._multi_lines(strands, [max_size])
        self._run('complexes', cmd_args, lines)

        # Read the output from file(s)
        if ordered:
            ocx_lines = self._read_tempfile('complexes.ocx').split('\n')
            # Process each lines
            output = []
            for line in ocx_lines:
                if line and not line.startswith('%'):
                    data = line.split('\t')
                    energy = float(data[-1])
                    complexes = [int(d) for d in data[2:2 + nstrands]]
                    permutation = int(data[1])
                    output.append({'energy': energy, 'complex': complexes,
                                   'permutation': permutation})

            key_lines = self._read_tempfile('complexes.ocx-key').split('\n')

            data_lines = [l for l in key_lines if not l.startswith('%')]
            data_lines.pop()
            for i, line in enumerate(data_lines):
                data = line.split('\t')
                keys = [int(d) for d in data[2:-1]]
                output[i]['order'] = keys

            if pairs:
                epairs_data = self._read_tempfile('complexes.ocx-epairs')
                pairslist = self._process_epairs(epairs_data)
                for i, pairs in enumerate(pairslist):
                    output[i]['epairs'] = self._pairs_to_np(pairs, dim)
                # TODO: add ocx-ppairs as well

            if mfe:
                mfe_data = self._read_tempfile('complexes.ocx-mfe')
                if degenerate:
                    return NotImplementedError('Not implemented for complexes')
                else:
                    mfe_output = self._process_mfe(mfe_data, complexes=True)
                    for i, mfedat in enumerate(mfe_output):
                        output[i]['mfe'] = mfedat['mfe']
                        output[i]['dotparens'] = mfedat['dotparens']
                        output[i]['pairlist'] = mfedat['pairlist']
        else:
            cx_lines = self._read_tempfile('complexes.cx').split('\n')
            # Remove empty last line
            cx_lines.pop()
            output = []
            for line in cx_lines:
                if not line.startswith('%'):
                    data = line.split('\t')
                    energy = float(data[-1])
                    complexes = [int(d) for d in data[1:1 + len(strands)]]
                    output.append({'energy': energy, 'complex': complexes})

            if pairs:
                # Process epairs
                epairs_data = self._read_tempfile('complexes.cx-epairs')
                pairslist = self._process_epairs(epairs_data)
                for i, pairs in enumerate(pairslist):
                    proba_mat = self._pairs_to_np(pairs, dim)
                    output[i]['epairs'] = proba_mat

        # Add strands (for downstream concentrations)
        for cx in output:
            cx['strands'] = [s.copy() for s in strands]

        return output

    @tempdirs.tempdir
    def complexes_timeonly(self, strands, max_size):
        '''Estimate the amount of time it will take to calculate all the
        partition functions for each circular permutation - estimate the time
        the actual \'complexes\' command will take to run.

        :param strands: Strands on which to run energy. Strands must be either
                       coral.DNA or coral.RNA).
        :type strands: list of coral.DNA or coral.RNA
        :param max_size: Maximum complex size to consider (maximum number of
                         strand species in complex).
        :type max_size: int
        :returns: The estimated time to run complexes' partition functions, in
                  seconds.
        :rtype: float

        '''
        cmd_args = ['-quiet', '-timeonly']
        lines = self._multi_lines(strands, [max_size])
        stdout = self._run('complexes', cmd_args, lines)
        return float(re.search('calculation\: (.*) seconds', stdout).group(1))

    @tempdirs.tempdir
    def concentrations(self, complexes, concs, ordered=False, pairs=False,
                       cutoff=0.001, temp=37.0):
        '''
        :param complexes: A list of the type returned by the complexes()
                          method.
        :type complexes: list
        :param concs: The concentration(s) of each strand species in the
                      initial complex. If they are all the same, a single
                      float can be used here.
        :type concs: list of floats or float
        :param ordered: Consider distinct ordered complexes - all distinct
                        circular permutations of each complex.
        :type ordered: bool
        :param pairs: Calculate base-pairing observables as with .pairs().
        :type pairs: bool
        :param cutoff: A setting when pairs is set to True - only probabilities
                       above this threshold will be returned.
        :type cutoff: float
        :param temp: Temperature in C.
        :type temp: float
        :returns: A list of dictionaries containing (at least) a
                  'concentrations' key. If 'pairs' is True, an 'fpairs' key
                  is added.
        :rtype: list

        '''
        # Check inputs
        nstrands = len(complexes[0]['strands'])
        try:
            if len(concs) != nstrands:
                raise ValueError('concs argument not same length as strands.')
        except TypeError:
            concs = [concs for i in range(len(complexes['strands']))]

        # Set up command-line arguments
        cmd_args = ['-quiet']
        if ordered:
            cmd_args.append('-ordered')

        # Write .con file
        with open(os.path.join(self._tempdir, 'concentrations.con')) as f:
            f.writelines(concs)

        # Write .cx or .ocx file
        header = ['%t Number of strands: {}'.format(nstrands),
                  '%\tid\tsequence']
        for i, strand in enumerate(complexes['strands']):
            header.append('%\t{}\t{}'.format(i + 1, strand))
        header.append('%\tT = {}'.format(temp))
        body = []
        for i, cx in enumerate(complexes):
            permutation = '\t'.join(complexes['complex'])
            line = '{}\t{}\t{}'.format(i + 1, permutation, complexes['energy'])
            body.append(line)

        if ordered:
            cxfile = os.path.join(self._tempdir, 'concentrations.ocx')
        else:
            cxfile = os.path.join(self._tempdir, 'concentrations.cx')

        with open(cxfile) as f:
            f.writelines(header + body)

        # Run 'concentrations'
        self._run('concentrations', cmd_args, None)

        # Parse the .eq (concentrations) file
        eq_lines = self._read_tempfile('concentrations.eq').split('\n')
        tsv_lines = [l for l in eq_lines if not l.startswith('%')]
        output = []
        for i, line in enumerate(tsv_lines):
            # It's a TSV
            data = line.split('\t')
            # Column 0 is an index
            # Columns 1-nstrands is the complex
            cx = [int(c) for c in data[1:nstrands]]
            # Column nstrands + 1 is the complex energy
            # Column nstrands + 2 is the equilibrium concentration
            eq = float(data[nstrands + 2])
            output[i] = {'complex': cx, 'concentration': eq}

        if pairs:
            # Read the .fpairs file
            pairs = self._read_tempfile('concentrations.fpairs')
            pairs_tsv = [l for l in pairs.split('\n') if not l.startswith('%')]
            # Remove first line (n complexes)
            dim = int(pairs_tsv.pop(0))
            pprob = [[int(p[0]), int(p[1]), float(p[2])] for p in pairs_tsv]
            # Convert to augmented numpy matrix
            fpairs_mat = self.pairs_to_np(pprob, dim)
            for i, out in enumerate(output):
                output[i]['fpairs'] = fpairs_mat

        return output

    @tempdirs.tempdir
    def distributions(self, complexes, counts, volume, maxstates=1e7,
                      ordered=False, temp=37.0):
        '''Runs the \'distributions\' NUPACK command. Note: this is intended
        for a relatively small number of species (on the order of ~20
        total strands for complex size ~14).

        :param complexes: A list of the type returned by the complexes()
                          method.
        :type complexes: list
        :param counts: A list of the exact number of molecules of each initial
                       species (the strands in the complexes command).
        :type counts: list of ints
        :param volume: The volume, in liters, of the container.
        :type volume: float
        :param maxstates: Maximum number of states to be enumerated, needed
                          as allowing too many states can lead to a segfault.
                          In NUPACK, this is referred to as lambda.
        :type maxstates: float
        :param ordered: Consider distinct ordered complexes - all distinct
                        circular permutations of each complex.
        :type ordered: bool
        :param temp: Temperature in C.
        :type temp: float
        :returns: A list of dictionaries containing (at least) a 'complexes'
                  key for the unique complex, an 'ev' key for the expected
                  value of the complex population and a 'probcols' list
                  indicating the probability that a given complex has
                  population 0, 1, ... max(pop) at equilibrium.
        :rtype: list
        :raises: LambdaError if maxstates is exceeded.

        '''
        # Check inputs
        nstrands = len(complexes[0]['strands'])
        if len(counts) != nstrands:
            raise ValueError('counts argument not same length as strands.')

        # Set up command-line arguments
        cmd_args = []
        if ordered:
            cmd_args.append('-ordered')

        # Write .count file
        countpath = os.path.join(self._tempdir, 'distributions.count')
        with open(countpath, 'w') as f:
            f.writelines([str(c) for c in counts] + [str(volume)])

        # Write .cx or .ocx file
        header = ['%t Number of strands: {}'.format(nstrands),
                  '%\tid\tsequence']
        for i, strand in enumerate(complexes['strands']):
            header.append('%\t{}\t{}'.format(i + 1, strand))
        header.append('%\tT = {}'.format(temp))
        body = []
        for i, cx in enumerate(complexes):
            permutation = '\t'.join(complexes['complex'])
            line = '{}\t{}\t{}'.format(i + 1, permutation, complexes['energy'])
            body.append(line)

        if ordered:
            cxfile = os.path.join(self._tempdir, 'distributions.ocx')
        else:
            cxfile = os.path.join(self._tempdir, 'distributions.cx')

        with open(cxfile) as f:
            f.writelines(header + body)

        # Run 'distributions'
        stdout = self._run('distributions', cmd_args, None)

        # Parse STDOUT
        stdout_lines = stdout.split('\n')
        if stdout_lines[0].startswith('Exceeded maximum number'):
            raise LambdaError('Exceeded maxstates combinations.')

        # pop_search = re.search('There are (*) pop', stdout_lines[0]).group(1)
        # populations = int(pop_search)
        # kT_search = re.search('of the box: (*) kT', stdout_lines[1]).group(1)
        # kT = float(kT_search)

        # Parse .dist file (comments header + TSV)
        dist_lines = self._read_tempfile('distributions.dist').split('\n')
        tsv_lines = [l for l in dist_lines if not l.startswith('%')]
        tsv_lines.pop()
        output = []
        for i, line in enumerate(tsv_lines):
            data = line.split('\t')
            # Column 0 is an index
            # Columns 1-nstrands are complexes
            cx = [int(d) for d in data[1:nstrands]]
            # Column nstrands + 1 is expected value of complex
            ev = float(data[nstrands + 1])
            # Columns nstrands + 2 and on are probability columns
            probcols = [float(d) for d in data[nstrands + 2:]]
            output[i]['complex'] = cx
            output[i]['ev'] = ev
            output[i]['probcols'] = probcols

        return output

    # Helper methods for preparing command input files
    def _multi_lines(self, strands, permutation):
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
    def _read_tempfile(self, filename):
        '''Read in and return file that's in the tempdir.

        :param filename: Name of the file to read.
        :type filename: str

        '''
        with open(os.path.join(self._tempdir, filename)) as f:
            return f.read()

    def _pairs_to_np(self, pairlist, dim):
        '''Given a set of pair probability lines, construct a numpy array.

        :param pairlist: a list of pair probability triples
        :type pairlist: list
        :returns: An upper triangular matrix of pair probabilities augmented
                  with one extra column that represents the unpaired
                  probabilities.
        :rtype: numpy.array

        '''
        mat = np.zeros((dim, dim + 1))
        for line in pairlist:
            i = int(line[0]) - 1
            j = int(line[1]) - 1
            prob = float(line[2])
            mat[i, j] = prob
        return mat

    def _process_mfe(self, data, complexes=False):
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
            if complexes:
                lines.pop(0)
            lines.pop(0)
            # Line 2 is the MFE
            mfe = float(lines.pop(0))
            # Line 3 is the dot-bracket structure
            dotparens = lines.pop(0)
            # If there are any more lines, they are a pair list format
            # structure
            pairlist = []
            for line in lines:
                pair = line.split('\t')
                pairlist.append([int(pair[0]) - 1, int(pair[1]) - 1])
            output.append({'mfe': mfe, 'dotparens': dotparens,
                           'pairlist': pairlist})

        return output

    def _process_epairs(self, filedata):
        commentline = '\n% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %\n'
        groups = filedata.split(commentline)
        groups.pop(0)
        groups.pop()
        output = []
        for group in groups[::2]:
            lines = group.split('\n')
            lines.pop(0)
            lines.pop(0)
            output.append([line.split('\t') for line in lines])
        return output

    # Helper methods for repetitive tasks
    def _set_material(self, strand_input, material, multi=False):
        if multi:
            if len(set([s.material for s in strand_input])) > 1:
                raise ValueError('Inputs must all be coral.DNA or all '
                                 'coral.RNA')
            if material is None:
                return strand_input[0].material
            else:
                return material
        else:
            if material is None:
                return strand_input.material
            else:
                return material

    def _prep_cmd_args(self, temp, dangles, material, pseudo, sodium,
                       magnesium, multi=False):
        cmd_args = []
        cmd_args += ['-T', temp]
        cmd_args += ['-dangles', dangles]
        cmd_args += ['-material', material]
        cmd_args += ['-sodium', sodium]
        cmd_args += ['-magnesium', magnesium]
        if pseudo:
            cmd_args.append('-pseudo')
        if multi:
            cmd_args.append('-multi')

        return cmd_args

    def _run(self, command, cmd_args, lines):
        prefix = command
        path = os.path.join(self._tempdir, '{}.in'.format(prefix))
        with open(path, 'w') as f:
            f.write('\n'.join(lines))

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
