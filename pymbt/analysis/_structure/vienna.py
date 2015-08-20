# -*- coding: utf-8
'''Vienna RNA module.'''
from subprocess import Popen, PIPE, STDOUT
from tempfile import mkdtemp
from os.path import isdir
from shutil import rmtree


class Vienna(object):
    '''Run Vienna RNA functions on a sequence.'''
    def __init__(self, seqs, dotbrackets=None, constraint_structures=None):
        '''
        :param seq: DNA or RNA sequences to evaluate.
        :type seq: list of pymbt.DNA, pymbt.RNA, or str
        :param dotbrackets: Dot bracket formatted structures of given seuqneces
        :type dotbrackets: list of str

        :returns: pymbt.analysis.Vienna instance.

        '''
        self._seqs = [str(seq) for seq in seqs]
        self._tempdir = ''
        if dotbrackets is None:
            self.dotbrackets = []
        else:
            self.dotbrackets = dotbrackets
        if constraint_structures is None:
            self.constraint_structures = []
        else:
            self.constraint_structures = constraint_structures

    def free_energy_of_structure(self, temp=30.0, index=None):
        '''Calculate the free energy(ies) of sequence(s) and structure(s)
        at a given temperature.
        :param temp: Temperature at which to run calculations (°C).
        :type temp: float
        :param index: Index of the sequence for which to calculate mfe.
        :type index: int
        :returns: Minimum Free Energy (mfe) list [kcal/mol].
        :rtype: list of floats

        '''
        # FIXME: will error if not initialized with dotbrackets
        def run_rnaeval(seq, dotbracket):
            process = Popen(['RNAeval', '-T', str(temp)], stdin=PIPE,
                            stdout=PIPE, stderr=STDOUT)
            rnaeval_input = seq + '\n' + dotbracket
            output = process.communicate(input=rnaeval_input)[0]
            lines = output.splitlines()
            lines = lines[-1].split('(')[-1].split(')')[0].strip()

            return float(lines)

        mfes = []
        if index is None:
            for seq, dotbracket in zip(self._seqs, self.dotbrackets):
                mfes.append(run_rnaeval(seq, dotbracket))
        else:
            for seq, dotbracket in zip(self._seqs[index],
                                       self.dotbrackets[index]):
                mfes.append(run_rnaeval(seq, dotbracket))
        return mfes

    def mfe(self, temp=50.0, return_structure=False, index=None,
            use_constraint_structures=False):
        '''Calculate the minimum free energy.
        :param temp: Temperature at which to run calculations.
        :type temp: float
        :param index: Specific sequence index on which to calcualte mfe.
        :type index: int
        :returns: Minimum Free Energy (mfe) list [kcal/mol].
        :rtype: list of floats

        '''
        mfes = []
        dotbrackets = []

        def run_rnafold(seq, constraint=''):
            arguments = ['RNAfold', '-T', str(temp), '-noPS', '-C']

            process = Popen(arguments, stdin=PIPE,
                            stdout=PIPE, stderr=STDOUT)
            rnafold_input = '\n'.join([seq, constraint])
            output = process.communicate(input=rnafold_input)[0]
            lines = output.splitlines()
            mfe = lines[-1].split('(')[-1].split(')')[0].strip()
            dotbracket = lines[-1].split(' ')[0].strip()

            return mfe, dotbracket

        if use_constraint_structures:
            if index is None:
                for seq, const in zip(self._seqs, self.constraint_structures):
                    mfe, dotbracket = run_rnafold(seq, constraint=const)
                    mfes.append(float(mfe))
                    dotbrackets.append(dotbracket)
            else:
                const = self.sequence_constraints[index]
                mfe, dotbracket = run_rnafold(self._seqs[index],
                                              constraint=const)
                mfes.append(float(mfe))
                dotbrackets.append(dotbracket)
        else:
            if index is None:
                for seq in self._seqs:
                    mfe, dotbracket = run_rnafold(seq)
                    mfes.append(float(mfe))
                    dotbrackets.append(dotbracket)
            else:
                mfe, dotbracket = run_rnafold(self._seqs[index])
                mfes.append(float(mfe))
                dotbrackets.append(dotbracket)

        if return_structure:
            return zip(mfes, dotbrackets)
        else:
            return mfes

    def unbound_probability(self, temp=50.0, index=None):
        '''Calculate per-pair probability of being unbound (secondary
        structure).
        :param temp: Temperature at which to run calculations.
        :type temp: float
        :param index: Specific sequence index on which to calcualte pairs.
        :type index: int
        :returns: Pair probability for every base in the sequences.
        :rtype: list of list of floats.

        '''
        def run_rnafold(seq):
            process = Popen(['RNAfold', '-p', '-T', str(temp)], stdin=PIPE,
                            stdout=PIPE, stderr=STDOUT, cwd=self._tempdir)
            process.communicate(input=seq)[0]
            with open('{}/dot.ps'.format(self._tempdir), 'r') as dot:
                text = dot.read()
                text = text[text.index('%data'):]
            split = text.split('\n')
            data = [x for x in split if x.endswith('ubox')]
            data = [x.rstrip(' ubox') for x in data]
            data = [x.split() for x in data]
            data = [(int(a), int(b), float(c)) for a, b, c in data]
            unbound = [1.0] * len(seq)
            for base1, base2, prob_sqr in data:
                probability = prob_sqr**2
                unbound[base1 - 1] -= probability
                unbound[base2 - 1] -= probability
            return unbound

        self._check_tempdir()
        unbounds = []
        if index is None:
            for seq in self._seqs:
                unbound = run_rnafold(seq)
                unbounds.append(unbound)
        else:
            for seq in self._seqs[index]:
                unbound = run_rnafold(seq)
                unbounds.append(unbound)

        self._close()
        return unbounds

    def hybridization_mfe(self, temp=37.0, simtype='duplex', indices=None):
        '''Calculate the mfe of two RNA seqeunces complexing
        :param temp: Temperature at which to run calculations.
        :type temp: float
        :param indices: list (len 2) of sequence indices on which to
                       calculate pairs.
        :type indices: list of ints
        :param simtype: Vienna simulation package to use. "duplex" dose not
                        allow self RNA interactions. "cofold" allows self RNA
                        interactions to occur.
        :type simtype: str
        :returns: Binding energy between two given sequence
        :rtype: (MFE, [dotbracket1, dotbracket2], [[start1,stop1],
                [[start2,stop2]])
        :rtype: list of list of floats.

        '''
        # TODO: refactor
        if indices is None:
            tempseqs = self._seqs
        else:
            tempseqs = [self._seqs[indices[0]], self._seqs[indices[1]]]
        if simtype == 'cofold':
            process = Popen(['RNAcofold', '-T', str(temp)], stdin=PIPE,
                            stdout=PIPE, stderr=STDOUT, cwd=self._tempdir)
            stringtoinput = tempseqs[0] + '&' + tempseqs[1]
            output = process.communicate(input=stringtoinput)[0]
            lines = output.splitlines()
            mfe = float(lines[-1].split('(')[-1].split(')')[0].strip())
            # TODO: make these produce calculated values
            brackets = []
            bindlocations = []
        elif simtype == 'duplex':
            process = Popen(['RNAduplex', '-T', str(temp)], stdin=PIPE,
                            stdout=PIPE, stderr=STDOUT, cwd=self._tempdir)
            stringtoinput = tempseqs[0] + '\n' + tempseqs[1]
            output = process.communicate(input=stringtoinput)[0]
            mfe = float(output.split('(')[-1].split(')')[0])
            brackets = [output.split('&')[0],
                        output.split('&')[1].split(' ')[0]]
            bind_location_1 = [int(i) for i in
                               output.split('  ')[1].strip().split(',')]
            bind_location_2 = [int(i) for i in
                               output.split('  ')[3].strip().split(',')]
            bindlocations = [bind_location_1, bind_location_2]
        return [(mfe, brackets, bindlocations)]

    def _close(self):
        '''Close the temporary dir (keeps /tmp clean).'''
        rmtree(self._tempdir)

    def _check_tempdir(self):
        '''If temp dir has been removed, create a new one.'''
        if not isdir(self._tempdir):
            self._tempdir = mkdtemp()
