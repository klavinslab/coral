'''Sanger sequencing alignment tools.'''
import coral as cr

# FIXME: sequencing that goes past 'end' of a circular reference
# is reported as an insertion
# TODO: consensus / master sequence for plotting / report / analysis


class Sanger(object):
    '''Align and analyze Sanger sequencing results.'''

    def __init__(self, reference, results, method='needle',
                 method_kwargs=None):
        '''
        :param reference: Reference sequence.
        :type reference: :class:`coral.DNA`
        :param results: Sequencing result string. A list of DNA objects is also
                        valid.
        :type results: list of coral.DNA sequences
        :param method: Alignment method to use. Options are:
                         \'needle\': Uses coral.analysis.needle_msa
                         \'MAFFT\': Uses coral.analysis.MAFFT
        :type method: str
        :param method_kwargs: Optional keyword arguments to send to the
                              alignment function.
        :type method_kwargs: dict
        :returns: instance of coral.analysis.Sanger (contains alignment and
                  provides analysis/visualization methods

        '''
        # Alignment params / thresholds
        # self._gap_open = -25
        # self._gap_extend = 0
        # self._score_threshold = 100

        # Sequences and calculations that get reused
        self.reference = reference
        self.results = results
        self.method = method
        if method_kwargs is None:
            self.method_kwargs = {}
        else:
            self.method_kwargs = method_kwargs

        self._remove_n()
        self.names = []
        for i, result in enumerate(self.results):
            if not result.name:
                self.names.append('Result {}'.format(i))
            else:
                self.names.append(result.name)

        # Align
        self.align()

    def align(self):
        if self.method == 'needle':
            self.alignment = cr.analysis.needle_msa(self.reference,
                                                    self.results,
                                                    **self.method_kwargs)
        elif self.method == 'MAFFT':
            self.alignment = cr.analysis.MAFFT([self.reference] + self.results,
                                               **self.method_kwargs)
        else:
            raise ValueError('Only \'needle\' or \'MAFFT\' methods allowed.')

        self.aligned_reference = self.alignment[0].copy()
        self.aligned_results = []
        for result in self.alignment[1:]:
            self.aligned_results.append(result.copy())

    def nonmatches(self):
        '''Report mismatches, indels, and coverage.'''
        # For every result, keep a dictionary of mismatches, insertions, and
        # deletions
        report = []
        for result in self.aligned_results:
            report.append(self._analyze_single(self.aligned_reference, result))

        return report

    def plot(self):
        '''Make a summary plot of the alignment and highlight nonmatches.'''
        import matplotlib.pyplot as plt
        import matplotlib.patches as patches

        # Constants to use throughout drawing
        n = len(self.results)
        nbases = len(self.aligned_reference)
        barheight = 0.4

        # Vary height of figure based on number of results
        figheight = 3 + 3 * (n - 1)
        fig = plt.figure(figsize=(9, figheight))
        ax1 = fig.add_subplot(111)

        # Plot bars to represent coverage area
        # Reference sequence
        ax1.add_patch(patches.Rectangle((0, 0), nbases, barheight,
                                        facecolor='black'))
        # Results
        for i, report in enumerate(self.nonmatches()):
            j = i + 1
            start, stop = report['coverage']
            patch = patches.Rectangle((start, j), stop - start, barheight,
                                      facecolor='darkgray')
            ax1.add_patch(patch)

            # Draw a vertical line for each type of result
            plt.vlines(report['mismatches'], j, j + barheight,
                       colors='b')
            plt.vlines(report['insertions'], j, j + barheight,
                       colors='r')

            # Terminal trailing deletions shouldn't be added
            deletions = []
            crange = range(*report['coverage'])
            deletions = [idx for idx in report['deletions'] if idx in crange]
            plt.vlines(deletions, j, j + barheight,
                       colors='g')

        ax1.set_xlim((0, nbases))
        ax1.set_ylim((-0.3, n + 1))
        ax1.set_yticks([i + barheight / 2 for i in range(n + 1)])
        ax1.set_yticklabels(['Reference'] + self.names)

        # Add legend
        mismatch_patch = patches.Patch(color='blue', label='Mismatch')
        insertion_patch = patches.Patch(color='red', label='Insertion')
        deletion_patch = patches.Patch(color='green', label='Deletion')
        plt.legend(handles=[mismatch_patch, insertion_patch, deletion_patch],
                   loc=1, ncol=3, mode='expand', borderaxespad=0.)

        plt.show()

    def _analyze_single(self, reference, result):
        '''Report mistmatches and indels for a single (aligned) reference and
        result.'''
        # TODO: Recalculate coverage based on reference (e.g. sequencing result
        # longer than template
        reference_str = str(reference)
        result_str = str(result)
        report = {'mismatches': [], 'insertions': [], 'deletions': []}
        for i, (ref, res) in enumerate(zip(reference_str, result_str)):
            if ref != res:
                # It's a mismatch or indel
                if ref == '-':
                    report['insertions'].append(i)
                if res == '-':
                    report['deletions'].append(i)
                else:
                    report['mismatches'].append(i)

        start = len(result_str) - len(result_str.lstrip('-'))
        stop = len(result_str) - len(result_str.rstrip('-'))
        report['coverage'] = [start, stop]

        return report

    def _remove_n(self):
        '''Remove terminal Ns from sequencing results.'''
        for i, result in enumerate(self.results):
            largest = max(str(result).split('N'), key=len)
            start = result.locate(largest)[0][0]
            stop = start + len(largest)
            if start != stop:
                self.results[i] = self.results[i][start:stop]
