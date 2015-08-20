'''Evaluate windows of a sequence for in-context structure.'''
try:
    from matplotlib import pylab
except ImportError:
    print "Failed to import matplotlib. Plotting structures won't work."
import pymbt.analysis


class StructureWindows(object):
    '''Evaluate windows of structure and plot the results.'''
    def __init__(self, dna):
        '''
        :param dna: DNA sequence to analyze.
        :type dna: pymbt.DNA

        '''
        self.template = dna
        self.walked = []
        self.core_starts = []
        self.core_ends = []
        self.scores = []

    def windows(self, window_size=60, context_len=90, step=10):
        '''Walk through the sequence of interest in windows of window_size,
        evaluate free (unbound) pair probabilities.

        :param window_size: Window size in base pairs.
        :type window_size: int
        :param context_len: The number of bases of context to use when
                            analyzing each window.
        :type context_len: int
        :param step: The number of base pairs to move for each new window.
        :type step: int

        '''
        self.walked = _context_walk(self.template, window_size, context_len,
                                    step)
        self.core_starts, self.core_ends, self.scores = zip(*self.walked)
        return self.walked

    def plot(self):
        '''Plot the results of the run method.'''
        if self.walked:
            fig = pylab.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(self.core_starts, self.scores, 'bo-')
            pylab.xlabel('Core sequence start position (base pairs).')
            pylab.ylabel('Score - Probability of being unbound.')
            pylab.show()
        else:
            raise Exception("Run calculate() first so there's data to plot!")


def _context_walk(dna, window_size, context_len, step):
    '''Generate context-dependent 'non-boundedness' scores for a DNA sequence.

    :param dna: Sequence to score.
    :type dna: pymbt.DNA
    :param window_size: Window size in base pairs.
    :type window_size: int
    :param context_len: The number of bases of context to use when analyzing
                        each window.
    :type context_len: int
    :param step: The number of base pairs to move for each new window.
    :type step: int

    '''
    # Generate window indices
    window_start_ceiling = len(dna) - context_len - window_size
    window_starts = range(context_len - 1, window_start_ceiling, step)
    window_ends = [start + window_size for start in window_starts]

    # Generate left and right in-context subsequences
    l_starts = [step * i for i in range(len(window_starts))]
    l_seqs = [dna[start:end] for start, end in zip(l_starts, window_ends)]
    r_ends = [x + window_size + context_len for x in window_starts]
    r_seqs = [dna[start:end].reverse_complement() for start, end in
              zip(window_starts, r_ends)]

    # Combine and calculate nupack pair probabilities
    seqs = l_seqs + r_seqs
    pairs_run = pymbt.analysis.nupack_multi(seqs, 'dna', 'pairs', {'index': 0})
    # Focus on pair probabilities that matter - those in the window
    pairs = [run[-window_size:] for run in pairs_run]
    # Score by average pair probability
    lr_scores = [sum(pair) / len(pair) for pair in pairs]

    # Split into left-right contexts again and sum for each window
    l_scores = lr_scores[0:len(seqs) / 2]
    r_scores = lr_scores[len(seqs) / 2:]
    scores = [(l + r) / 2 for l, r in zip(l_scores, r_scores)]

    # Summarize and return window indices and score
    summary = zip(window_starts, window_ends, scores)

    return summary
