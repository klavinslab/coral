'''Sanger sequencing alignment tools.'''
import coral.analysis

# FIXME: sequencing that goes past 'end' of a circular reference
# is reported as an insertion
# TODO: consensus / master sequence for plotting / report / analysis


class Sanger(object):
    '''Align and analyze Sanger sequencing results.'''
    def __init__(self, reference, results, material='DNA'):
        '''
        :param reference: Reference sequence.
        :type reference: :class:`coral.DNA`
        :param results: Sequencing result string. A list of DNA objects is also
                        valid.
        :type results: list of coral.DNA sequences
        :param material: Material of the alignment sequences - 'DNA' or
                         'Peptide'.
        :type material: str
        :returns: instance of coral.analysis.Sanger (contains alignment and
                  provides analysis/visualization methods

        '''
        # Alignment params / thresholds
        self._gap_open = -25
        self._gap_extend = 0
        self._score_threshold = 100

        # Only one input? Put in list so it can be processed the same.
        if type(results) != list:
            results = [results]
        # Sequences and calculations that get reused
        self._reference = reference
        self._results_input = results
        self._processed = [self._remove_n(x) for x in results]
        self.names = [seq.name for seq in results]

        self._material = material
        if self._material == 'DNA':
            self._matrix = 'DNA_simple'
        elif self._material == 'Peptide':
            self._matrix = 'BLOSUM62'
        else:
            raise ValueError('material arugment must be DNA or Peptide')

        # Align
        print '(Aligning...)'
        self.alignments, self.scores = self._align()
        # Sequencing coverage
        self.coverage = self._find_coverage()
        self._differences = self._find_differences()
        self._difference_n = [sum([len(x) for x in difference]) for difference
                              in self._differences]
        self.mismatches, self.insertions, self.deletions = self._differences

    def report(self):
        '''Report deletions, mismatches, and insertions.'''
        difference_names = ['Mismatches', 'Insertions', 'Deletions']
        print
        print 'Summary: '
        print '--------'
        print
        for name, num in zip(difference_names, self._difference_n):
            print '  {}: {}'.format(name, num)
        print

        for name, difference in zip(difference_names, self._differences):
            difference_n = sum([len(x) for x in difference])
            if difference_n:
                print '## {}'.format(name)
                for i, result_difference in enumerate(difference):
                    if result_difference:
                        result_name = self.names[i]
                        print '  {}'.format(result_name)
                        ref_i = self.alignments[i][0]
                        res_i = self.alignments[i][1]
                        for start, end in result_difference:
                            print
                            _sequences_display(ref_i, res_i, start, end)
                            print

    def plot(self, text=True):
        '''Plot visualization of the alignment results using matplotlib.

        :param text: Control whether text is plotted on top of sequence bars.
        :type text: bool

        '''
        try:
            from matplotlib import pyplot
            from matplotlib import cm
        except ImportError:
            raise ImportError('Optional dependency matplotlib not installed.')

        # Calculations
        # Reference bar information
        reference_x = 0
        reference_y = 14.25
        reference_width = len(self.alignments[0][0])
        reference_height = 1

        # Bin the features so they don't overlap when plotted
        features = self._reference.features
        feature_ranges = [(feature.start, feature.stop) for feature in
                          features]
        feature_bins = disjoint_bins(feature_ranges)
        feature_nbin = len(feature_bins)

        # Bin the alignments so they don't overlap when plotted
        alignment_bins = disjoint_bins(self.coverage)

        # Calculate discrepancy coordinates
        discrepancy_coords = [[], [], []]
        differences = [self.mismatches, self.insertions, self.deletions]
        for i, result_bin in enumerate(alignment_bins):
            for index in result_bin:
                for j, discrepancy_type in enumerate(differences):
                    for discrepancy in discrepancy_type[index]:
                        coord_y = i
                        coord_x = (discrepancy[0] + discrepancy[1]) // 2
                        coords = (coord_x, coord_y)
                        discrepancy_coords[j].append(coords)

        # Plotting
        # Plot calculations
        # Controls spacing between bars, size of bars, etc
        size = 10

        # Matplotlib commands
        # Plot a black reference bar at the bottom
        fig = pyplot.figure(figsize=(12, 9), dpi=90)
        sub1 = fig.add_subplot(111)
        sub1.broken_barh([(reference_x, reference_width)],
                         (reference_y, reference_height),
                         facecolors='black', edgecolors='none')

        dotted_x = (0, len(self.alignments[0][0]))
        height = (feature_nbin + 1.5) * size
        dotted_y = (height, height)
        sub1.plot(dotted_x, dotted_y, linestyle='dotted')

        # Plot the reference features on top of the bar
        # Plot the features by bin:
        for i, feature_bin in enumerate(feature_bins):
            for index in feature_bin:
                feature = features[index]
                name = feature.name

                width = feature.stop - feature.start
                height = size - 1
                y_index = (i + 1) * size
                mid = (feature.start + feature.stop) // 2

                pos = float(i) / len(features)
                sub1.broken_barh([(feature.start, width)],
                                 (y_index, height),
                                 facecolors=cm.Set3(pos),
                                 edgecolors='black')
                if text:
                    sub1.text(mid, y_index + size / 2, name, rotation=90,
                              size='smaller')

        # Plot sequencing results by bin
        for i, result_bin in enumerate(alignment_bins):
            for index in result_bin:
                start, stop = self.coverage[index]
                name = self.names[index]

                width = stop - start
                height = size - 1
                y_index = (i + feature_nbin + 2) * size
                text_x = start + (stop - start) // 8

                sub1.broken_barh([(start, width)], (y_index, height),
                                 facecolors='pink', edgecolors='black')
                if text:
                    sub1.text(text_x, y_index + size // 3, wrap_name(name),
                              rotation=0, size='smaller')

        # Plot mismatches, insertions, deletions
        sub1.plot(1000, 25)
        shapes = ['o', 'x', 'v']
        labels = ['mismatch', 'insertion', 'deletion']
        for coords, shape, label in zip(discrepancy_coords, shapes, labels):
            x_coords = [x[0] for x in coords]
            y_coords = [(x[1] + feature_nbin + 2) * size + 2 for x in coords]
            sub1.scatter(x_coords, y_coords, marker=shape, color='k',
                         label=label)
        sub1.legend()

        # Plot labeling, etc
        sub1.set_xlim(0, reference_width)
        sub1.set_xlabel('Base pairs from origin')
        sub1.set_yticks([15, 15 + size * (feature_nbin + 1)])
        sub1.set_yticklabels(['Reference', 'Results'])
        sub1.grid(True)
        sub1.xaxis.grid(False)
        sub1.yaxis.grid(False)

        pyplot.title('Alignment gap summary', fontstyle='normal')
        pyplot.show()

    def _remove_n(self, seq):
        '''Find largest non-N segment.

        :param seq: Sequence that contains Ns to remove
        :type seq: coral.DNA

        '''
        seq_str = str(seq)
        largest = max([x for x in seq_str.split('N')], key=len)
        seq_start = seq_str.index(largest)
        seq_stop = seq_start + len(largest)
        processed = seq[seq_start:seq_stop]
        return processed

    def _align(self):
        '''Align sequences using needle.'''
        # Align
        refs = [self._reference.copy() for x in self._processed]
        if len(self._processed) > 1:
            needle = coral.analysis.needle_multi(refs, self._processed,
                                                 gap_open=self._gap_open,
                                                 gap_extend=self._gap_extend,
                                                 matrix=self._matrix)
        else:
            needle_single = coral.analysis.needle(refs[0], self._processed[0],
                                                  gap_open=self._gap_open,
                                                  gap_extend=self._gap_extend,
                                                  matrix=self._matrix)
            needle = [needle_single]
        # Split into alignments and scores
        # TODO: use zip here to make it simpler
        alignments = [(str(ref), str(res)) for ref, res, score in
                      needle]
        scores = [result[2] for result in needle]
        # If a result scores too low, try reverse complement
        # TODO: Find their indices and use multiprocessing
        # TODO: If score is too low, add warning but don't show in alignment
        failed = []
        # reversed so we can pop bad alignments (not show them in the plot)
        for i, score in reversed(list(enumerate(scores))):
            if score < self._score_threshold:
                swapped_result = self._processed[i].reverse_complement()
                new_needle = coral.analysis.needle(self._reference,
                                                   swapped_result,
                                                   gap_open=self._gap_open,
                                                   gap_extend=self._gap_extend,
                                                   matrix=self._matrix)
                alignments[i] = (str(new_needle[0]), str(new_needle[1]))
                score = new_needle[2]
            if score < self._score_threshold:
                failed.append(self._processed[i].name)
                # Remove the failed alignment!
                alignments.pop(i)
        if failed:
            msg = 'The following results fell below the score threshold: '
            print msg + '{}'.format(failed)
        # Trim to reference - if reference is shorter than results
        # (e.g. reference is just plasmid insert, results start earlier)
        for i, (reference, result) in enumerate(alignments):
            f_count = 0
            r_count = 0
            # Count beginning gaps (hyphens), then trim
            for char in reference:
                if char != '-':
                    break
                else:
                    f_count += 1
            for char in reversed(reference):
                if char != '-':
                    break
                else:
                    r_count += 1

            if r_count:
                trimmed_ref = reference[f_count:]
                trimmed_res = result[f_count:]
            else:
                trimmed_ref = reference[f_count:]
                trimmed_res = result[f_count:]
            alignments[i] = (trimmed_ref, trimmed_res)

        return alignments, scores

    def _find_coverage(self):
        '''Find sequencing coverage using trailing and leading gaps.'''
        left_coverage = [len(res) - len(res.lstrip('-')) for ref, res in
                         self.alignments]
        right_coverage = [len(res.rstrip('-')) for ref, res in self.alignments]
        return zip(left_coverage, right_coverage)

    def _find_differences(self):
        '''Calculate mismatches, insertions, and deletions.'''
        all_mismatches = [[] for i in range(len(self.alignments))]
        all_insertions = [[] for i in range(len(self.alignments))]
        all_deletions = [[] for i in range(len(self.alignments))]
        # Got through each alignment
        for i, (reference, result) in enumerate(self.alignments):
            # Go through each alignment base by base
            coverage = self.coverage[i]
            ref_trim = reference[slice(*coverage)]
            res_trim = result[slice(*coverage)]
            for j, (refbase, resbase) in enumerate(zip(ref_trim, res_trim)):
                # If only the reference base is a '-', is insertion
                position = j + coverage[0]
                if refbase == '-':
                    all_insertions[i].append((position, position))
                # If only the reference base is a '-', is deletion
                elif resbase == '-':
                    all_deletions[i].append((position, position))
                # If bases don't match, is mismatch
                elif refbase != resbase:
                    all_mismatches[i].append((position, position))

        # Group mismatches/indels if they appear one after another
        all_mismatches = [_group_differences(mismatches) for mismatches in
                          all_mismatches]
        all_insertions = [_group_differences(insertions) for insertions in
                          all_insertions]
        all_deletions = [_group_differences(deletions) for deletions in
                         all_deletions]
        return all_mismatches, all_insertions, all_deletions

    def __repr__(self):
        '''Representation of a Sanger sequencing object.'''
        disc_sum = sum(self._difference_n)
        head = 'An alignment with {} discrepancies.'.format(disc_sum)
        if disc_sum:
            mismatches = 'mismatches: {}, '.format(self._difference_n[0])
            insertions = 'insertions: {}, '.format(self._difference_n[1])
            deletions = 'deletions: {}'.format(self._difference_n[2])
            diffs = '\n' + mismatches + insertions + deletions
        else:
            diffs = ''
        return head + diffs


def _group_differences(difference_list):
    '''Group adjoining insertions or deletions, respectively.

    :param difference_list: list of difference indices (2-tuple).
    :type difference_list: list
    :returns: a list of the same format as the input, but with adjoining
              insertions or deletions grouped together

    '''
    # By setting to (-2, -2), skips first base so that others can
    # 'look back' to see whether the current deletion is last + 1
    previous = (-2, -2)
    grouped_differences = []
    for difference in difference_list:
        if difference[0] - previous[1] == 1:
            grouped_differences[-1] = (previous[0], difference[1])
        else:
            grouped_differences.append(difference)
        previous = grouped_differences[-1]
    return grouped_differences


def _sequences_display(seq1, seq2, start, stop, context=10):
    '''Display two sequences, highlighting the regions where they differ.

    :param seq1: First sequence to compare.
    :type seq1: str
    :param seq2: Second sequence to compare.
    :type seq2: str
    :param start: Where to start displaying the sequence.
    :type start: int
    :param stop: Where to stop displaying the sequence.
    :type stop: int
    :param context: Extra context to add on either side of the displayed
                    sequences.
    :type context: int

    '''
    # TODO: if seqs aren't the same size, should display a little differently
    # TODO: display is incomplete - doesn't show all mismatches / deletions /
    # insertions at once, just the ones displayed at the moment.
    # Should instead calculate a single 'overview' set of sequences:
    # top, bottom, and middle, where middle is where all the seqs match.
    # Then just subset this
    indent = 4
    # Figure out how much context can be included (seq might start/end earlier)
    l_context = (max(start - context - 1, 0), start - 1)
    r_context = (stop + 1, min(start + context + 1, len(seq1)))

    def gen_levels(seq1, seq2, context_tuple):
        '''Generate a text column: left side or right side of difference.

        :param seq1: same as seq1 of parent
        :type seq1: str
        :param seq2: same as seq2 of parent
        :type seq2: str
        :context_tuple: l_context or r_context
        :type context_tuple: tuple

        '''
        top, bottom = [seq[slice(*context_tuple)] for seq in [seq1, seq2]]
        context_len = (context_tuple[1] - context_tuple[0])
        middle = '|' * context_len
        highlight = ' ' * context_len
        return top, bottom, middle, highlight

    # Generate left and right columns of text to display
    l_top, l_bottom, l_middle, l_highlight = gen_levels(seq1, seq2, l_context)
    r_top, r_bottom, r_middle, r_highlight = gen_levels(seq1, seq2, r_context)

    # Generate core column - core is where differences are
    core_slice = slice(start, stop + 1)
    core_len = stop - start + 1
    core_top = seq1[core_slice]
    core_bottom = seq2[core_slice]
    core_middle = ' ' * core_len
    core_highlight = '*' * core_len

    # Combine each level together for printing
    top = l_top + core_top + r_top
    middle = l_middle + core_middle + r_middle
    bottom = l_bottom + core_bottom + r_bottom
    highlight = l_highlight + core_highlight + r_highlight

    indent = ' ' * indent
    if start != stop:
        print '{0}Positions {1} to {2}:'.format(indent, start, stop)
    else:
        print '{0}Position {1}:'.format(indent, start)

    print indent + top
    print indent + middle
    print indent + bottom
    print indent + highlight


def disjoint_bins(ranges_list):
    '''Construct disjoint bins given a list of 1-D ranges (tuples).

    :param range_tuple_list: A list of tuples containing range values.
    :type range_tuple_list: list
    :returns: a list of bins containing indices of the input list. e.g. if
              the third (index 2) range is in the 2nd bin, the number 2 is in
              the 2nd list.
    :rtype: list of lists of ints

    '''
    # Keep track of the original order for reporting later
    ranges_list = [(x[0], x[1], i) for i, x in enumerate(ranges_list)]
    ranges_list = sorted(ranges_list, key=lambda starts: starts[0])

    remaining = ranges_list[:]
    binned = []
    while True:
        current_bin = []
        next_bin = []
        while remaining:
            last = remaining.pop(0)
            current_bin.append(last)
            next_bin += [x for x in remaining if x[0] < last[1]]
            remaining = [x for x in remaining if x[0] >= last[1]]
        binned.append(current_bin)

        if not remaining and not next_bin:
            break
        else:
            remaining = next_bin

    bins = [[x[2] for x in range_bin] for range_bin in binned]

    return bins


def wrap_name(str_in, wrap_len=15):
    '''Wrap plotted text to avoid overlaps (not perfect).

    :param str_in: Input string.
    :type str_in: str
    :param wrap_len: Length at which to wrap lines.
    :type wrap_len: str
    :returns: Wrapped text
    :rtype: str

    '''
    wrap_positions = range(0, len(str_in), wrap_len)
    out = [str_in[i:i + wrap_len] for i in wrap_positions]
    return '\n'.join(out)
