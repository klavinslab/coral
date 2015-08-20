'''DNA object classes.'''
import collections
import os
import re
import shutil
import subprocess
import tempfile
import pymbt.analysis
import pymbt.reaction
import pymbt.seqio
from ._sequence import process_seq, reverse_complement
from ._sequence import NucleotideSequence


class DNA(NucleotideSequence):
    '''DNA sequence.'''
    def __init__(self, dna, bottom=None, topology='linear', stranded='ds',
                 features=None, run_checks=True, id=None, name=''):
        '''
        :param dna: Input sequence (DNA).
        :type dna: str
        :param bottom: Manual input of bottom-strand sequence. Enables both
                       mismatches and initializing ssDNA.
        :type bottom: str
        :param topology: Topology of DNA - 'linear' or 'circular'.
        :type topology: str
        :param stranded: Strandedness of DNA - 'ss' for single-stranded or
                         'ds' for double-stranded.
        :type stranded: str
        :param features: List of annotated features.
        :type features: list
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :param id: An optional (unique) id field for your DNA sequence.
        :type id: str
        :param name: Optional name field for your DNA sequence.
        :type name: str
        :returns: pymbt.DNA instance.
        :rtype: pymbt.DNA
        :raises: ValueError if an element of `features` isn't of type
                 pymbt.Feature.
                 ValueError if top and bottom strands have different lengths.
                 ValueError if top and bottom strands are not complementary.

        '''
        # Convert to uppercase, run alphabet check
        super(DNA, self).__init__(dna, 'dna', features=features,
                                  run_checks=run_checks)
        # Set topology
        self.topology = topology
        # Set strandedness
        self.stranded = stranded
        # If bottom was specified, check it + add it
        if bottom:
            self._bottom = bottom
            if run_checks:
                self._bottom = process_seq(bottom, 'dna')
                if len(self._bottom) != len(self._sequence):
                    msg = 'Top and bottom strands are difference lengths.'
                    raise ValueError(msg)
        else:
            self._bottom = ''.join(['-' for x in self._sequence])
            # NOTE: inefficient to assign blanks the rev comp, but cleaner code
            if stranded == 'ds':
                self._bottom = str(reverse_complement(self._sequence, 'dna'))
        # Set id
        self.id = id
        # Set name
        self.name = name

    def ape(self, ape_path=None):
        '''Open in ApE.'''
        cmd = 'ApE'
        if ape_path is None:
            # Check for ApE in PATH
            ape_executables = []
            for path in os.environ['PATH'].split(os.pathsep):
                exepath = os.path.join(path, cmd)
                ape_executables.append(os.access(exepath, os.X_OK))
            if not any(ape_executables):
                raise Exception('Ape not in PATH. Use ape_path kwarg.')
        else:
            cmd = ape_path
        # Check whether ApE exists in PATH
        tmp = tempfile.mkdtemp()
        if self.name is not None and self.name:
            filename = os.path.join(tmp, '{}.ape'.format(self.name))
        else:
            filename = os.path.join(tmp, 'tmp.ape')
        pymbt.seqio.write_dna(self, filename)
        process = subprocess.Popen([cmd, filename])
        # Block until window is closed
        try:
            process.wait()
            shutil.rmtree(tmp)
        except KeyboardInterrupt:
            shutil.rmtree(tmp)

    def bottom(self):
        '''Return the raw string of the Crick (bottom) strand.

        :returns: The Crick strand.
        :rtype: str

        '''
        return self._bottom

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely-editable copy of the current sequence.
        :rtype: pymbt.DNA

        '''
        # Significant performance improvements by skipping alphabet check
        features_copy = [feature.copy() for feature in self.features]
        return type(self)(self._sequence, bottom=self._bottom,
                          topology=self.topology, stranded=self.stranded,
                          features=features_copy, id=self.id, name=self.name,
                          run_checks=False)

    def circularize(self):
        '''Circularize linear DNA.

        :returns: A circularized version of the current sequence.
        :rtype: pymbt.DNA

        '''
        if self.top()[-1] == '-' and self.bottom()[0] == '-':
            raise ValueError('Cannot circularize - termini disconnected.')
        if self.bottom()[-1] == '-' and self.top()[0] == '-':
            raise ValueError('Cannot circularize - termini disconnected.')

        copy = self.copy()
        copy.topology = 'circular'
        return copy

    def extract(self, name, remove_subfeatures=False):
        return super(DNA, self).extract(name, 'N',
                                        remove_subfeatures=remove_subfeatures)

    def flip(self):
        '''Flip the DNA - swap the top and bottom strands.

        :returns: Flipped DNA (bottom strand is now top strand, etc.).
        :rtype: pymbt.DNA

        '''
        copy = self.copy()
        copy._sequence, copy._bottom = copy._bottom, copy._sequence
        return copy

    def gc(self):
        '''Find the frequency of G and C in the current sequence.'''

        return len([base for base in self if str(base) == 'C' or
                    str(base) == 'G'])

    def insert(self, sequence, index):
        inserted = super(DNA, self).insert(sequence, index)
        inserted.topology = self.topology
        return inserted

    def is_rotation(self, other):
        if len(self) != len(other):
            return False
        for i in range(len(self)):
            if self.rotate(i) == other:
                return True
        # If all else fails, check reverse complement
        rc = self.reverse_complement()
        for i in range(len(self)):
            if rc.rotate(i) == other:
                return True
        return False

    def linearize(self, index=0):
        '''Linearize circular DNA at an index.

        :param index: index at which to linearize.
        :type index: int
        :returns: A linearized version of the current sequence.
        :rtype: pymbt.DNA
        :raises: ValueError if the input is linear DNA.

        '''
        if self.topology == 'linear':
            raise ValueError('Cannot relinearize linear DNA.')
        copy = self.copy()
        copy.topology = 'linear'
        copy = copy[index:] + copy[:index]
        return copy

    def locate(self, pattern):
        '''Find sequences matching a pattern.

        :param pattern: Sequence for which to find matches.
        :type pattern: str
        :returns: A list of top and bottom strand indices of matches.
        :rtype: list of lists of indices (ints)
        :raises: ValueError if the pattern is longer than either the input
                 sequence (for linear DNA) or twice as long as the input
                 sequence (for circular DNA).

        '''
        # TODO: If linear, should use the methods in BaseSequence
        if self.topology == 'circular':
            if len(pattern) > 2 * len(self):
                raise ValueError('Pattern too long.')
        else:
            if len(pattern) > len(self):
                raise ValueError('Pattern too long.')

        pattern = str(pattern).upper()
        regex = '(?=' + pattern + ')'

        if self.topology == 'circular':
            r = len(pattern) - 1
            l = len(self) - r + 1
            top = self._sequence[l:] + self._sequence + self._sequence[:r]
            bottom = self._bottom[l:] + self._bottom + self._bottom[:r]
        else:
            top = self._sequence
            bottom = self._bottom

        top_starts = [index.start() for index in re.finditer(regex, top)]
        bottom_starts = [index.start() for index in re.finditer(regex, bottom)]

        # Adjust indices if doing circular search
        if self.topology == 'circular' and len(pattern) > 1:
            top_starts = [start - r + 1 for start in top_starts]
            bottom_starts = [start - r + 1 for start in bottom_starts]

        return [top_starts, bottom_starts]

    def mw(self):
        '''Calculate the molecular weight.

        :returns: The molecular weight of the current sequence.
        :rtype: float

        '''
        counter = collections.Counter((self._sequence + self._bottom).lower())
        mw_a = counter['a'] * 313.2
        mw_t = counter['t'] * 304.2
        mw_g = counter['g'] * 289.2
        mw_c = counter['c'] * 329.2
        return mw_a + mw_t + mw_g + mw_c

    def rotate(self, index):
        '''Orient DNA to index (only applies to circular DNA).

        :param index: DNA position at which to re-zero the DNA.
        :type index: int
        :returns: The current sequence reoriented at `index`.
        :rtype: pymbt.DNA
        :raises: ValueError if applied to linear sequence or `index` is
                 negative.

        '''
        if self.topology == 'linear' and index != 0:
            raise ValueError('Cannot rotate linear DNA')
        if index < 0:
            raise ValueError('Rotation index must be positive')
        else:
            return (self[index:] + self[0:index]).circularize()

    def rotate_by_feature(self, featurename):
        '''Reorient the DNA based on a feature it contains (circular DNA only).

        :param featurename: A uniquely-named feature.
        :type featurename: str
        :returns: The current sequence reoriented at the start index of a
                  unique feature matching `featurename`.
        :rtype: pymbt.DNA
        :raises: ValueError if there is no feature of `featurename` or
                 more than one feature matches `featurename`.

        '''
        # REFACTOR: Parts are redundant with .extract()
        matched = []
        for feature in self.features:
            if feature.name == featurename:
                matched.append(feature.copy())
        count = len(matched)
        if count == 1:
            return self.rotate(matched[0].start)
        elif count > 1:
            raise ValueError('More than one feature has that name.')
        else:
            raise ValueError('No such feature in the sequence.')

    def reverse_complement(self):
        '''Reverse complement the DNA.

        :returns: A reverse-complemented instance of the current sequence.
        :rtype: pymbt.DNA

        '''
        # TODO: put into NucleotideSequence class
        copy = self.copy()
        # Note: if sequence is double-stranded, swapping strand is basically
        # (but not entirely) the same thing - gaps affect accuracy.
        copy._sequence = reverse_complement(copy._sequence, 'dna')
        copy._bottom = reverse_complement(copy._bottom, 'dna')

        # Fix features (invert)
        for feature in copy.features:
            # Swap strand
            if feature.strand == 1:
                feature.strand = 0
            else:
                feature.strand = 1
            # Swap start and stop
            feature.start, feature.stop = (feature.stop, feature.start)
            # Adjust start/stop to feature len
            feature.start = len(copy) - feature.start
            feature.stop = len(copy) - feature.stop

        return copy

    def tm(self, parameters='cloning'):
        '''Find the melting temperature.

        :param parameters: The tm method to use (cloning, santalucia98,
                       breslauer86)
        :type parameters: str

        '''
        return pymbt.analysis.tm(self, parameters=parameters)

    def to_ss(self):
        '''Produce single stranded version of the current sequence.

        :returns: The current sequence, converted to ssDNA.
        :rtype: pymbt.DNA

        '''
        copy = self.copy()

        # Do nothing if already single-stranded
        if self.stranded == 'ss':
            return copy

        copy._bottom = '-' * len(copy)
        for top, bottom in zip(copy.top(), reversed(copy.bottom())):
            if top == bottom == '-':
                raise ValueError('Coercing to single-stranded would ' +
                                 'introduce a double stranded break.')
        copy.stranded = 'ss'

        return copy

    def to_ds(self):
        '''Produce double stranded version of the current sequence.

        :returns: The current sequence, converted to dsDNA.
        :rtype: pymbt.DNA

        '''
        # TODO: protect .stranded attribute if requiring setter method
        copy = self.copy()
        # Do nothing if already set
        if self.stranded == 'ds':
            return copy

        # Find strand that's all gaps (if ss this should be the case)
        reverse_seq = self.reverse_complement()
        if all([char == '-' for char in self._sequence]):
            copy._sequence = reverse_seq._bottom
        elif all([char == '-' for char in self._bottom]):
            copy._bottom = reverse_seq._sequence
        copy.stranded = 'ds'

        return copy

    def top(self):
        '''Return the raw string of the Watson (top) strand.

        :returns: The Watson strand.
        :rtype: str

        '''
        return self._sequence

    def transcribe(self):
        '''Transcribe into RNA.

        :returns: An RNA sequence transcribed from the current DNA sequence.
        :rtype: pymbt.RNA

        '''
        return pymbt.reaction.transcribe(self)

    def __add__(self, other):
        '''Add DNA together.

        :param other: instance to be added to.
        :type other: compatible sequence object (currently only DNA).
        :returns: Concatenated DNA sequence.
        :rtype: pymbt.DNA
        :raises: Exception if either sequence is circular.
                 Exception if concatenating a sequence with overhangs would
                 create a discontinuity.

        '''
        if type(self) != type(other):
            try:
                other = type(self)(other)
            except AttributeError:
                raise TypeError('Cannot add {} to {}'.format(self, other))

        if self.topology == 'circular' or other.topology == 'circular':
            raise Exception('Can only add linear DNA.')

        discontinuity = [False, False]
        if len(self) != 0 and len(other) != 0:
            # If either is empty, let things proceed anyways
            discontinuity[0] = (self._sequence[-1] == '-' and
                                other._bottom[-1] == '-')
            discontinuity[1] = (self._bottom[0] == '-' and
                                other._sequence[0] == '-')

        for_discontinuity = discontinuity[0]
        rev_discontinuity = discontinuity[1]

        if for_discontinuity or rev_discontinuity:
            msg = 'Concatenated DNA would be discontinuous.'
            raise Exception(msg)

        if self.stranded == 'ds' or other.stranded == 'ds':
            stranded = 'ds'
        else:
            stranded = 'ss'

        tops = self._sequence + other._sequence
        bottoms = other._bottom + self._bottom
        self_features = [feature.copy() for feature in self.features]
        other_features = [feature.copy() for feature in other.features]
        for feature in other_features:
            feature.move(len(self))
        features = self_features + other_features

        new_instance = DNA(tops, bottom=bottoms, topology='linear',
                           stranded=stranded, run_checks=False,
                           features=features)

        return new_instance

    def __contains__(self, query):
        '''Defines `query in sequence` operator.

        :param query: query string or DNA sequence
        :type query: str or pymbt.DNA

        '''
        # query in forward sequence
        if super(DNA, self).__contains__(query, 'N'):
            return True
        # query in reverse complement
        elif super(DNA, self.reverse_complement()).__contains__(query, 'N'):
            return True
        # query in neither
        else:
            return False

    def __delitem__(self, index):
        '''Delete sequence at an index.

        :param index: index to delete
        :type index: int
        :returns: The current sequence with the base at `index` removed.
        :rtype: pymbt.DNA

        '''
        super(DNA, self).__delitem__(index)
        bottom_list = list(self._bottom[::-1])
        del bottom_list[index]
        self._bottom = ''.join(bottom_list)[::-1]

    def __getitem__(self, key):
        '''Index and slice sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object
        :returns: A subsequence matching the slice (`key`).
        :rtype: pymbt.DNA

        '''
        # Use BaseSequence method to assign top strand and figure out features
        if isinstance(key, slice):
            if all([k is None for k in [key.start, key.stop, key.step]]):
                # It's the copy slice operator ([:])
                return self.copy()
            else:
                # The key is a normal slice
                copy = super(DNA, self).__getitem__(key)

                # bottom_key = slice(-key.stop if key.stop is not None
                #                    else None,
                #                    -key.start if key.start is not None
                #                    else None,
                #                    key.step)
                # copy._bottom = copy._bottom[bottom_key]
                copy._bottom = copy._bottom[::-1][key][::-1]
        else:
            # The key is an integer
            copy = super(DNA, self).__getitem__(key)
            copy._bottom = copy._bottom[-key]

        copy.topology = 'linear'

        return copy

    def __eq__(self, other):
        '''Define equality - sequences, topology, and strandedness are the
        same.

        :returns: Whether current sequence's (Watson and Crick), topology,
                  and strandedness are equivalent to those of another sequence.
        :rtype: bool

        '''
        tops_equal = self._sequence == other._sequence
        bottoms_equal = self._bottom == other._bottom
        topology_equal = self.topology == other.topology
        stranded_equal = self.stranded == other.stranded
        if tops_equal and bottoms_equal and topology_equal and stranded_equal:
            return True
        else:
            return False

    def __repr__(self):
        '''String to print when object is called directly.'''
        parent = super(DNA, self).__repr__()
        display_bases = 40
        if len(self._sequence) < 90:
            bottom = self._bottom[::-1]
        else:
            rev_bottom = self._bottom[::-1]
            bottom = ''.join([rev_bottom[0:display_bases], ' ... ',
                              rev_bottom[-display_bases:]])
        first_line = '{} {}DNA:'.format(self.topology, self.stranded)
        to_print = '\n'.join([first_line, parent, bottom])
        return to_print

    def __setitem__(self, index, new_value):
        '''Sets value at index to new value.

        :param index: The index at which the sequence will be modified.
        :type index: int
        :param new_value: The new value at that index
        :type new_value: str or pymbt.DNA
        :returns: The current sequence with the sequence at `index` replaced
                  with `new_value`.
        :rtype: pymbt.DNA
        :raises: ValueError if `new_value` is '-'.

        '''
        new_value = str(new_value)
        if new_value == '-':
            raise ValueError('Cannot insert gap - split sequence instead.')
        # setitem on top strand
        super(DNA, self).__setitem__(index, new_value)
        # setitem on bottom strand
        if self.stranded == 'ds':
            sequence_list = list(self._bottom)[::-1]
            sequence_list[index] = str(DNA(new_value).reverse_complement())
            self._bottom = ''.join(sequence_list[::-1])
        else:
            self._bottom = '-' * len(self)


class RestrictionSite(object):
    '''Recognition site and properties of a restriction endonuclease.'''
    def __init__(self, recognition_site, cut_site, name=None):
        '''
        :param recognition_site: Input sequence.
        :type recognition_site: pymbt.DNA
        :param cut_site: 0-indexed indices where DNA is nicked (top, then
                         bottom strand). For an n-sized recognition site, there
                         are n + 1 positions at which to cut.
        :type cut_site: 2-tuple.
        :param name: Identifier of this restriction site
        :type name: str
        :returns: instance of pymbt.RestrictionSite

        '''
        self.recognition_site = recognition_site  # require DNA object
        # cutsite is indexed to leftmost base of restriction site
        self.cut_site = cut_site  # tuple of where top/bottom strands are cut
        # optional name
        self.name = name

    def is_palindrome(self):
        '''Report whether sequence is palindromic.

        :returns: Whether the restriction site is a palindrome.
        :rtype: bool

        '''
        return self.recognition_site.is_palindrome()

    def cuts_outside(self):
        '''Report whether the enzyme cuts outside its recognition site.
        Cutting at the very end of the site returns True.

        :returns: Whether the enzyme will cut outside its recognition site.
        :rtype: bool

        '''
        for index in self.cut_site:
            if index < 0 or index > len(self.recognition_site) + 1:
                return True
        return False

    def copy(self):
        '''Return copy of the restriction site.

        :returns: A safely editable copy of the current restriction site.
        :rtype: pymbt.RestrictionSite

        '''
        return RestrictionSite(self.recognition_site, self.cut_site,
                               self.name)

    def __repr__(self):
        '''Represent a restriction site.'''
        site = self.recognition_site
        cut_symbols = ('|', '|')
        if not self.cuts_outside():
            top_left = str(site[0:self.cut_site[0]])
            top_right = str(site[self.cut_site[0]:])
            top_w_cut = top_left + cut_symbols[0] + top_right

            bottom_left = site[0:self.cut_site[1]].reverse_complement()
            bottom_left = str(bottom_left)[::-1]
            bottom_right = site[self.cut_site[1]:].reverse_complement()
            bottom_right = str(bottom_right)[::-1]
            bottom_w_cut = bottom_left + cut_symbols[1] + bottom_right
        else:
            return '\n'.join([site.top() + ' {}'.format(self.cut_site),
                              site.bottom()])

        return '\n'.join([top_w_cut, bottom_w_cut])

    def __len__(self):
        '''Defines len operator.

        :returns: Length of the recognition site.
        :rtype: int

        '''
        return len(self.recognition_site)


class Primer(object):
    '''A DNA primer - ssDNA with tm, anneal, and optional overhang.'''
    def __init__(self, anneal, tm, overhang=None, name='', note=''):
        '''
        :param anneal: Annealing sequence
        :type anneal: pymbt.DNA
        :param overhang: Overhang sequence
        :type overhang: pymbt.DNA
        :param tm: melting temperature
        :type tm: float
        :param name: Optional name of the primer. Used when writing to csv with
                     seqio.write_primers.
        :type name: str
        :param note: Optional description to associate with the primer. Used
                     when writing to csv with seqio.write_primers.
        :type note: str
        :returns: pymbt.Primer instance.

        '''
        self.tm = tm
        self.anneal = anneal.to_ss()
        if overhang is not None:
            self.overhang = overhang.to_ss()
        else:
            self.overhang = DNA('', stranded='ss')
        self.name = name
        self.note = note

    def copy(self):
        '''Generate a Primer copy.

        :returns: A safely-editable copy of the current primer.
        :rtype: pymbt.DNA

        '''
        return type(self)(self.anneal, self.tm, overhang=self.overhang,
                          name=self.name, note=self.note)

    def primer(self):
        '''Produce full (overhang + annealing sequence) primer sequence.

        :returns: The DNA sequence of the primer.
        :rtype: pymbt.DNA

        '''
        return self.overhang + self.anneal

    def __repr__(self):
        '''Representation of a primer.'''
        if self.overhang:
            return 'Primer: {} Tm: {:.2f}'.format(self.overhang.top().lower() +
                                                  self.anneal.top(), self.tm)
        else:
            return 'Primer: {} Tm: {:.2f}'.format(self.anneal.top(), self.tm)

    def __str__(self):
        '''Coerce DNA object to string.

        :returns: A string of the full primer sequence.
        :rtype: str

        '''
        return str(self.primer())

    def __eq__(self, other):
        '''Define equality - sequences, topology, and strandedness are the
        same.

        :returns: Whether two primers have the same overhang and annealing
                  sequence.
        :rtype: bool

        '''
        anneal_equal = self.anneal == other.anneal
        overhang_equal = self.overhang == other.overhang
        if anneal_equal and overhang_equal:
            return True
        else:
            return False

    def __len__(self):
        '''Define len operator.

        :returns: The length of the full primer sequence.
        :rtype: int

        '''
        return len(self.primer())
