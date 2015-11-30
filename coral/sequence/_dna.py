'''DNA object classes.'''
import collections
import os
import shutil
import subprocess
import tempfile
import coral.analysis
import coral.reaction
import coral.seqio
from ._sequence import NucleicAcidSequence


class IPythonDisplayImportError(ImportError):
    '''Failed to import IPython display modules - display requires IPython'''


class DNA(object):
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
        :param name: Optional name field for your DNA sequence.
        :type name: str
        :returns: coral.DNA instance.
        :rtype: coral.DNA
        :raises: ValueError if an element of `features` isn't of type
                 coral.Feature.
                 ValueError if top and bottom strands have different lengths.
                 ValueError if top and bottom strands are not complementary.

        '''
        # TODO: accept sequences in general by running str() on it
        self._top = NucleicAcidSequence(dna, 'dna', run_checks=run_checks)

        if stranded == 'ss':
            self._bottom = NucleicAcidSequence('-' * len(self._top), 'dna',
                                               run_checks=False)
        else:
            if bottom is None:
                self._bottom = self._top.reverse_complement()
            else:
                if len(bottom) != len(self._top):
                    msg = 'Top and bottom strands are difference lengths.'
                    raise ValueError(msg)
                self._bottom = NucleicAcidSequence(bottom, 'dna',
                                                   run_checks=run_checks)
        if features is None:
            self.features = []
        else:
            self.features = features
        self.topology = topology
        self.stranded = stranded
        self.name = name

    def ape(self, ape_path=None):
        '''Open in ApE if `ApE` is in your command line path.'''
        # TODO: simplify - make ApE look in PATH only
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
        coral.seqio.write_dna(self, filename)
        process = subprocess.Popen([cmd, filename])
        # Block until window is closed
        try:
            process.wait()
            shutil.rmtree(tmp)
        except KeyboardInterrupt:
            shutil.rmtree(tmp)

    def bottom(self):
        '''Return the NucleicAcidSequence object of the Crick (bottom) strand.

        :returns: The Crick strand.
        :rtype: str

        '''
        return self._bottom

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely-editable copy of the current sequence.
        :rtype: coral.DNA

        '''
        # Significant performance improvements by skipping alphabet check
        features_copy = [feature.copy() for feature in self.features]
        return type(self)(str(self._top), bottom=str(self._bottom),
                          topology=self.topology, stranded=self.stranded,
                          features=features_copy, name=self.name,
                          run_checks=False)

    def circularize(self):
        '''Circularize linear DNA.

        :returns: A circularized version of the current sequence.
        :rtype: coral.DNA

        '''
        if self._top[-1].seq == '-' and self._bottom[0].seq == '-':
            raise ValueError('Cannot circularize - termini disconnected.')
        if self._bottom[-1].seq == '-' and self._top[0].seq == '-':
            raise ValueError('Cannot circularize - termini disconnected.')

        copy = self.copy()
        copy.topology = 'circular'
        return copy

    def display(self):
        '''Display a visualization of the sequence in an IPython notebook.'''
        try:
            from IPython.display import HTML
            import uuid
        except ImportError:
            raise IPythonDisplayImportError

        sequence_json = self.json()

        d3cdn = '//d3js.org/d3.v3.min.js'
        div_id = 'sequence_{}'.format(uuid.uuid1())

        cur_dir = os.path.abspath(os.path.dirname(__file__))
        d3_plasmid_path = os.path.join(cur_dir, 'd3-plasmid.js')
        with open(d3_plasmid_path) as f:
            d3_plasmid_js = f.read()

        html = '<div id={div_id}></div>'.format(div_id=div_id)
        js_databind = '''
        <script>
        require(["{d3_cdn}"], function(lib) {{
            window.data = {data};'''.format(div_id=div_id, d3_cdn=d3cdn,
                                            data=sequence_json)

        js_viz = '''
            d3sequence(window.data, "{div_id}")
        }});
        </script>
        '''.format(div_id=div_id)

        return HTML(html + js_databind + d3_plasmid_js + js_viz)

    def json(self):
        import json
        dna_json = {}
        dna_json['name'] = self.name
        dna_json['material'] = 'DNA'
        dna_json['topology'] = self.topology
        dna_json['sequence'] = str(self.top())
        dna_json['bottom'] = str(self.bottom())

        features_json = []
        for feature in self.features:
            feature_json = {}
            feature_json['start'] = feature.start
            feature_json['stop'] = feature.stop
            feature_json['type'] = feature.feature_type
            feature_json['name'] = feature.name
            if 'ApEinfo_fwdcolor' in feature.qualifiers:
                feature_json['color'] = feature.qualifiers['ApEinfo_fwdcolor']

            features_json.append(feature_json)
        dna_json['features'] = features_json

        return json.dumps(dna_json)

    def extract(self, feature, remove_subfeatures=False):
        '''Extract a feature from the sequence.

        :param feature: Feature object.
        :type feature: coral.sequence.Feature
        :param remove_subfeatures: Remove all features in the extracted
                                   sequence aside from the input feature.
        :type remove_subfeatures: bool
        :returns: A subsequence from start to stop of the feature.

        '''
        extracted = self[feature.start:feature.stop]
        # Turn gaps into Ns or Xs
        for gap in feature.gaps:
            for i in range(*gap):
                extracted[i] = self._any_char
        if remove_subfeatures:
            # Keep only the feature specified
            extracted.features = [feature]
        # Update feature locations
        # copy them
        for feature in extracted.features:
            feature.move(-feature.start)
        return extracted

    def flip(self):
        '''Flip the DNA - swap the top and bottom strands.

        :returns: Flipped DNA (bottom strand is now top strand, etc.).
        :rtype: coral.DNA

        '''
        copy = self.copy()
        copy._top, copy._bottom = copy._bottom, copy._top
        return copy

    def gc(self):
        '''Find the frequency of G and C in the current sequence.'''

        gc_n = len([base for base in self if str(base) == 'C' or
                    str(base) == 'G'])
        return float(gc_n) / len(self)

    def is_palindrome(self):
        if self._top.is_palindrome() and self._bottom.is_palindrome():
            return True
        else:
            return False

    def is_rotation(self, other):
        if len(self) != len(other):
            return False

        for i in range(len(self)):
            if self.rotate(i) == other:
                return True

        # Check reverse complement
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
        :rtype: coral.DNA
        :raises: ValueError if the input is linear DNA.

        '''
        if self.topology == 'linear':
            raise ValueError('Cannot relinearize linear DNA.')
        copy = self.copy()
        copy.topology = 'linear'
        copy = copy[index:] + copy[:index]
        return copy

    def locate(self, pattern):
        '''Find sequences matching a pattern. For a circular sequence, the
        search extends over the origin.

        :param pattern: str or NucleicAcidSequence for which to find matches.
        :type pattern: str or coral.DNA
        :returns: A list of top and bottom strand indices of matches.
        :rtype: list of lists of indices (ints)
        :raises: ValueError if the pattern is longer than either the input
                 sequence (for linear DNA) or twice as long as the input
                 sequence (for circular DNA).

        '''
        if self.topology == 'circular':
            if len(pattern) >= 2 * len(self):
                raise ValueError('Search pattern longer than searchable ' +
                                 'sequence.')
            top = self._top + self._top[:len(pattern) - 1]
            bottom = self._bottom + self._bottom[:len(pattern) - 1]
        else:
            if len(pattern) > len(self):
                raise ValueError('Search pattern longer than searchable ' +
                                 'sequence.')
            top = self._top
            bottom = self._bottom

        top_matches = top.locate(pattern)
        bottom_matches = bottom.locate(pattern)

        # Adjust indices if doing circular search
        def reindex(location_indices):
            for index in location_indices:
                if index > len(self):
                    index -= len(self)
            return location_indices

        if self.topology == 'circular' and len(pattern) > 1:
            top_matches = reindex(top_matches)
            bottom_matches = reindex(bottom_matches)

        return [top_matches, bottom_matches]

    def mw(self):
        '''Calculate the molecular weight.

        :returns: The molecular weight of the current sequence.
        :rtype: float

        '''
        bases = str(self._top) + str(self._bottom)
        counter = collections.Counter(bases.lower())
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
        :rtype: coral.DNA
        :raises: ValueError if applied to linear sequence or `index` is
                 negative.

        '''
        if self.topology == 'linear' and index != 0:
            raise ValueError('Cannot rotate linear DNA')
        else:
            return (self[index:] + self[:index]).circularize()

    def rotate_by_feature(self, feature):
        '''Reorient the DNA based on a feature it contains (circular DNA only).

        :param feature: A feature.
        :type feature: coral.Feature
        :returns: The current sequence reoriented at the start index of a
                  unique feature matching `featurename`.
        :rtype: coral.DNA
        :raises: ValueError if there is no feature of `featurename` or
                 more than one feature matches `featurename`.

        '''
        return self.rotate(feature.start)

    def reverse_complement(self):
        '''Reverse complement the DNA.

        :returns: A reverse-complemented instance of the current sequence.
        :rtype: coral.DNA

        '''
        copy = self.copy()
        # Note: if sequence is double-stranded, swapping strand is basically
        # (but not entirely) the same thing - gaps affect accuracy.
        copy._top = self._top.reverse_complement()
        copy._bottom = self._bottom.reverse_complement()

        # Adjust feature locations (invert)
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

    def select_features(self, term, by='name', fuzzy=False):
        '''Select features from the features list based on feature name,
           gene, or locus tag.
           :param term: Search term.
           :type term: str
           :param by: Feature attribute to search by. Options are 'name',
                      'gene', and 'locus_tag'.
           :type by: str
           :param fuzzy: If True, search becomes case-insensitive and will also
                         find substrings - e.g. if fuzzy search is enabled, a
                         search for 'gfp' would return a hit for a feature
                         named 'GFP_seq'.
           :type fuzzy: bool
           :returns: A list of features matched by the search.
           :rtype: list
        '''
        features = []
        if fuzzy:
            fuzzy_term = term.lower()
            for feature in self.features:
                if fuzzy_term in feature.__getattribute__(by).lower():
                    features.append(feature)
        else:
            for feature in self.features:
                if feature.__getattribute__(by) == term:
                    features.append(feature)

        return features

    def tm(self, parameters='cloning'):
        '''Find the melting temperature.

        :param parameters: The tm method to use (cloning, santalucia98,
                       breslauer86)
        :type parameters: str

        '''
        return coral.analysis.tm(self, parameters=parameters)

    def to_ss(self):
        '''Produce single stranded version of the current sequence.

        :returns: The current sequence, converted to ssDNA.
        :rtype: coral.DNA

        '''
        copy = self.copy()

        # Do nothing if already single-stranded
        if self.stranded == 'ss':
            return copy

        copy._bottom = NucleicAcidSequence('-' * len(copy), 'dna',
                                           run_checks=False)
        for top, bottom in zip(copy._top, reversed(copy._bottom)):
            if top == bottom == '-':
                raise ValueError('Coercing to single-stranded would ' +
                                 'introduce a double stranded break.')
        copy.stranded = 'ss'

        return copy

    def to_ds(self):
        '''Produce double stranded version of the current sequence.

        :returns: The current sequence, converted to dsDNA.
        :rtype: coral.DNA

        '''
        copy = self.copy()
        # Do nothing if already set
        if self.stranded == 'ds':
            return copy

        # Find strand that's all gaps (if ss this should be the case)
        reverse_seq = self.reverse_complement()
        if all([char == '-' for char in self._top]):
            copy._top = reverse_seq._bottom
        elif all([char == '-' for char in self._bottom]):
            copy._bottom = reverse_seq._top
        copy.stranded = 'ds'

        return copy

    def top(self):
        '''Return the NucleicAcidSequence object of the Watson (top) strand.

        :returns: The Watson strand.
        :rtype: str

        '''
        # TODO: Make .top and .bottom properties so they can be gotten/set
        # but also have docstrings
        return self._top

    def transcribe(self):
        '''Transcribe into RNA.

        :returns: An RNA sequence transcribed from the current DNA sequence.
        :rtype: coral.RNA

        '''
        return coral.reaction.transcribe(self)

    def __add__(self, other):
        '''Add DNA together.

        :param other: instance to be added to.
        :type other: compatible sequence object (currently only DNA).
        :returns: Concatenated DNA sequence.
        :rtype: coral.DNA
        :raises: Exception if either sequence is circular.
                 Exception if concatenating a sequence with overhangs would
                 create a discontinuity.

        '''
        if self.topology == 'circular' or other.topology == 'circular':
            raise Exception('Can only add linear DNA.')

        discontinuity = [False, False]
        if len(self) != 0 and len(other) != 0:
            # If either is empty, let things proceed anyways
            discontinuity[0] = (self._top[-1] == '-' and
                                other._bottom[-1] == '-')
            discontinuity[1] = (self._bottom[0] == '-' and
                                other._top[0] == '-')

        for_discontinuity = discontinuity[0]
        rev_discontinuity = discontinuity[1]

        if for_discontinuity or rev_discontinuity:
            msg = 'Concatenated DNA would be discontinuous.'
            raise Exception(msg)

        if self.stranded == 'ds' or other.stranded == 'ds':
            stranded = 'ds'
        else:
            stranded = 'ss'

        tops = str(self._top + other._top)
        bottoms = str(other._bottom + self._bottom)

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
        :type query: str or coral.DNA

        '''
        if query in self._top or query in self._bottom:
            return True
        else:
            return False

    def __delitem__(self, key):
        '''Delete sequence at an key.

        :param key: key to delete
        :type key: int
        :returns: The current sequence with the base at `key` removed.
        :rtype: coral.DNA

        '''
        del self._top[key]
        del self._bottom[-key]

    def __getitem__(self, key):
        '''Index and slice sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object
        :returns: A subsequence matching the slice (`key`).
        :rtype: coral.DNA

        '''
        # Adjust features
        def in_slice(feature):
            if key.start is not None and feature.start < key.start:
                return False
            elif key.stop is not None and feature.stop > key.stop:
                return False
            else:
                return True

        saved_features = []
        if self.features:
            if isinstance(key, slice):
                # If a slice, remove stuff that isn't in the slide and
                # adjust feature starts/stops
                if key.step == 1 or key.step is None:
                    for feature in self.features:
                        if in_slice(feature):
                            if key.start:
                                feature.move(-key.start)
                            saved_features.append(feature.copy())
            else:
                for feature in self.features:
                    if feature.start == feature.stop == key:
                        saved_features.append(feature.copy())
                        saved_features[-1].move(key)

        # Run __getitem__ on top and bottom sequences
        new_top = self._top.__getitem__(key)
        # TODO: improve efficiency of bottom's __getitem__ (this is a slow
        # strategy)
        new_bottom = self._bottom[::-1][key][::-1]

        copy = type(self)(str(new_top), bottom=str(new_bottom),
                          topology='linear', stranded=self.stranded,
                          features=saved_features, name=self.name,
                          run_checks=False)

        copy.topology = 'linear'

        return copy

    def __eq__(self, other):
        '''Define equality - sequences and strandedness are the.

        :returns: Whether current sequence's (Watson and Crick), topology,
                  and strandedness are equivalent to those of another sequence.
        :rtype: bool

        '''
        tops_equal = self._top == other._top
        bottoms_equal = self._bottom == other._bottom
        stranded_equal = self.stranded == other.stranded
        if tops_equal and bottoms_equal and stranded_equal:
            return True
        else:
            return False

    def __len__(self):
        return len(self._top)

    def __mul__(self, n):
        new_top = self._top * n
        new_bottom = self._bottom * n
        features_copy = [feature.copy() for feature in self.features]

        return type(self)(str(new_top), str(new_bottom),
                          topology=self.topology, stranded=self.stranded,
                          features=features_copy, name=self.name,
                          run_checks=False)

    def __radd__(self, other):
        if other == 0 or other is None:
            # For compatibility with sum()
            return self
        elif type(self) != type(other):
            try:
                other = type(self)(other)
            except AttributeError:
                raise TypeError('Cannot add {} to {}'.format(self, other))
        return self + other

    def __repr__(self):
        '''String to print when object is called directly.'''
        top = self._top.__repr__()
        bottom = self._bottom.__repr__()[::-1]

        return '{}\n{}'.format(top, bottom)

    def __setitem__(self, index, new_value):
        '''Sets value at index to new value.

        :param index: The index at which the sequence will be modified.
        :type index: int
        :param new_value: The new value at that index
        :type new_value: str or coral.DNA
        :returns: The current sequence with the sequence at `index` replaced
                  with `new_value`.
        :rtype: coral.DNA
        :raises: ValueError if `new_value` is '-'.

        '''
        if new_value == '-':
            raise ValueError('Cannot insert gap - split sequence instead.')

        self._top[index] = new_value
        self._bottom = self._bottom[::-1][index][::-1]

    def __str__(self):
        return str(self._top)


class RestrictionSite(object):
    '''Recognition site and properties of a restriction endonuclease.'''
    def __init__(self, recognition_site, cut_site, name=None):
        '''
        :param recognition_site: Input sequence.
        :type recognition_site: coral.DNA
        :param cut_site: 0-indexed indices where DNA is nicked (top, then
                         bottom strand). For an n-sized recognition site, there
                         are n + 1 positions at which to cut.
        :type cut_site: 2-tuple.
        :param name: Identifier of this restriction site
        :type name: str
        :returns: instance of coral.RestrictionSite

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
        :rtype: coral.RestrictionSite

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
        :type anneal: coral.DNA
        :param overhang: Overhang sequence
        :type overhang: coral.DNA
        :param tm: melting temperature
        :type tm: float
        :param name: Optional name of the primer. Used when writing to csv with
                     seqio.write_primers.
        :type name: str
        :param note: Optional description to associate with the primer. Used
                     when writing to csv with seqio.write_primers.
        :type note: str
        :returns: coral.Primer instance.

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
        :rtype: coral.DNA

        '''
        return type(self)(self.anneal, self.tm, overhang=self.overhang,
                          name=self.name, note=self.note)

    def primer(self):
        '''Produce full (overhang + annealing sequence) primer sequence.

        :returns: The DNA sequence of the primer.
        :rtype: coral.DNA

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
