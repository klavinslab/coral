'''DNA object classes.'''
import os
import shutil
import subprocess
import tempfile
import coral.analysis
import coral.reaction
import coral.seqio
from ._sequence import Feature
from ._nucleicacid import NucleicAcid


class IPythonDisplayImportError(ImportError):
    '''Failed to import IPython display modules - display requires IPython'''


class ssDNA(NucleicAcid):
    '''ssDNA sequence.'''

    def __init__(self, sequence, circular=False, run_checks=True):
        super(ssDNA, self).__init__(sequence, 'dna', circular=circular,
                                    run_checks=run_checks, any_char='N')
        self.ds = False

    def copy(self):
        return type(self)(self.seq, circular=self.circular, run_checks=False)

    def to_ds(self):
        return DNA(self.seq, circular=self.circular, run_checks=False)


class DNA(object):
    '''dsDNA sequence.'''

    def __init__(self, dna, circular=False, features=None, run_checks=True,
                 bottom=None, name=None):
        '''
        :param dna: Input sequence (DNA).
        :type dna: str
        :param circular: The topology of the DNA - True for circular DNA, False
                         for linear DNA.
        :type circular: bool
        :param features: List of annotated features.
        :type features: list
        :param run_checks: Check inputs / formats (disabling increases speed):
                           alphabet check
                           case
        :type run_checks: bool
        :param bottom: String for the bottom sequence for manually setting
                       overhangs, fast copying.
        :type bottom: str
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
        self.material = 'dna'
        dna = dna.strip()
        self.top = ssDNA(dna, circular=circular, run_checks=run_checks)
        if bottom is None:
            # If bottom isn't auto-generated, expectation is that it will be
            # added manually
            self.bottom = self.top.reverse_complement()
        else:
            self.bottom = ssDNA(bottom, circular=circular, run_checks=False)

        if features is None:
            self.features = []
        else:
            self.features = features
        self.circular = circular
        self.ds = True
        if name is None:
            self.name = ''
        else:
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

    def copy(self):
        '''Create a copy of the current instance.

        :returns: A safely-editable copy of the current sequence.
        :rtype: coral.DNA

        '''
        # Significant performance improvements by skipping alphabet check
        features_copy = [feature.copy() for feature in self.features]
        copy = type(self)(self.top.seq, circular=self.circular,
                          features=features_copy, name=self.name,
                          bottom=self.bottom.seq, run_checks=False)
        return copy

    def circularize(self):
        '''Circularize linear DNA.

        :returns: A circularized version of the current sequence.
        :rtype: coral.DNA

        '''
        if self.top[-1].seq == '-' and self.bottom[0].seq == '-':
            raise ValueError('Cannot circularize - termini disconnected.')
        if self.bottom[-1].seq == '-' and self.top[0].seq == '-':
            raise ValueError('Cannot circularize - termini disconnected.')

        copy = self.copy()
        copy.circular = True
        copy.top.circular = True
        copy.bottom.circular = True
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
        require([\'{d3_cdn}\'], function(lib) {{
            window.data = {data};'''.format(div_id=div_id, d3_cdn=d3cdn,
                                            data=sequence_json)

        js_viz = '''
            d3sequence(window.data, \'{div_id}\')
        }});
        </script>
        '''.format(div_id=div_id)

        return HTML(html + js_databind + d3_plasmid_js + js_viz)

    def json(self):
        import json
        dna_json = {}
        dna_json['name'] = self.name
        dna_json['material'] = self.material
        dna_json['circular'] = self.circular
        dna_json['sequence'] = self.top.seq
        dna_json['bottom'] = self.bottom.seq

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

    def excise(self, feature):
        '''Removes feature from circular plasmid and linearizes. Automatically
        reorients at the base just after the feature. This operation is
        complementary to the .extract() method.

        :param feature_name: The feature to remove.
        :type feature_name: coral.Feature

        '''
        rotated = self.rotate_to_feature(feature)
        excised = rotated[feature.stop - feature.start:]

        return excised

    def extract(self, feature, remove_subfeatures=False):
        '''Extract a feature from the sequence. This operation is complementary
        to the .excise() method.

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
        copy.top, copy.bottom = copy.bottom, copy.top
        copy.features = [_flip_feature(f, len(self)) for f in copy.features]
        return copy

    def gc(self):
        '''Find the frequency of G and C in the current sequence. Calculates
        the exact frequency of Gs and Cs in each strand, so overhangs are
        counted only once.'''

        return (self.top.gc() + self.bottom.gc()) / 2

    def is_palindrome(self):
        if self.top.is_palindrome() and self.bottom.is_palindrome():
            return True
        else:
            return False

    def is_rotation(self, other):
        top_rotation = self.top.is_rotation(other.top)
        bottom_rotation = self.bottom.is_rotation(other.bottom)
        if top_rotation and bottom_rotation:
            return True
        else:
            return False

    def linearize(self, index=0):
        '''Linearize circular DNA at an index.

        :param index: index at which to linearize.
        :type index: int
        :returns: A linearized version of the current sequence.
        :rtype: coral.DNA
        :raises: ValueError if the input is linear DNA.

        '''
        if not self.circular:
            raise ValueError('Cannot relinearize linear DNA.')
        copy = self.copy()
        # Snip at the index
        if index:
            return copy[index:] + copy[:index]
        copy.circular = False
        copy.top.circular = False
        copy.bottom.circular = False

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
        top_matches = self.top.locate(pattern)
        bottom_matches = self.bottom.locate(pattern)

        return [top_matches, bottom_matches]

    def mw(self):
        '''Calculate the molecular weight.

        :returns: The molecular weight of the current sequence.
        :rtype: float

        '''
        return self.top.mw() + self.bottom.mw()

    def rotate(self, n):
        '''Rotate Sequence by n bases.

        :param n: Number of bases to rotate.
        :type n: int
        :returns: The current sequence reoriented at `index`.
        :rtype: coral.DNA
        :raises: ValueError if applied to linear sequence or `index` is
                 negative.

        '''
        if not self.circular and n != 0:
            raise ValueError('Cannot rotate linear DNA')
        else:
            copy = self.copy()
            copy.top = self.top.rotate(n)
            copy.bottom = self.bottom.rotate(-n)
            copy.features = []
            for feature in self.features:
                feature_copy = feature.copy()
                feature_copy.move(n)
                # Adjust the start/stop if we move over the origin
                feature_copy.start = feature_copy.start % len(self)
                feature_copy.stop = feature_copy.stop % len(self)
                copy.features.append(feature_copy)

            return copy.circularize()

    def rotate_to(self, index):
        '''Orient DNA to index (only applies to circular DNA).

        :param index: DNA position at which to re-zero the DNA.
        :type index: int
        :returns: The current sequence reoriented at `index`.
        :rtype: coral.DNA
        :raises: ValueError if applied to linear sequence or `index` is
                 negative.

        '''
        return self.rotate(-index)

    def rotate_to_feature(self, feature):
        '''Reorient the DNA based on a feature it contains (circular DNA only).

        :param feature: A feature.
        :type feature: coral.Feature
        :returns: The current sequence reoriented at the start index of a
                  unique feature matching `featurename`.
        :rtype: coral.DNA
        :raises: ValueError if there is no feature of `featurename` or
                 more than one feature matches `featurename`.

        '''
        return self.rotate_to(feature.start)

    def reverse_complement(self):
        '''Reverse complement the DNA.

        :returns: A reverse-complemented instance of the current sequence.
        :rtype: coral.DNA

        '''
        copy = self.copy()
        # Note: if sequence is double-stranded, swapping strand is basically
        # (but not entirely) the same thing - gaps affect accuracy.
        copy.top = self.top.reverse_complement()
        copy.bottom = self.bottom.reverse_complement()

        # Remove all features - the reverse complement isn't flip!
        copy.features = []

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

    def to_feature(self, name=None, feature_type='misc_feature'):
        '''Create a feature from the current object.

        :param name: Name for the new feature. Must be specified if the DNA
                     instance has no .name attribute.
        :type name: str
        :param feature_type: The type of feature (genbank standard).
        :type feature_type: str

        '''
        if name is None:
            if not self.name:
                raise ValueError('name attribute missing from DNA instance'
                                 ' and arguments')
            name = self.name
        return Feature(name, start=0, stop=len(self),
                       feature_type=feature_type)

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
        if self.circular or other.circular:
            raise Exception('Can only add linear DNA.')

        discontinuity = [False, False]
        if len(self) != 0 and len(other) != 0:
            # If either is empty, let things proceed anyways
            discontinuity[0] = (self.top[-1] == '-' and
                                other.bottom[-1] == '-')
            discontinuity[1] = (self.bottom[0] == '-' and
                                other.top[0] == '-')

        for_discontinuity = discontinuity[0]
        rev_discontinuity = discontinuity[1]

        if for_discontinuity or rev_discontinuity:
            msg = 'Concatenated DNA would be discontinuous.'
            raise Exception(msg)

        tops = self.top.seq + other.top.seq
        bottoms = other.bottom.seq + self.bottom.seq

        self_features = [feature.copy() for feature in self.features]
        other_features = [feature.copy() for feature in other.features]
        for feature in other_features:
            feature.move(len(self))
        features = self_features + other_features

        new_instance = DNA(tops, circular=False, run_checks=False,
                           bottom=self.bottom.seq, features=features)
        new_instance.bottom = ssDNA(bottoms)

        return new_instance

    def __contains__(self, query):
        '''Defines `query in sequence` operator.

        :param query: query string or DNA sequence
        :type query: str or coral.DNA

        '''
        if query in self.top or query in self.bottom:
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
        del self.top[key]
        del self.bottom[-key]

    def __getitem__(self, key):
        '''Index and slice sequences.

        :param key: int or slice object for subsetting.
        :type key: int or slice object
        :returns: A subsequence matching the slice (`key`).
        :rtype: coral.DNA

        '''
        # Adjust features
        def in_slice(feature, key, circular=False):
            '''Classify a coral.Feature object as within a slice or not.

            :param feature: The feature to test.
            :type feature: coral.Feature
            :param key: A slice key.
            :type key: slice
            :param circular: Sets whether the parent sequence is circular.
            :type circular: bool

            '''
            if circular:
                # FIXME: once circular slicing w/ negative indices is
                # implemented, need to account for that here.
                # Decide whether a feature's location is within a slice's
                # coordinates
                if key.start is None:
                    if key.stop is None:
                        # Both are none - copying whole sequence using [:]
                        return True
                    else:
                        # The slice looks like [:stop], i.e. [0:stop]
                        if feature.stop > key.stop:
                            # feature extends beyond key stop
                            return False
                        else:
                            if feature.start > key.stop:
                                # If feature extends over origin, remove it
                                return False
                            else:
                                # Feature is between 0 and key.stop
                                return True
                else:
                    if key.stop is None:
                        # The slice looks like [start:] i.e.[start:length]
                        key_start = key.start % len(self)
                        if feature.start < key_start:
                            return False
                        else:
                            if feature.stop < key_start:
                                # Feature extends over origin - remove it
                                return False
                            else:
                                # Feature is between key.start and end of seq
                                return True
                    else:
                        # The slice looks like [key.start:key.stop]
                        if feature.start < key.start or \
                           feature.stop > key.stop:
                            return False
                        else:
                            return True
            if key.start is None:
                if key.stop is None:
                    # Copying whole sequence with [:]
                    return True
                else:
                    # Slice looks like [:key.stop]
                    if feature.stop > key.stop:
                        # Feature ends after key.stop
                        return False
                    else:
                        # Feature ends before key.stop
                        return True
            else:
                if key.stop is None:
                    # Slice looks like [key.start:]
                    if feature.start < key.start:
                        return False
                    else:
                        return True
                else:
                    # The slice looks like [key.start:key.stop]
                    if feature.start < key.start or \
                       feature.stop > key.stop:
                        return False
                    else:
                        return True

        saved_features = []
        if self.features:
            if isinstance(key, slice):
                # If a slice, remove stuff that isn't in the slice and
                # adjust feature starts/stops
                if key.step == 1 or key.step is None:
                    for feature in self.features:
                        if in_slice(feature, key, self.circular):
                            feature_copy = feature.copy()
                            if key.start:
                                feature_copy.move(-(key.start % len(self)))
                            saved_features.append(feature_copy)
                else:
                    # Don't copy any features - a non-1 step size should
                    # not leave any features instace.
                    pass
            else:
                for feature in self.features:
                    if feature.start == feature.stop == key:
                        feature_copy = feature.copy()
                        feature_copy.move(key)
                        saved_features.append(feature_copy)

        # Run __getitem__ on top and bottom sequences
        new_top = self.top.__getitem__(key)
        # TODO: improve efficiency of bottom's __getitem__ (this is a slow
        # strategy)
        new_bottom = self.bottom[::-1][key][::-1]

        copy = type(self)(new_top.seq, circular=False,
                          features=saved_features, bottom=self.bottom.seq,
                          run_checks=False)
        copy.bottom = new_bottom

        return copy

    def __eq__(self, other):
        '''Define equality - sequences and strandedness are the.

        :returns: Whether current sequence's (Watson and Crick), topology,
                  and strandedness are equivalent to those of another sequence.
        :rtype: bool

        '''
        tops_equal = self.top == other.top
        bottoms_equal = self.bottom == other.bottom
        if tops_equal and bottoms_equal:
            return True
        else:
            return False

    def __hash__(self):
        # Enables the use of functions like set() - hash unique attributes
        return hash(self.top.seq + self.bottom.seq + str(self.circular))

    def __len__(self):
        return len(self.top)

    def __mul__(self, n):
        copy = self.copy()
        copy.top = self.top * n
        copy.bottom = self.bottom * n
        copy.features = []
        for i in range(n):
            for feature in self.features:
                new_feature = feature.copy()
                new_feature.start += i * len(self)
                new_feature.stop += i * len(self)
                copy.features = []

        return copy

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
        top = self.top.__repr__()
        bottom = self.bottom.__repr__()[::-1]

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

        self.top[index] = new_value
        self.bottom = self.bottom[::-1][index][::-1]

    def __str__(self):
        return self.top.seq


def _flip_feature(self, feature, parent_len):
    '''Adjust a feature's location when flipping DNA.

    :param feature: The feature to flip.
    :type feature: coral.Feature
    :param parent_len: The length of the sequence to which the feature belongs.
    :type parent_len: int

    '''
    copy = feature.copy()
    # Put on the other strand
    if copy.strand == 0:
        copy.strand = 1
    else:
        copy.strand = 0
    # Adjust locations - guarantee that start is always less than end
    copy.start = parent_len - copy.start
    copy.stop = parent_len - copy.stop
    copy.start, copy.stop = copy.stop, copy.start
    return copy


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
            return '\n'.join([site.top + ' {}'.format(self.cut_site),
                              site.bottom])

        return '\n'.join([top_w_cut, bottom_w_cut])

    def __len__(self):
        '''Defines len operator.

        :returns: Length of the recognition site.
        :rtype: int

        '''
        return len(self.recognition_site)


class Primer(object):
    '''A DNA primer - ssDNA with tm, anneal, and optional overhang.'''

    def __init__(self, anneal, tm=None, tm_method='cloning', overhang=None,
                 name='', note=''):
        '''
        :param anneal: Annealing sequence
        :type anneal: coral.ssDNA
        :param overhang: Overhang sequence
        :type overhang: coral.ssDNA
        :param tm: Melting temperature (allow manually setting the melting
                   temp, auto-calculation won't happen).
        :type tm: float
        :param tm_method: Method to use when calling cr.analysis.tm.
        :type tm_method: str
        :param name: Optional name of the primer. Used when writing to csv with
                     seqio.write_primers.
        :type name: str
        :param note: Optional description to associate with the primer. Used
                     when writing to csv with seqio.write_primers.
        :type note: str
        :returns: coral.Primer instance.

        '''
        if tm is not None:
            self.tm = tm
        else:
            self.tm = anneal.tm()

        # TODO: make into property so tm updates automatically
        if anneal.ds:
            self.anneal = anneal.top
        else:
            self.anneal = anneal

        if overhang is not None:
            if overhang.ds:
                self.overhang = overhang.top
            else:
                self.overhang = overhang
        else:
            self.overhang = ssDNA('')

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
        '''Produce (ss) DNA of the overhang + anneal sequence.

        :returns: The (ss) DNA sequence of the primer.
        :rtype: coral.ssDNA

        '''
        return self.overhang + self.anneal

    def to_ds(self):
        '''Produces (ds) DNA of the overhang + anneal sequence.

        :returns: The DNA sequence of the primer.
        :rtype: coral.DNA

        '''
        return self.primer().to_ds()

    def __repr__(self):
        '''Representation of a primer.'''
        anneal = self.anneal.seq.upper()
        if self.overhang:
            overhang = self.overhang.seq.lower()
            return 'Primer: {} Tm: {:.2f}'.format(overhang + anneal, self.tm)
        else:
            return 'Primer: {} Tm: {:.2f}'.format(anneal, self.tm)

    def __str__(self):
        '''Coerce Primer object to string.

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
