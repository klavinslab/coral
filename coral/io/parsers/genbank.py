'''Read and write DNA sequences.'''
import coral as cr
import re


def parse_genbank(handle):
    '''Parse DNA from a genbank file string.

    :param handle: File handle (i.e. the output of open()).
    :type lines: file
    :returns: DNA sequence.
    :rtype: coral.DNA

    '''
    read = handle.read()

    # In genbank format, only keywords can be at the beginning of a line and
    # are made up of all-caps characters
    # keywords = re.findall('^([A-Z]+)', read, flags=re.M)
    keywords_search = re.findall('^([A-Z]+)(.+?)(?=\n[A-Z]+|\n//)', read,
                                 flags=re.S | re.M)
    keywords = {key: key + val for key, val in keywords_search}

    # The ORIGIN keyword (sequence) is always last
    sequence = cr.DNA('')
    for line in keywords['ORIGIN'].split('\n'):
        if line:
            sequence += cr.DNA(''.join(line.split()[1:]))

    # Read in LOCUS information
    locus = keywords['LOCUS']
    sequence.name = locus[12:28].strip()
    # bp = locus[1]
    if 'ss' in locus[44:46]:
        sequence = sequence.to_ss()

    if 'circular' in locus[53:63]:
        sequence.circular = True
    else:
        sequence.circular = False

    # # Read in DEFINITION
    # definition = keyword_dict['DEFINITION'].strip().split()[1]
    # sequence.definition = definition

    # # Read in FEATURES
    # Features get grouped by a key (6 chars in) and end either on the next
    # key or the next major all-caps keyword section
    feature_groups = re.findall('\n(\s{5}[\w\'-]+.*?)(?=\n\s{5}[\w\'-]+|\Z)',
                                keywords['FEATURES'], flags=re.S)

    for feature_group in feature_groups:
        feature = process_feature(feature_group)
        if feature is not None:
            sequence.features.append(feature)

    return sequence


def process_location(loc_str):
    '''Given the string format of a genbank feature location, return the start,
    stop, strand, and gaps information.

    :param location_str: location information for a given genbank feature, in
                         raw string format.
    :type location_str: str
    :returns: start, stop, strand, and gaps information as a tuple.
    :rtype: tuple of integer, integer, integer, and a list of lists integer
            pairs (or None)

    '''
    if loc_str.startswith('complement'):
        # The complement () operator seems to always be the outermost operator
        strand = 1
        loc_str = loc_str.lstrip('complement(')[:-1]
    else:
        strand = 0

    def get_loc(loc_str):
        '''If the location has no order or join operator, this returns its
        start + stop information.'''
        single_number = re.search('^[0-9]+$', loc_str)
        if single_number:
            # The location value can be a single integer
            start = int(loc_str) - 1
            stop = start + 1
            return start, stop
        between = re.search('^[0-9]+\^[0-9]+$', loc_str)
        if between:
            # The feature is 'between' two locations, e.g. an endonuclease cut
            # site. Coral does not support features that have no length and
            # will instead simply list this as the latter position - e.g. if
            # the cut site is between bases 22 and 23,the feature will be
            # located at position 23.
            start = int(re.search('[0-9]+\^', loc_str)) - 1
            stop = start + 1
            return start, stop
        range_single = re.search('^([0-9]+)\.([0-9]+)$', loc_str)
        if range_single:
            # The feature is somewhere in the range specified. This is not
            # supported by coral and the feature will simply be described by
            # the entire range
            start = int(range_single.group(1)) - 1
            stop = int(range_single.group(2))
            return start, stop
        contiguous = re.search('^([0-9]+)\.\.([0-9]+)$', loc_str)
        if contiguous:
            start = int(contiguous.group(1)) - 1
            stop = int(contiguous.group(2))
            return start, stop

        raise ValueError('Feature did not fall into any location '
                         'category: {}'.format(loc_str))

    order = re.search('^order\((.*)\)', loc_str)
    join = re.search('^join\((.*)\)', loc_str)
    if order or join:
        if order:
            contents = order.group(1)
        if join:
            contents = join.group(1)
        loci = contents.split(',')
        locs = [get_loc(l) for l in loci]
        start = min(locs, key=lambda x: x[0])[0]
        stop = max(locs, key=lambda x: x[1])[1]

        locs_sorted = sorted(locs, key=lambda x: x[0])
        gaps = []
        for region1, region2 in zip(locs_sorted[:-1], locs_sorted[1:]):
            gaps.append([region1[1] - 1, region2[0] + 1])
    else:
        start, stop = get_loc(loc_str)
        gaps = None

    return start, stop, strand, gaps


def process_feature(feature_str):
    '''Create a coral.Feature instance from a genbank feature section (e.g.
    a single CDS feature and its qualifiers.

    :param feature_str: The raw string from a genbank file representing a
                        feature, with newlines intact, e.g. CDS 4\n/label="1".
    :type feature_str: str
    :returns: coral.Feature or None if feature has no /label qualifier.
    :rtype: coral.Feature or NoneType

    '''
    line0 = feature_str.split('\n')[0].split()
    feature_type = line0[0]
    start, stop, strand, gaps = process_location(line0[1])

    # Qualifiers begin at position 22 and begin with \. They end either
    # when another qualifier appears (\), a new feature appears, or the
    # next major KEYWORD section begins
    qual_regex = '(\s{21}/[\w\'-]+.*?)(?=\n\s{21}/[\w\'-]+|\Z)'
    unprocessed_quals = re.findall(qual_regex, feature_str, flags=re.S)

    qualifiers = {}
    for qual_lines in unprocessed_quals:
        # Newlines within a qualifiers should be removed along with left-pad
        # whitespace
        qual_line = ' '.join([l.strip() for l in qual_lines.split('\n')])
        # Qualifiers have a name (e.g. /name) and optionally a value following
        # an equals sign (e.g. \name="value"). Values do not need to have
        # quotes around them and are used inconsistently. Qualifiers with no
        # values act like boolean flags (e.g. mark as /pseudo for pseudogene).

        # Does the qualifier have a value? This only occurs if there's
        # a \name= string
        name_value = re.search('^/(.*)(?==)=(.*)', qual_line)
        if name_value:
            qual_name = name_value.group(1)
            qual_value = name_value.group(2).strip('"')
        else:
            qual_name = re.search('^/(.*)', qual_line).group(1)
            qual_value = True
        qualifiers[qual_name] = qual_value

    # Only save features that have labels (i.e. names)
    if 'label' in qualifiers:
        feature_name = qualifiers['label']
        gene = qualifiers['gene'] if 'gene' in qualifiers else None
        locus_tag = None
        if 'locus_tag' in qualifiers:
            locus_tag = qualifiers['locus_tag']
        feature = cr.Feature(feature_name, start, stop,
                             feature_type=feature_type, gene=gene,
                             locus_tag=locus_tag, qualifiers=qualifiers,
                             strand=strand, gaps=gaps)

        return feature
