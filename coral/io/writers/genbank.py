'''Write genbank sequences.'''
import datetime
import math
import os
import textwrap


def write_genbank(sequence, handle):
    '''Write coral.DNA, coral.ssDNA, or coral.RNA objects to a genbank file.
    The spec used for the output format is the official Genbank file format
    spec available at ftp://ftp.ncbi.nih.gov/genbank/gbrel.txt.

    :param sequence: The sequence to write.
    :type sequence: coral.DNA, coral.ssDNA, or coral.RNA
    :param handle: File handle (i.e. the output of open()).
    :type handle: file

    '''
    lines = []
    # Construct the LOCUS keyword. Requirements:
    # 1. LOCUS identifier starts at position 13 and is at most 16 chars long
    # 2. Last digit for number of bp end on position 40.
    locus = []
    locus.append('LOCUS       ')
    if sequence.name.strip():
        locus.append(str(sequence.name[:16]))
    else:
        locus.append(os.path.splitext(os.path.basename(handle.name))[0][:16])
    locus.append(' ' * (39 - (12 + min(len(sequence), 16))))
    locus.append(str(len(sequence)) + ' bp ')
    if sequence.ds:
        locus.append('ds-')
    else:
        locus.append('ss-')
    locus.append(sequence.material.upper() + '   ')
    if sequence.circular:
        locus.append('circular')
    else:
        locus.append('linear  ')
    locus.append('    ')
    today = datetime.date.today()
    locus.append(today.strftime('%d-') + today.strftime('%b').upper() +
                 today.strftime('-%Y'))
    lines.append(''.join(locus))

    # Construct the DEFINITION keyword, which is supposed to be used for the
    # organism. In our case, we'll just re-use the name
    lines.append('DEFINITION  ' + sequence.name)

    # Construct the ACCESSION keyword. This is intentionally left blank until
    # someone has a use case for this in synbio
    lines.append('ACCESSION   .')

    # Construct the KEYWORDS keyword - also intentionally left blank
    lines.append('KEYWORDS    .')

    # Construct the SOURCE keyword - again used for the organism. This is not
    # expected for the vast majority of designs, and is also left blank
    lines.append('SOURCE      .')

    # Construct the FEATURES keyword - this one is actually used.
    lines.append('FEATURES             Location/Qualifiers')

    def make_qualifier(key, value):
        indent = ' ' * 21
        wraplen = 79 - 21
        qualifier = '/{}="{}"'.format(key, value)
        return [indent + l for l in textwrap.wrap(qualifier, wraplen)]

    for feature in sequence.features:
        # Make the feature header (feature type + location)
        if feature.gaps:
            loc = 'join({}'.format(feature.start)
            for gap in feature.gaps:
                loc += '..' + str(gap[0] + 1) + ',' + str(gap[1])
            loc += '..{})'.format(feature.stop)
        else:
            if feature.start == feature.stop + 1:
                loc = str(feature.stop)
            else:
                loc = str(feature.start + 1) + '..' + str(feature.stop)
        if feature.strand == 1:
            loc = 'complement({})'.format(loc)
        lines.append(''.join(['     ', feature.feature_type,
                              ' ' * (16 - len(feature.feature_type)), loc]))

        # Add qualifiers
        lines += make_qualifier('label', feature.name)
        for key, value in feature.qualifiers.iteritems():
            lines += make_qualifier(key, value)

    # Create ORIGIN keyword - stores the sequence
    lines.append('ORIGIN')
    chunks = textwrap.wrap(str(sequence), 10)
    for i in range(int(math.ceil(len(chunks) / 6.))):
        index = str(i * 60 + 1)
        lines.append(' ' * (9 - len(index)) + index + ' ' +
                     ' '.join(chunks[i * 6:i * 6 + 6]))
    # End genbank file
    lines.append('//')

    handle.write('\n'.join(lines))
