'''Read and write DNA sequences.'''
import csv
import os
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.SeqFeature import CompoundLocation
import pymbt
import pymbt.constants.genbank


class PrimerAnnotationError(ValueError):
    pass


class FeatureNameError(ValueError):
    pass


def read_dna(path):
    '''Read DNA from file. Uses BioPython and coerces to pymbt format.

    :param path: Full path to input file.
    :type path: str
    :returns: DNA sequence.
    :rtype: pymbt.DNA

    '''
    filename, ext = os.path.splitext(os.path.split(path)[-1])

    genbank_exts = ['.gb', '.ape']
    fasta_exts = ['.fasta', '.fa', '.fsa', '.seq']
    abi_exts = ['.abi', '.ab1']

    if any([ext == extension for extension in genbank_exts]):
        file_format = 'genbank'
    elif any([ext == extension for extension in fasta_exts]):
        file_format = 'fasta'
    elif any([ext == extension for extension in abi_exts]):
        file_format = 'abi'
    else:
        raise ValueError('File format not recognized.')

    seq = SeqIO.read(path, file_format)
    dna = pymbt.DNA(str(seq.seq))
    if seq.name == '.':
        dna.name = filename
    else:
        dna.name = seq.name

    # Features
    for feature in seq.features:
        try:
            dna.features.append(_seqfeature_to_pymbt(feature))
        except FeatureNameError:
            pass
    dna.features = sorted(dna.features, key=lambda feature: feature.start)
    # Used to use data_file_division, but it's inconsistent (not always the
    # molecule type)
    dna.topology = 'linear'
    with open(path) as f:
        first_line = f.read().split()
        for word in first_line:
            if word == 'circular':
                dna.topology = 'circular'

    return dna


def read_sequencing(directory):
    '''Read .seq and .abi/.ab1 results files from a dir.

    :param directory: Path to directory containing sequencing files.
    :type directory: str
    :returns: A list of DNA sequences.
    :rtype: pymbt.DNA list

    '''
    dirfiles = os.listdir(directory)
    seq_exts = ['.seq', '.abi', '.ab1']
    # Exclude files that aren't sequencing results
    seq_paths = [x for x in dirfiles if os.path.splitext(x)[1] in seq_exts]
    paths = [os.path.join(directory, x) for x in seq_paths]
    sequences = [read_dna(x).to_ss() for x in paths]

    return sequences


def write_dna(dna, path):
    '''Write DNA to a file (genbank or fasta).

    :param dna: DNA sequence to write to file
    :type dna: pymbt.DNA
    :param path: file path to write. Has to be genbank or fasta file.
    :type path: str

    '''
    # Check if path filetype is valid, remember for later
    ext = os.path.splitext(path)[1]
    if ext == '.gb' or ext == '.ape':
        filetype = 'genbank'
    elif ext == '.fa' or ext == '.fasta':
        filetype = 'fasta'
    else:
        raise ValueError('Only genbank or fasta files are supported.')

    # Convert features to Biopython form
    # Information lost on conversion:
    #     specificity of feature type
    #     strandedness
    #     topology
    features = []
    for feature in dna.features:
        features.append(_pymbt_to_seqfeature(feature))
    # Biopython doesn't like 'None' here
    bio_id = dna.id if dna.id else ''
    # Maximum length of name is 16
    seq = SeqRecord(Seq(str(dna), alphabet=ambiguous_dna), id=bio_id,
                    name=dna.name[0:16], features=features,
                    description=dna.name)
    seq.annotations['data_file_division'] = dna.topology

    if filetype == 'genbank':
        SeqIO.write(seq, path, 'genbank')
    elif filetype == 'fasta':
        SeqIO.write(seq, path, 'fasta')


def write_primers(primer_list, path, names=None, notes=None):
    '''Write a list of primers out to a csv file. The first three columns are
    compatible with the current IDT order form (name, sequence, notes). By
    default there are no notes, which is an optional parameter.

    :param primer_list: A list of primers.
    :type primer_list: pymbt.Primer list
    :param path: A path to the csv you want to write.
    :type path: str
    :param names: A list of strings to name each oligo. Must be the same length
                  as the primer_list.
    :type names: str list
    :param notes: A list of strings to provide a note for each oligo. Must be
                  the same length as the primer_list.
    :type notes: str list

    '''
    # Check for notes and names having the right length, apply them to primers
    if names is not None:
        if len(names) != len(primer_list):
            names_msg = 'Mismatch in number of notes and primers.'
            raise PrimerAnnotationError(names_msg)
        for i, name in enumerate(names):
            primer_list[i].name = name
    if notes is not None:
        if len(notes) != len(primer_list):
            notes_msg = 'Mismatch in number of notes and primers.'
            raise PrimerAnnotationError(notes_msg)
        for i, note in enumerate(notes):
            primer_list[i].note = note

    # Write to csv
    with open(path, 'w') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['name', 'sequence', 'notes'])
        for primer in primer_list:
            string_rep = str(primer.overhang).lower() + str(primer.anneal)
            writer.writerow([primer.name, string_rep, primer.note])


def _process_feature_type(feature_type, bio_to_pymbt=True):
    '''Translate genbank feature types into usable ones (currently identical).
    The feature table is derived from the official genbank spec (gbrel.txt)
    available at http://www.insdc.org/documents/feature-table

    :param feature_type: feature to convert
    :type feature_type: str
    :param bio_to_pymbt: from pymbt to Biopython (True) or the other direction
                   (False)
    :param bio_to_pymbt: bool
    :returns: pymbt version of genbank feature_type, or vice-versa.
    :rtype: str

    '''

    err_msg = 'Unrecognized feature type: {}'.format(feature_type)
    if bio_to_pymbt:
        try:
            name = pymbt.constants.genbank.TO_PYMBT[feature_type]
        except KeyError:
            raise ValueError(err_msg)
    else:
        try:
            name = pymbt.constants.genbank.TO_BIO[feature_type]
        except KeyError:
            raise ValueError(err_msg)
    return name


def _seqfeature_to_pymbt(feature):
    '''Convert a Biopython SeqFeature to a pymbt.Feature.

    :param feature: Biopython SeqFeature
    :type feature: Bio.SeqFeature

    '''
    # Some genomic sequences don't have a label attribute
    # TODO: handle genomic cases differently than others. Some features lack
    # a label but should still be incorporated somehow.
    qualifiers = feature.qualifiers
    if 'label' in qualifiers:
        feature_name = qualifiers['label'][0]
    elif 'locus_tag' in qualifiers:
        feature_name = qualifiers['locus_tag'][0]
    else:
        raise FeatureNameError('Unrecognized feature name')
    # Features with gaps are special, require looking at subfeatures
    # Assumption: subfeatures are never more than one level deep
    if feature.location_operator == 'join':
        # Feature has gaps. Have to figure out start/stop from subfeatures,
        # calculate gap indices. A nested feature model may be required
        # eventually.
        # Reorder the sub_feature list by start location
        # Assumption: none of the subfeatures overlap so the last entry in
        # the reordered list also has the final stop point of the feature.
        # FIXME: Getting a deprecation warning about using sub_features
        # instead of feature.location being a CompoundFeatureLocation
        reordered = sorted(feature.location.parts,
                           key=lambda location: location.start)
        starts = [int(location.start) for location in reordered]
        stops = [int(location.end) for location in reordered]
        feature_start = starts.pop(0)
        feature_stop = stops.pop(-1)
        starts = [start - feature_start for start in starts]
        stops = [stop - feature_start for stop in stops]
        feature_gaps = list(zip(stops, starts))
    else:
        # Feature doesn't have gaps. Ignore subfeatures.
        feature_start = int(feature.location.start)
        feature_stop = int(feature.location.end)
        feature_gaps = []
    feature_type = _process_feature_type(feature.type)
    if feature.location.strand == -1:
        feature_strand = 1
    else:
        feature_strand = 0
    if 'gene' in qualifiers:
        gene = qualifiers['gene']
    else:
        gene = []
    if 'locus_tag' in qualifiers:
        locus_tag = qualifiers['locus_tag']
    else:
        locus_tag = []
    pymbt_feature = pymbt.Feature(feature_name, feature_start,
                                  feature_stop, feature_type,
                                  gene=gene, locus_tag=locus_tag,
                                  qualifiers=qualifiers,
                                  strand=feature_strand,
                                  gaps=feature_gaps)
    return pymbt_feature


def _pymbt_to_seqfeature(feature):
    '''Convert a pymbt.Feature to a Biopython SeqFeature.

    :param feature: pymbt Feature.
    :type feature: pymbt.Feature

    '''
    bio_strand = 1 if feature.strand == 1 else -1
    ftype = _process_feature_type(feature.feature_type, bio_to_pymbt=False)
    sublocations = []
    if feature.gaps:
        # There are gaps. Have to define location_operator and  add subfeatures
        location_operator = 'join'
        # Feature location means nothing for 'join' sequences?
        # TODO: verify
        location = FeatureLocation(ExactPosition(0), ExactPosition(1),
                                   strand=bio_strand)
        # Reconstruct start/stop indices for each subfeature
        stops, starts = zip(*feature.gaps)
        starts = [feature.start] + [start + feature.start for start in starts]
        stops = [stop + feature.start for stop in stops] + [feature.stop]
        # Build subfeatures
        for start, stop in zip(starts, stops):
            sublocation = FeatureLocation(ExactPosition(start),
                                          ExactPosition(stop),
                                          strand=bio_strand)
            sublocations.append(sublocation)
        location = CompoundLocation(sublocations, operator='join')
    else:
        # No gaps, feature is simple
        location_operator = ''
        location = FeatureLocation(ExactPosition(feature.start),
                                   ExactPosition(feature.stop),
                                   strand=bio_strand)
    seqfeature = SeqFeature(location, type=ftype,
                            qualifiers={'label': [feature.name]},
                            location_operator=location_operator)
    return seqfeature
