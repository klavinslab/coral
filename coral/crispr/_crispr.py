'''CRISPR-Cas design methods. The current goal is to capture the design methods
listed in the AddGene guide: http://www.addgene.org/CRISPR/guide/'''
import coral


def TemplateBoundsError(Exception):
    '''Raise if designed DNA would extend outside the template.'''


def locate_pams(sequence, pam='NGG'):
    '''Locate PAM sequences within a target sequence.

    :param sequence: Sequence to search for PAMs.
    :type sequence: coral.DNA
    :param pam: PAM sequence (
    :type pam: str or coral.DNA
    '''
    return coral.DNA(str(sequence)).locate(pam)


def grna_targets(sequence, length, min_tm=65, pam='NGG'):
    # Handle pam inputs
    pam = coral.DNA(str(pam))
    # Find all binding sites for the PAM
    pam_locations = locate_pams(sequence)

    # Design upstream target sequences
    def design_target(template, location, length):
        end = location
        start = location - length
        if start < 0:
            # length can't be satisfied
            raise TemplateBoundsError('gRNA extends outside of template.')
        target = template[start:end]
        return target

    targets = []
    grnas = []
    for i, strand in enumerate(pam_locations):
        for location in strand:
            if i == 0:
                target = design_target(sequence, location, length)
            else:
                target = design_target(sequence.reverse_complement(), location,
                                       length)
            targets.append(target)
            grna = target + pam
            grnas.append(grna)

    return grnas


def offtarget_check(grna, method='zhang'):
    pass
