'''Check for primer dimers using NUPACK.'''
import coral as cr


def dimers(primer1, primer2, concentrations=[5e-7, 3e-11], temp=None):
    '''Calculate expected fraction of primer dimers.

    :param primer1: Forward primer.
    :type primer1: coral.DNA
    :param primer2: Reverse primer.
    :type primer2: coral.DNA
    :param template: DNA template.
    :type template: coral.DNA
    :param concentrations: list of concentrations for primers and the
                           template. Defaults are those for PCR with 1kb
                           template.
    :type concentrations: list
    :param temp: Temperature at which to do the simulation (e.g. the Tm) in C.
    :type temp: float
    :returns: Fraction of dimers versus the total amount of primer added.
    :rtype: float

    '''
    # It is not reasonable (yet) to use a long template for doing these
    # computations directly, as NUPACK does an exhaustive calculation and
    # would take too long without a cluster.
    # Instead, this function compares primer-primer binding to
    # primer-complement binding

    # If temp is not supplied, use primer Tm average
    if temp is None:
        temp = (primer1.tm + primer2.tm) / 2.

    # Simulate binding of template vs. primers
    nupack = cr.structure.NUPACK()
    # Need output of Nupack's complexes as input to 'concentrations'
    primer1seq = primer1.primer()
    primer2seq = primer2.primer()
    species = [primer1seq, primer2seq, primer1seq.reverse_complement(),
               primer2seq.reverse_complement()]

    complexes = nupack.complexes(species, 2, temp=temp)

    # Include reverse complement concentration
    primer_concs = [concentrations[0]] * 2
    template_concs = [concentrations[1]] * 2
    concs = primer_concs + template_concs
    nupack_concs = nupack.concentrations(complexes, concs)
    dimer_conc = nupack_concs[5]['concentration']
    print species

    return dimer_conc / concs[0]
