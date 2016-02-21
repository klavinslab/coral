'''Check for primer dimers using NUPACK.'''
import coral as cr


def dimers(primer1, primer2, concentrations=[5e-7, 3e-11]):
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
    :returns: Fraction of dimers versus the total amount of primer added.
    :rtype: float

    '''
    # It is not reasonable (yet) to use a long template for doing these
    # computations directly, as NUPACK does an exhaustive calculation and
    # would take too long without a cluster.
    # Instead, this function compares primer-primer binding to
    # primer-complement binding

    # Simulate binding of template vs. primers
    nupack = cr.structure.NUPACK([primer1.primer(), primer2.primer(),
                                  primer1.primer().reverse_complement(),
                                  primer2.primer().reverse_complement()])
    # Include reverse complement concentration
    primer_concs = [concentrations[0]] * 2
    template_concs = [concentrations[1]] * 2
    concs = primer_concs + template_concs
    nupack_concs = nupack.concentrations(2, conc=concs)
    dimer_conc = nupack_concs['concentrations'][5]

    return dimer_conc / concs[0]
