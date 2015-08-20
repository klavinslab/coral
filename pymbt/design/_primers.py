'''Primer design tools.'''
import pymbt
import warnings


def primer(dna, tm=65, min_len=10, tm_undershoot=1, tm_overshoot=3,
           end_gc=False, tm_parameters='cloning', overhang=None,
           structure=False):
    '''Design primer to a nearest-neighbor Tm setpoint.

    :param dna: Sequence for which to design a primer.
    :type dna: pymbt.DNA
    :param tm: Ideal primer Tm in degrees C.
    :type tm: float
    :param min_len: Minimum primer length.
    :type min_len: int
    :param tm_undershoot: Allowed Tm undershoot.
    :type tm_undershoot: float
    :param tm_overshoot: Allowed Tm overshoot.
    :type tm_overshoot: float
    :param end_gc: Obey the 'end on G or C' rule.
    :type end_gc: bool
    :param tm_parameters: Melting temp calculator method to use.
    :type tm_parameters: string
    :param overhang: Append the primer to this overhang sequence.
    :type overhang: str
    :param structure: Evaluate primer for structure, with warning for high
                      structure.
    :type structure: bool
    :returns: A primer.
    :rtype: pymbt.Primer
    :raises: ValueError if the input sequence is lower than the Tm settings
             allow.
             ValueError if a primer ending with G or C can't be found given
             the Tm settings.

    '''
    # Check Tm of input sequence to see if it's already too low
    seq_tm = pymbt.analysis.tm(dna, parameters=tm_parameters)
    if seq_tm < (tm - tm_undershoot):
        msg = 'Input sequence Tm is lower than primer Tm setting'
        raise ValueError(msg)
    # Focus on first 90 bases - shouldn't need more than 90bp to anneal
    dna = dna[0:90]

    # Generate primers from min_len to 'tm' + tm_overshoot
    # TODO: this is a good place for optimization. Only calculate as many
    # primers as are needed. Use binary search.
    primers_tms = []
    last_tm = 0
    bases = min_len
    while last_tm <= tm + tm_overshoot and bases != len(dna):
        next_primer = dna[0:bases]
        last_tm = pymbt.analysis.tm(next_primer, parameters=tm_parameters)
        primers_tms.append((next_primer, last_tm))
        bases += 1

    # Trim primer list based on tm_undershoot and end_gc
    primers_tms = [(primer, melt) for primer, melt in primers_tms if
                   melt >= tm - tm_undershoot]
    if end_gc:
        primers_tms = [(primer, melt) for primer, melt in primers_tms if
                       primer.top().endswith(('C', 'G'))]
    if not primers_tms:
        raise ValueError('No primers could be generated using these settings')

    # Find the primer closest to the set Tm, make it single stranded
    tm_diffs = [abs(melt - tm) for primer, melt in primers_tms]
    best_index = tm_diffs.index(min(tm_diffs))
    best_primer, best_tm = primers_tms[best_index]
    best_primer = best_primer.to_ss()

    # Apply overhang
    if overhang:
        overhang = overhang.to_ss()

    output_primer = pymbt.Primer(best_primer, best_tm, overhang=overhang)

    def _structure(primer):
        '''Check annealing sequence for structure.

        :param primer: Primer for which to evaluate structure
        :type primer: sequence.Primer

        '''
        # Check whole primer for high-probability structure, focus in on
        # annealing sequence, report average
        nupack = pymbt.analysis.Nupack(primer.primer())
        pairs = nupack.pairs(0)
        anneal_len = len(primer.anneal)
        pairs_mean = sum(pairs[-anneal_len:]) / anneal_len
        if pairs_mean < 0.5:
            warnings.warn('High probability structure', Warning)
        return pairs_mean
    if structure:
        _structure(output_primer)
    return output_primer


def primers(dna, tm=65, min_len=10, tm_undershoot=1, tm_overshoot=3,
            end_gc=False, tm_parameters='cloning', overhangs=None,
            structure=False):
    '''Design primers for PCR amplifying any arbitrary sequence.

    :param dna: Input sequence.
    :type dna: pymbt.DNA
    :param tm: Ideal primer Tm in degrees C.
    :type tm: float
    :param min_len: Minimum primer length.
    :type min_len: int
    :param tm_undershoot: Allowed Tm undershoot.
    :type tm_undershoot: float
    :param tm_overshoot: Allowed Tm overshoot.
    :type tm_overshoot: float
    :param end_gc: Obey the 'end on G or C' rule.
    :type end_gc: bool
    :param tm_parameters: Melting temp calculator method to use.
    :type tm_parameters: string
    :param overhangs: 2-tuple of overhang sequences.
    :type overhangs: tuple
    :param structure: Evaluate each primer for structure, with warning for high
                      structure.
    :type structure: bool
    :returns: A list primers (the output of primer).
    :rtype: list

    '''
    if not overhangs:
        overhangs = [None, None]
    templates = [dna, dna.reverse_complement()]
    primer_list = []
    for template, overhang in zip(templates, overhangs):
        primer_i = primer(template, tm=tm, min_len=min_len,
                          tm_undershoot=tm_undershoot,
                          tm_overshoot=tm_overshoot, end_gc=end_gc,
                          tm_parameters=tm_parameters,
                          overhang=overhang, structure=structure)
        primer_list.append(primer_i)
    return primer_list
