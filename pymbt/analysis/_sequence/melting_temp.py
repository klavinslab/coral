# -*- coding: utf-8
'''Calculate the thermodynamic melting temperatures of nucleotide sequences.'''
from math import log, log10
from . import tm_params

# TODO: Owczarzy et al 2004 has better salt correction
# TODO: Remove sugimoto? It's missing important details (like salt correction)
# TODO: Write adjusted santalucia method. Did a fit of 20000 sequences
# comparing santalucia98 to finnzymes' modified breslauer method. As expected
# they correlate heavily with offset of -6.014 and slope of 1.335 -
# i.e. to get approximately the same result, do santalucia98, multiply by 1.335
# , and subtract 6.
# Double check those stats
# TODO: Make new hybrid method - combine santalucia unified with owczarzy
# corrections, compare to finnzymes.

# Methods that are currently verified to work using reference sequences:
#   'cloning'
#   'santalucia98'


def tm(seq, dna_conc=50, salt_conc=50, parameters='cloning'):
    '''Calculate nearest-neighbor melting temperature (Tm).

    :param seq: Sequence for which to calculate the tm.
    :type seq: pymbt.DNA
    :param dna_conc: DNA concentration in nM.
    :type dna_conc: float
    :param salt_conc: Salt concentration in mM.
    :type salt_conc: float
    :param parameters: Nearest-neighbor parameter set. Available options:
                       'breslauer': Breslauer86 parameters
                       'sugimoto': Sugimoto96 parameters
                       'santalucia96': SantaLucia96 parameters
                       'santalucia98': SantaLucia98 parameters
                       'cloning': breslauer without corrections
                       'cloning_sl98': santalucia98 fit to 'cloning'
    :type parameters: str
    :returns: Melting temperature (Tm) in °C.
    :rtype: float
    :raises: ValueError if parameter argument is invalid.

    '''
    if parameters == 'breslauer':
        params = tm_params.BRESLAUER
    elif parameters == 'sugimoto':
        params = tm_params.SUGIMOTO
    elif parameters == 'santalucia96':
        params = tm_params.SANTALUCIA96
    elif parameters == 'santalucia98' or parameters == 'cloning_sl98':
        params = tm_params.SANTALUCIA98
    elif parameters == 'cloning':
        params = tm_params.CLONING
    else:
        raise ValueError('Unsupported parameter set.')

    # Thermodynamic parameters
    pars = {'delta_h': params['delta_h'], 'delta_s': params['delta_s']}
    pars_error = {'delta_h': params['delta_h_err'],
                  'delta_s': params['delta_s_err']}

    # Error corrections - done first for use of reverse_complement parameters
    if parameters == 'breslauer':
        deltas = breslauer_corrections(seq, pars_error)
    elif parameters == 'sugimoto':
        deltas = breslauer_corrections(seq, pars_error)
    elif parameters == 'santalucia96':
        deltas = breslauer_corrections(seq, pars_error)
    elif parameters == 'santalucia98' or parameters == 'cloning_sl98':
        deltas = santalucia98_corrections(seq, pars_error)
    elif parameters == 'cloning':
        deltas = breslauer_corrections(seq, pars_error)
        deltas[0] += 3.4
        deltas[1] += 12.4

    # Sum up the nearest-neighbor enthalpy and entropy
    seq = str(seq).upper()

    new_delt = _pair_deltas(seq, pars)
    deltas[0] += new_delt[0]
    deltas[1] += new_delt[1]

    # Unit corrections
    salt_conc /= 1e3
    dna_conc /= 1e9
    deltas[0] *= 1e3

    # Universal gas constant (R)
    R = 1.9872

    # Supposedly this is what dnamate does, but the output doesn't match theirs
#    melt = (-deltas[0] / (-deltas[1] + R * log(dna_conc / 4.0))) +
#                          16.6 * log(salt_conc) - 273.15
#    return melt
    # Overall equation is supposedly:
    # sum{dH}/(sum{dS} + R ln(dna_conc/b)) - 273.15
    # with salt corrections for the whole term (or for santalucia98,
    # salt corrections added to the dS term.
    # So far, implementing this as described does not give results that match
    # any calculator but Biopython's

    if parameters == 'breslauer' or parameters == 'cloning':
        numerator = -deltas[0]
        # Modified dna_conc denominator
        denominator = (-deltas[1]) + R * log(dna_conc / 16.0)
        # Modified Schildkraut-Lifson equation adjustment
        salt_adjustment = 16.6 * log(salt_conc) / log(10.0)
        melt = numerator / denominator + salt_adjustment - 273.15
    elif parameters == 'santalucia98' or 'cloning_sl98':
        # TODO: dna_conc should be divided by 2.0 when dna_conc >> template
        # (like PCR)
        numerator = -deltas[0]
        # SantaLucia 98 salt correction
        salt_adjustment = 0.368 * (len(seq) - 1) * log(salt_conc)
        denominator = -deltas[1] + salt_adjustment + R * log(dna_conc / 4.0)
        melt = -deltas[0] / denominator - 273.15
    elif parameters == 'santalucia96':
        # TODO: find a way to test whether the code below matches another
        # algorithm. It appears to be correct, but need to test it.
        numerator = -deltas[0]
        denominator = -deltas[1] + R * log(dna_conc / 4.0)
        # SantaLucia 96 salt correction
        salt_adjustment = 12.5 * log10(salt_conc)
        melt = numerator / denominator + salt_adjustment - 273.15
    elif parameters == 'sugimoto':
        # TODO: the stuff below is untested and probably wrong
        numerator = -deltas[0]
        denominator = -deltas[1] + R * log(dna_conc / 4.0)
        # Sugimoto parameters were fit holding salt concentration constant
        # Salt correction can be chosen / ignored? Remove sugimoto set since
        # it's so similar to santalucia98?
        salt_correction = 16.6 * log10(salt_conc)
        melt = numerator / denominator + salt_correction - 273.15

    if parameters == 'cloning_sl98':
        # Corrections to make santalucia98 method approximate cloning method.
        # May be even better for cloning with Phusion than 'cloning' method
        melt *= 1.27329212575
        melt += -2.55585450119

    return melt


def _pair_deltas(seq, pars):
    '''Add up nearest-neighbor parameters for a given sequence.

    :param seq: DNA sequence for which to sum nearest neighbors
    :type seq: str
    :param pars: parameter set to use
    :type pars: dict
    :returns: nearest-neighbor delta_H and delta_S sums.
    :rtype: tuple of floats

    '''
    delta0 = 0
    delta1 = 0
    for i in range(len(seq) - 1):
        curchar = seq[i:i + 2]
        delta0 += pars['delta_h'][curchar]
        delta1 += pars['delta_s'][curchar]
    return delta0, delta1


def breslauer_corrections(seq, pars_error):
    '''Sum corrections for Breslauer '84 method.

    :param seq: sequence for which to calculate corrections.
    :type seq: str
    :param pars_error: dictionary of error corrections
    :type pars_error: dict
    :returns: Corrected delta_H and delta_S parameters
    :rtype: list of floats

    '''
    deltas_corr = [0, 0]
    contains_gc = 'G' in str(seq) or 'C' in str(seq)
    only_at = str(seq).count('A') + str(seq).count('T') == len(seq)
    symmetric = seq == seq.reverse_complement()
    terminal_t = str(seq)[0] == 'T' + str(seq)[-1] == 'T'

    for i, delta in enumerate(['delta_h', 'delta_s']):
        if contains_gc:
            deltas_corr[i] += pars_error[delta]['anyGC']
        if only_at:
            deltas_corr[i] += pars_error[delta]['onlyAT']
        if symmetric:
            deltas_corr[i] += pars_error[delta]['symmetry']
        if terminal_t and delta == 'delta_h':
            deltas_corr[i] += pars_error[delta]['terminalT'] * terminal_t
    return deltas_corr


def santalucia98_corrections(seq, pars_error):
    '''Sum corrections for SantaLucia '98 method (unified parameters).

    :param seq: sequence for which to calculate corrections.
    :type seq: str
    :param pars_error: dictionary of error corrections
    :type pars_error: dict
    :returns: Corrected delta_H and delta_S parameters
    :rtype: list of floats

    '''
    deltas_corr = [0, 0]
    first = str(seq)[0]
    last = str(seq)[-1]

    start_gc = first == 'G' or first == 'C'
    start_at = first == 'A' or first == 'T'
    end_gc = last == 'G' or last == 'C'
    end_at = last == 'A' or last == 'T'
    init_gc = start_gc + end_gc
    init_at = start_at + end_at

    symmetric = seq == seq.reverse_complement()

    for i, delta in enumerate(['delta_h', 'delta_s']):
        deltas_corr[i] += init_gc * pars_error[delta]['initGC']
        deltas_corr[i] += init_at * pars_error[delta]['initAT']
        if symmetric:
            deltas_corr[i] += pars_error[delta]['symmetry']
    return deltas_corr
