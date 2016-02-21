'''Nearest-neighbor method Tm calculation parameters.

Nearest-neighbor parameters don't publish full NN parameters. Assumptions:
    AA = TT
    GG = CC
    CA = TG
    CT = AG
    GA = TC
    GT = AC

'''


BRESLAUER = {
    'delta_h': {
        'AA': 9.1,
        'TT': 9.1,
        'AT': 8.6,
        'TA': 6.0,
        'CA': 5.8,
        'TG': 5.8,
        'GT': 6.5,
        'AC': 6.5,
        'CT': 7.8,
        'AG': 7.8,
        'GA': 5.6,
        'TC': 5.6,
        'CG': 11.9,
        'GC': 11.1,
        'GG': 11.0,
        'CC': 11.0},
    'delta_h_err': {
        'anyGC': 0.0,
        'onlyAT': 0.0,
        'symmetry': 0.0,
        'terminalT': 0.0},
    'delta_s': {
        'AA': 24.0,
        'TT': 24.0,
        'AT': 23.9,
        'TA': 16.9,
        'CA': 12.9,
        'TG': 12.9,
        'GT': 17.3,
        'AC': 17.3,
        'CT': 20.8,
        'AG': 20.8,
        'GA': 13.5,
        'TC': 13.5,
        'CG': 27.8,
        'GC': 26.7,
        'GG': 26.6,
        'CC': 26.6},
    'delta_s_err': {
        'anyGC': 16.77,
        'onlyAT': 20.13,
        'symmetry': 1.34,
        'terminalT': 0.0}}


SANTALUCIA96 = {
    'delta_h': {
        'AA': 8.4,
        'TT': 8.4,
        'AT': 6.5,
        'TA': 6.3,
        'CA': 7.4,
        'TG': 7.4,
        'GT': 8.6,
        'AC': 8.6,
        'CT': 6.1,
        'AG': 6.1,
        'GA': 7.7,
        'TC': 7.7,
        'CG': 10.1,
        'GC': 11.1,
        'GG': 6.7,
        'CC': 6.7},
    'delta_h_err': {
        'anyGC': 0.0,
        'onlyAT': 0.0,
        'symmetry': 0.0,
        'terminalT': -0.4},
    'delta_s': {
        'AA': 23.6,
        'TT': 23.6,
        'AT': 18.8,
        'TA': 18.5,
        'CA': 19.3,
        'TG': 19.3,
        'GT': 23.0,
        'AC': 23.0,
        'CT': 16.1,
        'AG': 16.1,
        'GA': 20.3,
        'TC': 20.3,
        'CG': 25.5,
        'GC': 28.4,
        'GG': 15.6,
        'CC': 15.6},
    'delta_s_err': {
        'anyGC': 5.9,
        'onlyAT': 9.0,
        'symmetry': 1.4,
        'terminalT': 0.0}}


SUGIMOTO = {
    'delta_h': {
        'AA': 8.0,
        'TT': 8.0,
        'AT': 5.6,
        'TA': 6.6,
        'CA': 8.2,
        'TG': 8.2,
        'GT': 9.4,
        'AC': 9.4,
        'CT': 6.6,
        'AG': 6.6,
        'GA': 8.8,
        'TC': 8.8,
        'CG': 11.8,
        'GC': 10.5,
        'GG': 10.9,
        'CC': 10.9},
    'delta_h_err': {
        'anyGC': -0.6,
        'onlyAT': -0.6,
        'symmetry': 0.0,
        'terminalT': 0.0},
    'delta_s': {
        'AA': 21.9,
        'TT': 21.9,
        'AT': 15.2,
        'TA': 18.4,
        'CA': 21.0,
        'TG': 21.0,
        'GT': 25.5,
        'AC': 25.5,
        'CT': 16.4,
        'AG': 16.4,
        'GA': 23.5,
        'TC': 23.5,
        'CG': 29.0,
        'GC': 26.4,
        'GG': 28.4,
        'CC': 28.4},
    'delta_s_err': {
        'anyGC': 9.0,
        'onlyAT': 9.0,
        'symmetry': 1.4,
        'terminalT': 0.0}}


SANTALUCIA98 = {
    'delta_h': {
        'AA': 7.9,
        'TT': 7.9,
        'AT': 7.2,
        'TA': 7.2,
        'CA': 8.5,
        'TG': 8.5,
        'GT': 8.4,
        'AC': 8.4,
        'CT': 7.8,
        'AG': 7.8,
        'GA': 8.2,
        'TC': 8.2,
        'CG': 10.6,
        'GC': 9.8,
        'GG': 8.0,
        'CC': 8.0},
    'delta_h_err': {
        'initGC': -0.1,
        'initAT': -2.3,
        'symmetry': 0.0},
    'delta_s': {
        'AA': 22.2,
        'TT': 22.2,
        'AT': 20.4,
        'TA': 21.3,
        'CA': 22.7,
        'TG': 22.7,
        'GT': 22.4,
        'AC': 22.4,
        'CT': 21.0,
        'AG': 21.0,
        'GA': 22.2,
        'TC': 22.2,
        'CG': 27.2,
        'GC': 24.4,
        'GG': 19.9,
        'CC': 19.9},
    'delta_s_err': {
        'initGC': 2.8,
        'initAT': -4.1,
        'symmetry': 1.4}}


CLONING = {
    'delta_h': {
        'AA': 9.1,
        'TT': 9.1,
        'AT': 8.6,
        'TA': 6.0,
        'CA': 5.8,
        'TG': 5.8,
        'GT': 6.5,
        'AC': 6.5,
        'CT': 7.8,
        'AG': 7.8,
        'GA': 5.6,
        'TC': 5.6,
        'CG': 11.9,
        'GC': 11.1,
        'GG': 11.0,
        'CC': 11.0},
    'delta_h_err': {
        'anyGC': 0.0,
        'onlyAT': 0.0,
        'symmetry': 0.0,
        'terminalT': 0.0},
    'delta_s': {
        'AA': 24.0,
        'TT': 24.0,
        'AT': 23.9,
        'TA': 16.9,
        'CA': 12.9,
        'TG': 12.9,
        'GT': 17.3,
        'AC': 17.3,
        'CT': 20.8,
        'AG': 20.8,
        'GA': 13.5,
        'TC': 13.5,
        'CG': 27.8,
        'GC': 26.7,
        'GG': 26.6,
        'CC': 26.6},
    'delta_s_err': {
        'onlyAT': 0.0,
        'anyGC': 0.0,
        'symmetry': 0.0,
        'terminalT': 0.0}}
