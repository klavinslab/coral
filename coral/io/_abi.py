'''Read and write DNA sequences.'''
import coral as cr
import numpy as np
import os
from . import parsers
from .exceptions import UnsupportedFileError


def read_abi(path, trim=True, attach_trace=True):
    '''Read a single ABI/AB1 Sanger sequencing file.

    :param path: Full path to input file.
    :type path: str
    :param trim: Determines whether the sequence will be trimmed using Richard
                 Mott's algorithm (trims based on quality).
    :type trim: bool
    :param attach_trace: Determines whether to attach the trace result as a
                         .trace attribute of the returned sequence and the
                         trace peak locations as a .tracepeaks attribute. The
                         trace attribute is a 2D numpy array with 4 columns in
                         the order GATC.
    :type attach_trace: bool
    :returns: DNA sequence.
    :rtype: coral.DNA

    '''
    filename, ext = os.path.splitext(os.path.split(path)[-1])

    abi_exts = ['.abi', '.ab1']

    if ext in abi_exts:
        with open(path) as f:
            abi = parsers.ABI(f)
    else:
        raise UnsupportedFileError('File format not recognized.')

    seq = abi.seq_remove_ambig(abi.seq)

    # Attach the trace results to the seq
    if attach_trace:
        order = abi.data['baseorder'].upper()
        trace = [abi.data['raw' + str(order.index(b) + 1)] for b in 'GATC']
        trace = np.array(trace)
        tracepeaks = np.array(abi.data['tracepeaks'])

    if trim:
        try:
            sequence = cr.DNA(abi.trim(seq))
        except ValueError:
            # A ValueError is raised if the sequence is too short
            pass

        trim_start = seq.index(str(sequence))
        # Adjust trace data based on trimming
        idx = (trim_start, trim_start + len(sequence))
        peaks = tracepeaks[idx[0]:idx[1]]
        sequence.trace = trace[peaks[0]:peaks[-1], :]
        sequence.tracepeaks = peaks
    else:
        sequence = cr.DNA(seq)

    sequence.name = abi.name

    return sequence


def read_abis(directory, trim=True, attach_trace=True):
    '''Read all ABI sequences files in a directory.

    :param directory: Path to directory containing sequencing files.
    :type directory: str
    :param trim: Determines whether the sequence will be trimmed using Richard
                 Mott's algorithm (trims based on quality).
    :type trim: bool
    :param attach_trace: Determines whether to attach the trace result as a
                         .trace attribute of the returned sequence. The trace
                         attribute is a 2D numpy array with 4 columns in the
                         order GATC.
    :type attach_trace: bool
    :returns: A list of DNA sequences.
    :rtype: coral.DNA list

    '''
    dirfiles = os.listdir(directory)
    abis = []
    for dirfile in dirfiles:
        path = os.path.join(directory, dirfile)
        try:
            abis.append(read_abi(path, trim=trim, attach_trace=attach_trace))
        except UnsupportedFileError:
            pass

    return abis
