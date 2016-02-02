'''Read and write DNA sequences.'''
import coral as cr


def parse_fasta(lines):
    # FASTA files are repeats of two elements: a description (starting with >)
    # and a potentially multi-line sequence.

    # Consume the lines until done
    seqs = []
    seq = cr.DNA('')
    description = lines.pop(0)
    for line in lines:
        if line.startswith('>'):
            seq.description = description
            seqs.append(seq)
            description = line
        else:
            seq += cr.DNA(line)

    return seqs
