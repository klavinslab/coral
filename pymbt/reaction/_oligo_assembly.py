"""Simulate building a construct by assembling oligos with PCA."""
# FIXME: Would not catch the case where e.g. the first and second oligos
# bound each other almost perfectly, ruling out the third oligo from binding
# i.e. the assemble_oligos function does not test for conflicting overlaps


class AssemblyError(Exception):
    """Raise exception if assembly can't complete for any reason."""
    pass


def assemble_oligos(dna_list, reference=None):
    """Given a list of DNA sequences, assemble into a single construct.
    :param dna_list: List of DNA sequences - they must be single-stranded.
    :type dna_list: pymbt.DNA list
    :param reference: Expected sequence - once assembly completed, this will
                      be used to reorient the DNA (assembly could potentially
                      occur from either side of a linear DNA construct if
                      oligos are in a random order). If this fails, an
                      AssemblyError is raised.
    :type reference: pymbt.DNA
    :raises: AssemblyError if it can't assemble for any reason.
    :returns: A single assembled DNA sequence
    :rtype: pymbt.DNA

    """
    # FIXME: this protocol currently only supports 5' ends on the assembly
    # Find all matches for every oligo. If more than 2 per side, error.
    # Self-oligo is included in case the 3' end is self-complementary.
    # 1) Find all unique 3' binders (and non-binders).
    match_3 = [bind_unique(seq, dna_list, right=True) for i, seq in
               enumerate(dna_list)]
    # 2) Find all unique 5' binders (and non-binders).
    match_5 = [bind_unique(seq, dna_list, right=False) for i, seq in
               enumerate(dna_list)]
    # Assemble into 2-tuple
    zipped = zip(match_5, match_3)
    #return zipped
    # 3) If none found, error out with 'oligo n has no binders'
    for i, oligo_match in enumerate(zipped):
        if not any(oligo_match):
            error = "Oligo {} has no binding partners.".format(i + 1)
            raise AssemblyError(error)
    # 4) There should be exactly 2 oligos that bind at 3' end but
    # not 5'.
    ends = []
    for i, (five, three) in enumerate(zipped):
        if five is None and three is not None:
            ends.append(i)
    # 5) If more than 2, error with 'too many ends'.
    if len(ends) > 2:
        raise AssemblyError("Too many (>2) end oligos found.")
    # 6) If more than 2, error with 'not enough ends'.
    if len(ends) < 2:
        raise AssemblyError("Not enough (<2) end oligos found.")
    # NOTE:If 1-4 are satisfied, unique linear assembly has been found (proof?)
    # 8) Start with first end and build iteratively
    last_index = ends[0]
    assembly = dna_list[last_index].to_ds()
    flip = True
    # This would be slightly less complicated if the sequences were tied to
    # their match info in a tuple
    # Append next region n - 1 times
    for i in range(len(dna_list) - 1):
        if flip:
            # Next oligo needs to be flipped before concatenation
            # Grab 3' match from last oligo's info
            current_index, matchlen = zipped[last_index][1]
            # Get new oligo sequence, make double-stranded for concatenation
            next_oligo = dna_list[current_index].to_ds()
            # Reverse complement for concatenation
            next_oligo = next_oligo.reverse_complement()
            # Don't reverse complement the next one
            flip = False
        else:
            # Grab 5' match from last oligo's info
            current_index, matchlen = zipped[last_index][0]
            # Get new oligo sequence, make double-stranded for concatenation
            next_oligo = dna_list[current_index].to_ds()
            # Reverse complement the next one
            flip = True
        # Trim overlap from new sequence
        next_oligo = next_oligo[(matchlen - 1):]
        # Concatenate and update last oligo's information
        assembly += next_oligo
        last_index = current_index
    if reference:
        if assembly == reference or assembly == reference.reverse_complement():
            return assembly
        else:
            raise AssemblyError("Assembly did not match reference")
    else:
        return assembly


def bind_unique(reference, query_list, min_overlap=12, right=True):
    """(5' or 3' region on reference sequence that uniquely matches the reverse
    complement of the associated (5' or 3') region of one sequence in a list of
    query sequences.

    :param reference: Reference sequence.
    :type reference: pymbt.DNA
    :param query_list: List of query sequences.
    :type query_list: pymbt.DNA list
    :param min_overlap: Minimum overlap for a match (in bp).
    :type min_overlap: int
    :param right: Check right side of sequence (3'). False results in 5' check.
    :type right: bool
    :returns: Tuple of the indices of any matches and the size of the match in
              bp.
    :rtype: tuple of ints
    :raises: AssemblyError if more than one match is found.

    """
    size = min_overlap
    found = []
    # Reverse complementing here provides massive speedup?
    rev_query = [seq.reverse_complement() for seq in query_list]
    while not found and not size > len(reference):
        for i, seq in enumerate(rev_query):
            if right:
                # FIXME: these getitems are the slowest part of assembly
                # Easiest speedup?
                if reference.endswith(seq[:size]):
                    found.append(i)
            else:
                if reference.startswith(seq[-size:]):
                    found.append(i)
        size += 1
    if len(found) > 1:
        raise AssemblyError("Ambiguous oligo binding")
    if not found:
        return None
    else:
        return found[0], size
