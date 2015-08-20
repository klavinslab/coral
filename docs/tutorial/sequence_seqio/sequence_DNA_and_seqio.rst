
``sequence.DNA``
----------------

``sequence.DNA`` is the core data structure of ``pymbt``. It contains
DNA sequences and has methods (data-associatd functions) for their
manipulation. Most design functions in ``pymbt`` return a
``sequence.DNA`` object or something that contains a ``sequence.DNA``
object (like ``sequence.Primer``).

For this tutorial we need to read in a genbank sequence and create a new
``sequence.DNA`` object, so we import the ``seqio`` and ``sequence``
modules of ``pymbt``.

.. code:: python

    from pymbt import seqio, sequence
Simple sequences - slicing
~~~~~~~~~~~~~~~~~~~~~~~~~~

Making a new ``sequence.DNA`` object is very easy - by default, it takes
a string to initialize.

.. code:: python

    example_dna = sequence.DNA('atgcatgc')
    example_dna



.. parsed-literal::

    linear dsDNA:
    atgcatgc
    tacgtacg



``sequence.DNA`` does checking on the input to make sure it's really
DNA. For now, it supports only unambiguous letters - no N, Y, M, etc.

.. code:: python

    try:
        # Try to make an invalid DNA sequence - none of these characters (including space) are valid DNA
        not_dna = sequence.DNA('hello world')
    except ValueError as e:
        # If you put an invalid sequence in, it raises a ValueError with a custom message. 
        # This prints out the message without the long traceback.
        print e

.. parsed-literal::

    Encountered a non-DNA character


``sequence.DNA`` has most of the special python container methods
defined - what this means is that in most cases, you can slice
``sequence.DNA`` objects just as if they were string or lists.

.. code:: python

    # Reverse a sequence
    example_dna[::-1]



.. parsed-literal::

    linear dsDNA:
    cgtacgta
    gcatgcat



.. code:: python

    # Extract the first three bases
    example_dna[0:3]



.. parsed-literal::

    linear dsDNA:
    atg
    tac



.. code:: python

    # Extract the last seven bases
    example_dna[-7:]



.. parsed-literal::

    linear dsDNA:
    tgcatgc
    acgtacg



.. code:: python

    # Grab every other base starting at index 0
    example_dna[::2]



.. parsed-literal::

    linear dsDNA:
    agag
    tctc



.. code:: python

    # Is the sequence 'AT' in our sequence? How about 'AC'?
    print "'AT' is in our sequence: {}.".format("at" in example_dna)
    print "'AC' is in our sequence: {}.".format("ac" in example_dna)

.. parsed-literal::

    'AT' is in our sequence: True.
    'AC' is in our sequence: False.


Several other common special methods and operators are defined for
sequences - you can concatenate DNA (so long as it isn't circular) using
``+``, repeat linear sequences using ``*`` with an integer, check for
equality with ``==`` and ``!=`` (note: features, not just sequences,
must be identical), check the length with ``len(dna_object)``, etc.

Simple sequences - methods
~~~~~~~~~~~~~~~~~~~~~~~~~~

In addition to slicing, ``sequence.DNA`` provides methods for common
molecular manipulations. For example, reverse complementing a sequence
is a single call:

.. code:: python

    example_dna.reverse_complement()



.. parsed-literal::

    linear dsDNA:
    gcatgcat
    cgtacgta



An extremely important method is the ``.copy()`` method. It may seem
redundant to have an entire function for copying a sequence - why not
just assign a ``sequence.DNA`` object to a new variable? As in most
high-level languages, python does not actually copy entire objects in
memory when assignment happens - it just adds another reference to the
same data. The short of it is that the very common operation of
generating a lot of new variants to a sequence, or copying a sequence,
requires the use of a ``.copy()`` method. For example, if you want to
generate a new list of variants where an 'a' is substituted one at a
time at each part of the sequence, using ``.copy()`` returns the correct
result (the first example) while directly accessing example\_dna has
horrible consequences (the edits build up, as they all modify the same
piece of data sequentially):

.. code:: python

    # Correct way:
    copy_list = [example_dna.copy() for i, x in enumerate(example_dna)]
    for i, seq in enumerate(example_dna):
        copy_list[i][i] = 'a'
    print [str(x) for x in copy_list]
    print
    
    # Incorrect way:
    copy = example_dna.copy()
    copy_list = [copy for i, x in enumerate(example_dna)]
    for i, seq in enumerate(example_dna):
        copy_list[i][i] = 'a'
    print [str(x) for x in copy_list]

.. parsed-literal::

    ['atgcatgc', 'aagcatgc', 'atacatgc', 'atgaatgc', 'atgcatgc', 'atgcaagc', 'atgcatac', 'atgcatga']
    
    ['aaaaaaaa', 'aaaaaaaa', 'aaaaaaaa', 'aaaaaaaa', 'aaaaaaaa', 'aaaaaaaa', 'aaaaaaaa', 'aaaaaaaa']


An important fact about ``sequence.DNA`` methods and slicing is that
none of the operations modify the object directly - if we look at
example\_dna, it has not been reverse-complemented itself. Running
``example_dna.reverse_complement()`` outputs a new sequence, so if you
want to save your chance you need to assign a variable:

.. code:: python

    revcomp_dna = example_dna.reverse_complement()
    print example_dna
    print
    print revcomp_dna

.. parsed-literal::

    atgcatgc
    
    gcatgcat


You can also access important attributes of a ``sequence.DNA`` object
directly. The following are examples of how to get important sequences
or information about a sequence.

.. code:: python

    example_dna.top()  # The top strand - a simple python string in the 5' -> 3' orientation.



.. parsed-literal::

    'atgcatgc'



.. code:: python

    example_dna.bottom()  # The bottom strand - another python string, also in the 5' -> 3' orientation.



.. parsed-literal::

    'gcatgcat'



.. code:: python

    # Sequences are double stranded, or 'ds' by default. 
    # This is a directly accessible attribute, not a method, so () is not required.
    example_dna.stranded



.. parsed-literal::

    'ds'



.. code:: python

    # To change the 'strandedness', use the set_stranded method
    example_dna.set_stranded('ss')



.. parsed-literal::

    linear ssDNA:
    atgcatgc
    --------



.. code:: python

    # To access the topology of the strand, look at the .topology attribute.
    # Sequences can be either linear or circular.
    example_dna.topology



.. parsed-literal::

    'linear'



.. code:: python

    # You can switch between topologies using the .circularize and .linearize methods
    circular_dna = example_dna.circularize()
    circular_dna



.. parsed-literal::

    circular dsDNA:
    atgcatgc
    tacgtacg



.. code:: python

    # Linearization is more complex - you can choose the index at which to linearize a circular sequence.
    # This simulates a precise double stranded break at the index of your choosing.
    # The following example shows the difference between linearizing at index 0 (default) versus index 2
    # (python 0-indexes, so index 2 = 3rd base, i.e. 'g' in 'atg')
    print circular_dna.linearize()
    print
    print circular_dna.linearize(2)

.. parsed-literal::

    atgcatgc
    
    gcatgcat


.. code:: python

    # Sometimes you just want to rotate the sequence around - i.e. switch the top and bottom strands. 
    # For this, use the .flip() method
    example_dna.flip()



.. parsed-literal::

    linear dsDNA:
    gcatgcat
    cgtacgta



Complex DNA sequences
~~~~~~~~~~~~~~~~~~~~~

More complex sequences (like plasmids) have many annotated pieces and
benefit from other methods. ``sequence.DNA`` has many methods for
accessing and modifying complex sequences.

The following sequence is a plasmid that integrates at the *S.
cerevisiae* HO locus via ends-out integration, inserting the GEV
transactivator from McIsaac et al. 2011:

.. code:: python

    pKL278 = seqio.read_dna('../files_for_tutorial/maps/pMODKan-HO-pACT1GEV.ape')
Sequences have ``.name`` and ``.id`` attributes that are empty string by
default. By convention, you should fill them with appropriate strings
for your use case - the name is a human-readable name while id should be
a unique number or string.

.. code:: python

    pKL278.name  # Raw genbank name field - truncated due to genbank specifications



.. parsed-literal::

    'pMODKan_HO_pACT1GE'



Large sequences have summary representations, useful for getting a
general idea of which sequence you're manipulating

.. code:: python

    pKL278  # The sequence representation - shows ~40 bases on each side.



.. parsed-literal::

    circular dsDNA:
    tcgcgcgtttcggtgatgacggtgaaaacctctgacacat ... ttaacctataaaaataggcgtatcacgaggccctttcgtc
    agcgcgcaaagccactactgccacttttggagactgtgta ... aattggatatttttatccgcatagtgctccgggaaagcag



Complex sequences usually have annotations to categorize functional or
important elements. This plasmid has a lot of features - it's a yeast
shuttle vector, so it has sequences for propagating in *E. coli*,
sequences for integrating into the *S. cerevisiae* genome, sequences for
selection after transformation, and an expression cassette (promoter,
gene, terminator). In addition, it has common primer sites and annotated
subsequences.

.. code:: python

    pKL278.features  # Man that's way too many features



.. parsed-literal::

    [pGEX_3_primer 'misc_feature' feature (28 to 51) on strand 1,
     pMOD_t1pre 'misc_feature' feature (132 to 154) on strand 0,
     PmeI(1) 'misc_feature' feature (154 to 162) on strand 0,
     HO Targeting 1 'misc_feature' feature (162 to 725) on strand 0,
     pMOD_t1suf 'misc_feature' feature (725 to 755) on strand 0,
     KANMX Wach et al 1994 (genome del. project) 'misc_feature' feature (755 to 1152) on strand 0,
     KanMX CDS 'misc_feature' feature (1152 to 1962) on strand 0,
     KanMX terminator 'misc_feature' feature (1962 to 2200) on strand 0,
     M13 Forward (-47) primer 'primer_bind' feature (2200 to 2224) on strand 0,
     pACT1 'misc_feature' feature (2224 to 2885) on strand 0,
     Extra sequence not found in Gottschling map 'misc_feature' feature (2921 to 2932) on strand 0,
     GAL4(1-93) DBD 'misc_feature' feature (2940 to 3218) on strand 0,
     Differs from Gottschling map (backbone) 'misc_feature' feature (3218 to 3219) on strand 0,
     hER HBD 'misc_feature' feature (3255 to 4140) on strand 0,
     HSV1 VP16 'misc_feature' feature (4140 to 4344) on strand 0,
     Differs from Gottschling Map 'misc_feature' feature (4235 to 4236) on strand 0,
     stop codon 'misc_feature' feature (4344 to 4347) on strand 0,
     L2 'misc_feature' feature (4347 to 4377) on strand 0,
     T + pBluescript KS linker 'misc_feature' feature (4377 to 4399) on strand 0,
     CYC1 'terminator' feature (4403 to 4643) on strand 0,
     pYESTrp_rev primer 'primer_bind' feature (4412 to 4431) on strand 1,
     T7 EEV primer 'primer_bind' feature (4643 to 4665) on strand 0,
     upstream HO targeting 'misc_feature' feature (4665 to 5571) on strand 0,
     PmeI 'misc_feature' feature (5571 to 5579) on strand 0,
     PmeI site 'misc_feature' feature (5571 to 5579) on strand 0,
     M13R 'misc_feature' feature (5579 to 5619) on strand 0,
     origin-extended 'misc_feature' feature (5804 to 5889) on strand 0,
     ori 'misc_feature' feature (5889 to 6744) on strand 0,
     is a g in normal maps. 'misc_feature' feature (6426 to 6427) on strand 0,
     bla 'misc_feature' feature (6744 to 7605) on strand 0,
     AmpR promoter 'misc_feature' feature (7605 to 7684) on strand 0,
     New Feature 'misc_feature' feature (7684 to 7704) on strand 0]



With all of these features, manual slicing is inconvenient. The
``.extract()`` method makes it easy to isolate features from a complex
sequence:

.. code:: python

    # The beta-lactamase coding sequence, essential for propagation in *E. coli* on Amp/Carb media.
    # Note that it is transcribed in the direction of the bottom strand (right to left on this sequence)
    pKL278.extract('bla')



.. parsed-literal::

    linear dsDNA:
    ttaccaatgcttaatcagtgaggcacctatctcagcgatc ... aaaagggaataagggcgacacggaaatgttgaatactcat
    aatggttacgaattagtcactccgtggatagagtcgctag ... ttttcccttattcccgctgtgcctttacaacttatgagta



The ``.features`` attribute is just a list of ``sequence.Feature``
objects - you can add or remove them at will using standard python list
methods (like ``.pop`` and ``.append``). The use of ``sequence.Feature``
will be covered in a different tutorial.

In addition, you can efficiently match patterns in your sequence using
``.locate()``, which searches for a string on both the top and bottom
strands, returning a tuple containing the indexes of the matches (top
and bottom strands). In the following case, there are 8 matches for the
top strand and 5 for the bottom strand. In the case of a palindromic
query, only the top strand is reported.

.. code:: python

    pKL278.locate('atgcc')  # All occurrences of the pattern atgcc on the top and bottom strands (both 5'->3')



.. parsed-literal::

    [[78, 286, 1380, 2431, 4177, 4315, 7261, 7556], [737, 3718, 3828, 4131, 6939]]



Other methods
~~~~~~~~~~~~~

There are additional methods that can't be (easily) demonstrated in this
tutorial.

The ``.ape()`` method will launch ApE with your sequence if it is found
in your PATH environment variable. This enables some convenient analyses
that are faster with a GUI like simulating a digest or viewing the
general layout of annotations.
