
Sequences
=========

``sequence.DNA``
----------------

``coral.DNA`` is the core data structure of ``coral``. If you are
already familiar with core python data structures, it mostly acts like a
container similar to lists or strings, but also provides further
object-oriented methods for DNA-specific tasks, like reverse
complementation. Most design functions in ``coral`` return a
``coral.DNA`` object or something that contains a ``coral.DNA`` object
(like ``coral.Primer``). In addition, there are related ``coral.RNA``
and ``coral.Peptide`` objects for representing RNA and peptide sequences
and methods for converting between them.

To get started with ``coral.DNA``, import ``coral``:

.. code:: ipython2

    import coral as cor

Your first sequence
~~~~~~~~~~~~~~~~~~~

Let's jump right into things. Let's make a sequence that's the first 30
bases of gfp from *A. victoria*. To initialize a sequence, you feed it a
string of DNA characters.

.. code:: ipython2

    example_dna = cor.DNA('atgagtaaaggagaagaacttttcactgga')
    display(example_dna)



.. parsed-literal::

    ATGAGTAAAGGAGAAGAACTTTTCACTGGA
    TACTCATTTCCTCTTCTTGAAAAGTGACCT


A few things just happened behind the scenes. First, the input was
checked to make sure it's DNA (A, T, G, and C). For now, it supports
only unambiguous letters - no N, Y, R, etc. Second, the internal
representation is converted to an uppercase string - this way, DNA is
displayed uniformly and functional elements (like annealing and overhang
regions of primers) can be delineated using case. If you input a non-DNA
sequence, a ``ValueError`` is raised.

For the most part, a ``sequence.DNA`` instance acts like a python
container and many string-like operations work.

.. code:: ipython2

    # Extract the first three bases
    display(example_dna[0:3])



.. parsed-literal::

    ATG
    TAC


.. code:: ipython2

    # Extract the last seven bases
    display(example_dna[-7:])



.. parsed-literal::

    CACTGGA
    GTGACCT


.. code:: ipython2

    # Reverse a sequence
    display(example_dna[::-1])



.. parsed-literal::

    AGGTCACTTTTCAAGAAGAGGAAATGAGTA
    TCCAGTGAAAAGTTCTTCTCCTTTACTCAT


.. code:: ipython2

    # Grab every other base starting at index 0
    display(example_dna[::2])



.. parsed-literal::

    AGGAAGGAACTTATG
    TCCTTCCTTGAATAC


.. code:: ipython2

    # Is the sequence 'AT' in our sequence? How about 'AC'?
    print "'AT' is in our sequence: {}.".format("AT" in example_dna)
    print "'ATT' is in our sequence: {}.".format("ATT" in example_dna)


.. parsed-literal::

    'AT' is in our sequence: True.
    'ATT' is in our sequence: False.


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

.. code:: ipython2

    example_dna.reverse_complement()




.. parsed-literal::

    TCCAGTGAAAAGTTCTTCTCCTTTACTCAT
    AGGTCACTTTTCAAGAAGAGGAAATGAGTA



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

.. code:: ipython2

    example_dna.copy()




.. parsed-literal::

    ATGAGTAAAGGAGAAGAACTTTTCACTGGA
    TACTCATTTCCTCTTCTTGAAAAGTGACCT



.. code:: ipython2

    # Incorrect way (editing shared + mutable sequence):
    example_dna = cor.DNA('atgagtaaaggagaagaacttttcactgga')
    variant_list = []
    for i, base in enumerate(example_dna):
        variant = example_dna
        variant.top[i] = 'A'
        variant.bottom[i] = 'T'
        variant_list.append(variant)
    print [str(x) for x in variant_list]
    
    print
    
    # Correct way (copy mutable sequence, then edit):
    example_dna = cor.DNA('atgagtaaaggagaagaacttttcactgga')
    variant_list = []
    for i, base in enumerate(example_dna):
        variant = example_dna.copy()
        variant.top[i] = 'A'
        variant.bottom[i] = 'T'
        variant_list.append(variant)
    print [str(x) for x in variant_list]


.. parsed-literal::

    ['AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA']
    
    ['ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'AAGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATAAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAATAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGAAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAAGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGAAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAAAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAAAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAAATTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACATTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTATTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTATCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTACACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTAACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCAATGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACAGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTAGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGAA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA']


An important fact about ``sequence.DNA`` methods and slicing is that
none of the operations modify the object directly (they don't mutate
their parent) - if we look at example\_dna, it has not been
reverse-complemented itself. Running
``example_dna.reverse_complement()`` outputs a new sequence, so if you
want to save your chance you need to assign a variable:

.. code:: ipython2

    revcomp_dna = example_dna.reverse_complement()
    display(example_dna)
    display(revcomp_dna)



.. parsed-literal::

    ATGAGTAAAGGAGAAGAACTTTTCACTGGA
    TACTCATTTCCTCTTCTTGAAAAGTGACCT



.. parsed-literal::

    TCCAGTGAAAAGTTCTTCTCCTTTACTCAT
    AGGTCACTTTTCAAGAAGAGGAAATGAGTA


You also have direct access important attributes of a ``sequence.DNA``
object. The following are examples of how to get important sequences or
information about a sequence.

.. code:: ipython2

    # The top strand - a simple python string in the 5' -> 3' orientation.
    example_dna.top




.. parsed-literal::

    ATGAGTAAAGGAGAAGAACTTTTCACTGGA



.. code:: ipython2

    # The bottom strand - another python string, also in the 5' -> 3' orientation.
    example_dna.bottom




.. parsed-literal::

    TCCAGTGAAAAGTTCTTCTCCTTTACTCAT



.. code:: ipython2

    # Sequences are double stranded, or 'ds' by default. 
    # This is a directly accessible attribute, not a method, so () is not required.
    example_dna.ds




.. parsed-literal::

    True



.. code:: ipython2

    # DNA can be linear or circular - check the boolean `circular` attribute.
    example_dna.circular




.. parsed-literal::

    False



.. code:: ipython2

    # You can switch between topologies using the .circularize and .linearize methods.
    # Circular DNA has different properties:
    #  1) it can't be concatenated to
    #  2) sequence searches using .locate will search over the current origin (e.g. from -10 to +10 for a 20-base sequence).
    circular_dna = example_dna.circularize()
    circular_dna.circular




.. parsed-literal::

    True



.. code:: ipython2

    # Linearization is more complex - you can choose the index at which to linearize a circular sequence.
    # This simulates a precise double stranded break at the index of your choosing.
    # The following example shows the difference between linearizing at index 0 (default) versus index 2
    # (python 0-indexes, so index 2 = 3rd base, i.e. 'g' in 'atg')
    print circular_dna.linearize()
    print
    print circular_dna.linearize(2)


.. parsed-literal::

    ATGAGTAAAGGAGAAGAACTTTTCACTGGA
    
    GAGTAAAGGAGAAGAACTTTTCACTGGAAT


.. code:: ipython2

    # Sometimes you just want to rotate the sequence around - i.e. switch the top and bottom strands. 
    # For this, use the .flip() method
    example_dna.flip()




.. parsed-literal::

    TCCAGTGAAAAGTTCTTCTCCTTTACTCAT
    AGGTCACTTTTCAAGAAGAGGAAATGAGTA


