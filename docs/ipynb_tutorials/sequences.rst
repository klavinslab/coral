
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

.. code:: python

    import coral as cor

Your first sequence
~~~~~~~~~~~~~~~~~~~

Let's jump right into things. Let's make a sequence that's the first 30
bases of gfp from *A. victoria*. To initialize a sequence, you feed it a
string of DNA characters.

.. code:: python

    example_dna = cor.DNA('atgagtaaaggagaagaacttttcactgga')
    display(example_dna)



.. parsed-literal::

    linear dsDNA:
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

.. code:: python

    # Extract the first three bases
    display(example_dna[0:3])



.. parsed-literal::

    linear dsDNA:
    ATG
    TAC


.. code:: python

    # Extract the last seven bases
    display(example_dna[-7:])



.. parsed-literal::

    linear dsDNA:
    CACTGGA
    GTGACCT


.. code:: python

    # Reverse a sequence
    display(example_dna[::-1])



.. parsed-literal::

    linear dsDNA:
    AGGTCACTTTTCAAGAAGAGGAAATGAGTA
    TCCAGTGAAAAGTTCTTCTCCTTTACTCAT


.. code:: python

    # Grab every other base starting at index 0
    display(example_dna[::2])



.. parsed-literal::

    linear dsDNA:
    AGGAAGGAACTTATG
    TCCTTCCTTGAATAC


.. code:: python

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

.. code:: python

    example_dna.reverse_complement()




.. parsed-literal::

    linear dsDNA:
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

.. code:: python

    # Incorrect way:
    copy = example_dna.copy()
    copy_list = [copy for i, x in enumerate(example_dna)]
    for i, seq in enumerate(example_dna):
        copy_list[i][i] = 'a'
    print [str(x) for x in copy_list]
    
    # Correct way:
    copy_list = [example_dna.copy() for i, x in enumerate(example_dna)]
    for i, seq in enumerate(example_dna):
        copy_list[i][i] = 'a'
    print [str(x) for x in copy_list]
    print


.. parsed-literal::

    ['AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA']
    ['ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'AAGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATAAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAATAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGAAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAAGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGAAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAAAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAAAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAAATTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACATTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTATTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTATCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTACACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTAACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA', 'ATGAGTAAAGGAGAAGAACTTTTCAATGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACAGGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTAGA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGAA', 'ATGAGTAAAGGAGAAGAACTTTTCACTGGA']
    


An important fact about ``sequence.DNA`` methods and slicing is that
none of the operations modify the object directly (they don't mutate
their parent) - if we look at example\_dna, it has not been
reverse-complemented itself. Running
``example_dna.reverse_complement()`` outputs a new sequence, so if you
want to save your chance you need to assign a variable:

.. code:: python

    revcomp_dna = example_dna.reverse_complement()
    display(example_dna)
    display(revcomp_dna)



.. parsed-literal::

    linear dsDNA:
    ATGAGTAAAGGAGAAGAACTTTTCACTGGA
    TACTCATTTCCTCTTCTTGAAAAGTGACCT



.. parsed-literal::

    linear dsDNA:
    TCCAGTGAAAAGTTCTTCTCCTTTACTCAT
    AGGTCACTTTTCAAGAAGAGGAAATGAGTA


You also have direct access important attributes of a ``sequence.DNA``
object. The following are examples of how to get important sequences or
information about a sequence.

.. code:: python

    # The top strand - a simple python string in the 5' -> 3' orientation.
    example_dna.top()




.. parsed-literal::

    'ATGAGTAAAGGAGAAGAACTTTTCACTGGA'



.. code:: python

    # The bottom strand - another python string, also in the 5' -> 3' orientation.
    example_dna.bottom()




.. parsed-literal::

    'TCCAGTGAAAAGTTCTTCTCCTTTACTCAT'



.. code:: python

    # Sequences are double stranded, or 'ds' by default. 
    # This is a directly accessible attribute, not a method, so () is not required.
    example_dna.stranded




.. parsed-literal::

    'ds'



.. code:: python

    # To change the 'strandedness', use to_ss() or to_ds()
    example_dna.to_ss()




.. parsed-literal::

    linear ssDNA:
    ATGAGTAAAGGAGAAGAACTTTTCACTGGA
    ------------------------------



.. code:: python

    # To access the topology of the strand, look at the .topology attribute.
    # Sequences can be either linear or circular.
    example_dna.topology




.. parsed-literal::

    'linear'



.. code:: python

    # You can switch between topologies using the .circularize and .linearize methods.
    # Circular DNA has different properties:
    #  1) it can't be concatenated to
    #  2) sequence searches using .locate will search over the current origin (e.g. from -10 to +10 for a 20-base sequence).
    circular_dna = example_dna.circularize()
    circular_dna




.. parsed-literal::

    circular dsDNA:
    ATGAGTAAAGGAGAAGAACTTTTCACTGGA
    TACTCATTTCCTCTTCTTGAAAAGTGACCT



.. code:: python

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


.. code:: python

    # Sometimes you just want to rotate the sequence around - i.e. switch the top and bottom strands. 
    # For this, use the .flip() method
    example_dna.flip()




.. parsed-literal::

    linear dsDNA:
    TCCAGTGAAAAGTTCTTCTCCTTTACTCAT
    AGGTCACTTTTCAAGAAGAGGAAATGAGTA



