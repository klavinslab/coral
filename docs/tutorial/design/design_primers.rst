
Primer Design
-------------

One of the first things anyone learns in a molecular biology lab is how
to design primers. The exact strategies vary a lot and are sometimes
polymerase-specific. ``pymbt`` uses the Klavins' lab approach of
targeting a specific melting temperature (Tm) and nothing else, with the
exact Tm targeted being between 65°C and 72°C, the choice being personal
preference. ``pymbt`` currently defaults to 72°C on the Phusion
(modified Breslauer 1986) Tm calculator.

``pymbt.design_primer`` is a function that takes in a ``sequence.DNA``
object and rapidly finds the 5' subsequence that is closest to the
desired Tm (within a user-definable error range). If the entire sequence
would make a primer with too low of a Tm, a descriptive error is
produced.

For this tutorial, let's design primers that will amplify the gene EYFP.

.. code:: python

    import pymbt as pbt
First we read in a plasmid from Havens et al. 2012 and isolate the EYFP
sequence.

.. code:: python

    plasmid = pbt.seqio.read_dna("../files_for_tutorial/maps/pGP4G-EYFP.ape")
    eyfp = plasmid.extract("EYFP")
    print len(eyfp)
    eyfp

.. parsed-literal::

    717




.. parsed-literal::

    linear dsDNA:
    ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGC ... CGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG
    TACCACTCGTTCCCGCTCCTCGACAAGTGGCCCCACCACG ... GCGGCGGCCCTAGTGAGAGCCGTACCTGCTCGACATGTTC



Designing primers is straightforward - you just call
``design.design_primer`` with a ``sequence.DNA`` object as the input.

.. code:: python

    # Forward and reverse, one at a time using design_primer()
    forward = pbt.design.primer(eyfp)
    reverse = pbt.design.primer(eyfp.reverse_complement())
    # Both at once using design_primers()
    forward, reverse = pbt.design.primers(eyfp)
    # design_primer has many options, including adding overhangs
    custom_forward = pbt.design.primer(eyfp, tm=65, min_len=12, 
                                       tm_undershoot=1, tm_overshoot=1, 
                                       end_gc=True, tm_parameters="santalucia98", 
                                       overhang=pbt.DNA("GGGGGATCGAT"))
    print forward
    print
    print custom_forward

.. parsed-literal::

    ATGGTGAGCAAGGGCG
    
    GGGGGATCGATATGGTGAGCAAGGGCGAGGAGCTGTTCAC


Designing primers and getting a string output is just the first step in
primer design - we want to know whether the primers actually *work* and
write them out to a file. The point of programming DNA is that you
*never* copy and paste!

To simulate a PCR using the rules of molecular biology, use
``pymbt.reaction.pcr``. The output is a subsequence of the template DNA
- the features may not match the plasmid exactly (due to being truncated
by the PCR), but the sequences match. If a primer would bind in multiple
places (exact matches to the template), the pcr function will fail and
give a useful message.

You can check for identical sequences using python's built in ==
operator.

.. code:: python

    amplicon = pbt.reaction.pcr(plasmid, forward, reverse)
    amplicon == eyfp



.. parsed-literal::

    True



Now that we have verified that our primers should at least amplify the
DNA that we want, let's write out our primers to file so they can be
submitted to an oligo synthesis company.

.. code:: python

    # First we give our primers names (the `.name` attribute is empty by default)
    forward.name = "EYFP_forward"
    reverse.name = "EYFP_reverse"
    # Then we write to file - a csv (comma separated value file)
    pbt.seqio.write_primers([forward, reverse], "./designed_primers.csv", ["Forward EYFP primer", "Reverse EYFP primer"])
The csv file can then be opened in a spreadsheet application like Excel
or processed by a downstream program. This is the format of the csv:

.. code:: python

    import csv
    with open("./designed_primers.csv", "r") as csv_file:
        reader = csv.reader(csv_file)
        lines = [line for line in reader]
    for line in lines:
        print line

.. parsed-literal::

    ['name', 'sequence', 'notes']
    ['Forward EYFP primer', 'ATGGTGAGCAAGGGCG', '']
    ['Reverse EYFP primer', 'CTTGTACAGCTCGTCCATGCC', '']


.. code:: python

    