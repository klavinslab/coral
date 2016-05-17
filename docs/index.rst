Examples
========

.. code-block:: python

    >>> # This creates Golden Gate cloning primers for any gene
    >>> # and then verifies the expected PCR product
    >>> prefix = cr.ssDNA('CCGGTCTCGATCG')
    >>> suffix = cr.ssDNA('CCGGTCTCTAGCA').reverse_complement()
    >>> overhangs = [prefix, suffix.reverse_complement()]
    >>> primers = cr.design.primers(my_gene,
                                    tm=65,
                                    overhangs=overhangs)
    >>> amplicon = cr.reaction.pcr(my_gene, prefix, suffix)
    >>> amplicon === prefix.to_ds() + my_gene + suffix.to_ds()
    True
