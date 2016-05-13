Coral: SynBio Programming Library
=================================

.. container:: row

    .. container:: jumbocontainer

        .. container:: jumbotron

            .. raw:: html

                <div class="img"></div>

            .. container:: col-md-3 logo

                .. raw:: html

                   <span class="helper"></span>

                .. image:: coral_256.png

            .. container:: col-md-9 codeexample

                .. code-block:: python

                    >>> # Capture your lab's design workflow in simple Python code
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

.. _`introduction`:

Introduction
============

Try Coral
=========

For most people, installing Coral will be as simple as

.. code-block:: none

   pip install numpy
   pip install coral

Windows users may need to install vcpython2.7 first http://aka.ms/vcpython27

For more detailed instructions as well as what optional packages can be
installed (such as those required to do structural analysis), see
:ref:`installation`.

