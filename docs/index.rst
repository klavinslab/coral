Coral: SynBio Programming Library
=================================

.. container:: row

    .. container:: col-md-8

        .. container:: row jumbo img-rounded

            .. image:: _static/brain_coral_blur_small.jpg

            .. container:: col-md-6 logo

                .. image:: coral_256.png

            .. container:: col-md-6 codeexample

                .. code-block:: python

                    primers = cr.cloning.design_primers(dna)
                    amplicon = cr.reaction.pcr(dna, primers)
                    amplicon == dna

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

