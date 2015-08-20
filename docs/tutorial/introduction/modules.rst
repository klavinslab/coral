
pymbt modules
-------------

pymbt has 7 modules: analysis, constants, database, design, reaction,
seqio, and sequence.

The modules have been split up by function - the activity that a user
wants to execute. For example, anything related to accessing scientific
databases is in the database module and activities related to designing
sequences are in the design module.

The modules are explicitly organized via their \_\_init\_\_.py files.
All this means is that anything available via pymbt.module.\* is usable
and hopefully useful. You can explore the functions and classes defined
for each module by reading more of the ipython documentation, sphinx
autodoc documentation, or interactively investigating modules in the
ipython notebook using tab completion and ? documentation. pymbt follows
the PEP 8 style guidelines on class and function names so that you can
differentiate between them - classes use CamelCase and functions use
lower\_case with underscores.

.. code:: python

    import pymbt as pbt  # alternative you can import each module by itself e.g. from pymbt import design
    dir(pbt)  # dir lists everything in a module/object. Ignore the double underscore items.



.. parsed-literal::

    ['DNA',
     'Feature',
     'Peptide',
     'Primer',
     'RNA',
     'RestrictionSite',
     '__builtins__',
     '__doc__',
     '__file__',
     '__name__',
     '__package__',
     '__path__',
     '_sequence',
     'analysis',
     'constants',
     'database',
     'design',
     'matplotlib',
     'reaction',
     'seqio',
     'simulation']



Top-level
~~~~~~~~~

In addition to the core modules, the top-level pymbt module provides the
core data structures used in pymbt - DNA, RNA, and Peptide (as well as
specialized classes like Primer).

.. code:: python

    dna = pbt.DNA("ATGC")
    print "DNA: {}".format(dna)
    # You can also run methods on the object - in this case, check if the DNA is palindromic
    print "Palindrome?: {}".format(dna.is_palindrome())
    print
    rna = pbt.RNA("AUGC")
    print "RNA: {}".format(rna)
    print
    pep = pbt.Peptide("mlnp")
    print "Peptide: {}".format(pep)

.. parsed-literal::

    DNA: ATGC
    Palindrome?: False
    
    RNA: AUGC
    
    Peptide: MLNP


As you can see above, to make DNA, RNA, or Peptide objects you just
invoke the correct sequence. command and give it a valid string as an
argument. Case does not matter, but precision does - only unambiguous
and valid DNA, RNA, or Peptide sequences are allowed. The sequence
module also contains special cases of DNA objects (Primer,
RestrictionSite, Feature), which are covered in detail later. You can
treat DNA, RNA, and Peptide objects much like strings or lists in
python, so addition, multiplication, slicing, and container logic are
all defined.

analysis
~~~~~~~~

The analysis module is focused on providing functions and classes for
analyzing DNA, RNA, and Peptides, focusing on information inherent to
the sequence (palindromes, repeats, melting temperatures), structural
information (Vienna RNA and NUPACK classes), and sequencing (Sanger
sequencing analysis).

.. code:: python

    # Example: finding the Tm of ATGCATGCATGCATGC according to the SantaLucia98 method.
    pbt.analysis.tm(dna * 4, parameters="santalucia98")



.. parsed-literal::

    48.03216557174494



constants
~~~~~~~~~

The constants module contains data - information that doesn't change
(i.e. is constant). This includes alphabets (sets of characters) that
define DNA, RNA, and peptides and other standards, such as the genbank
feature table.

database
~~~~~~~~

The database module is for accessing scientific databases. It currently
has limited functionality, talking only to the Rebase database of
restriction enzymes.

design
~~~~~~

The design module holds classes and functions for the design of new
constructs. The two most important functions are design\_primer and
gibson. The former designs primers for a given input sequence while the
latter designs Gibson primers for a whole series of input fragments.

reaction
~~~~~~~~

The reaction module simulates reactions relevant to cloning and basic
molecular genetics, including transcription, reverse transcription,
translation, exonuclease activity, extracting coding sequences,
digesting with restriction endonucleases, pcr, and Gibson assembly.

seqio
~~~~~

The seqio module is for sequence input/output - reading and writing
sequences. The module currently supports reading in individual sequences
(fasta or genbank) using read\_dna, reading in all the .ab1, .abi, and
.seq files in a directory using read\_sequencing, and writing DNA
objects to file (fasta or genbank).

.. code:: python

    