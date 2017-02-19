What is Coral?
==============

Coral is a Python library for designing DNA sequences from the ground up.

DNA is a string of molecules that form a genetic code made up of As, Ts, Gs,
and Cs, each letter representing a specific molecule in that code. Manipulating
this sequence is the basis of modern biotechnology: how you arrange the
sequence determines basically everything that happens to construct and operate
biological systems.

Biological engineers of all kinds spend a lot of their time deciding what
DNA sequences they want to build, and how they're going to build them. More
often than not, the decision on how to design a sequence comes down to applying
a some rules over and over again.

Unfortunately, every lab has its own way of doing things, including how to
design DNA sequences. This is where Coral comes in: Coral makes it easy to
write down your lab's process as Python code so that the steps can be easily
re-executed

Why use Coral?
==============

Coral programs are a record of your design process
--------------------------------------------------

How does your lab design a DNA primer (a short piece of DNA used for various
purposes))? How about a neighboring lab? What's the best method for your
application? This is a surprisingly difficult question to answer, because DNA
design information is handed down as "lab lore" from person to person. It's
like a game of telephone that costs hundreds to thousands of dollars every time
a mistake is made.

A Coral program is a text file that gives step-by-step instructions for how you
want to design your DNA sequence(s). Because it's text, you can send your
program to someone else so they can use your design method, include it in a
publication, or put it in version control to keep track of how it changes over
time.

Coral programs can automate the boring stuff
--------------------------------------------

If you've ever designed a DNA sequence, you've probably had to memorize and
apply a lot of rules. Flip this sequence here, make sure it's long enough, etc.
This quickly becomes tedious and eats up time better spent on other things.

Coral is written with design automation in mind, so you can take the tedious
design steps that you (or your colleagues) have locked in your heads and
quickly turn them into executable code.

Coral programs prevent common mistakes
--------------------------------------

Researchers spend an inordinate amount of time applying those rules,
effectively acting like a computer themselves. But we're human: we get bored,
we make mistakes. This is particularly easy to do with DNA design, because you
often have to apply many rules at once to a given task, and if you forget even
one of them the whole thing can be ruined. Even worse, you might not find out
about the mistake until you've already spent weeks trying to work with your
DNA, since it's often not obvious why your experiment isn't working at first.

Coral programs can prevent these kinds of mistakes in two ways: automation and
output validation. When your design process is a Coral program, it's impossible
for the computer to simply forget to apply a design step or rule.

Coral programs are extensible
-----------------------------

* Concisely express operations on sequences - spend your time on the logic of
  your application, not how to get BioPython/BioJava to be a design, rather than
  analysis, tool.
* Increase complexity through abstraction - can do searches through millions
  of sequences just as easily as picking one at random.
* No need to reinvent the wheel - cloning functions and structural analysis
  is built in.

Who should use Coral?
=====================

Because Coral programs can be used for research, engineering, and education,
Coral has a wide intended audience. Let's look at some typical potential users.

.. raw:: html

   <div class="row">
     <div class="col-md-6 text-center">
       <div class="panel panel-default">
         <div class="panel-heading">
           <i class="fa fa-eyedropper" aria-hidden="true"></i>
           <h4>Lab Person</h4>
         </div>
         <div class="panel-body">
           <p>You're an expert at inventing, building, and testing cool new biological systems. But when you're not doing repetitive tasks in the lab, you end up spending hours copying and pasting As, Ts, Gs, and Cs and hoping that you don't forget anything important.</p>
         </div>
       </div>
     </div>

     <div class="col-md-6 text-center">
       <div class="panel panel-default">
         <div class="panel-heading">
           <i class="fa fa-code" aria-hidden="true"></i>
           <h4>Coder</h4>
         </div>
         <div class="panel-body">
           <p>You have a coding or computer science background, and are interested in the hot new field of synthetic biology. But how do you even get started and how can you contribute, when your only resources are arcane rulebooks and word of mouth lab lore?</p>
         </div>
       </div>
     </div>

     <div class="col-md-6 text-center">
       <div class="panel panel-default">
         <div class="panel-heading">
           <i class="fa fa-graduation-cap" aria-hidden="true"></i>
           <h4>Professor</h4>
         </div>
         <div class="panel-body">
           <p>The expertise in your lab is constantly leaving: your grad students graduate and your post-docs get jobs. And every new student or post-doc in your lab has to learn or re-learn the basics of design for your lab, even though it just boils down to applying a few simple rules.</p>
         </div>
       </div>
     </div>

     <div class="col-md-6 text-center">
       <div class="panel panel-default">
         <div class="panel-heading">
           <i class="fa fa-industry" aria-hidden="true"></i>
           <h4>Industry</h4>
         </div>
         <div class="panel-body">
           <p>You develop and sell new technologies, so on top of all of the concerns of the lab researcher and professor, you also need to follow a schedule and come in under budget.</p>
         </div>
       </div>
     </div>
   </div>

Coral helps you
===============

* Keep a record of your design process

A Coral program is a text file that gives an exact specification of your design
process. Coral programs can be shared, modified, and re-executed.

* Automate the boring stuff

Coral programs can automate the boring stuff (like primer design) so you can
spend more of your time on the cool stuff. At the same time, Coral is flexible,
so you aren't stuck with a pre-existing design strategy: you can encode the
exact process that your lab prefers.

* Prevent human error

We all make mistakes, especially when doing boring stuff. Why make yourself
validate your sequences when Coral can do it for you? With Coral, you can
easily validate your designs through output checking. Adding overhangs with a
PCR? Run the reaction.pcr with your primers and template and compare it
programmatically to the expected product.

* Create complexity and scale



Save and reuse your design principles
-------------------------------------

* .. raw:: html

    <i class="fa fa-eyedropper"></i> Automate the boring parts of your design workflow.

* .. raw:: html

    <i class="fa fa-code"></i> Your biological designs will read like normal Python code. Learn biological sequence design from clear, concise operations.

* .. raw:: html

    <i class="fa fa-graduation-cap"></i> The design knowledge of your labs remains even when your grad students and post-docs leave.

* .. raw:: html

    <i class="fa fa-industry"></i> Spend less time on low-level details while amassing a databank of design processes.

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

License
=======

Coral is licensed under the permissive and industry-friendly MIT license.
