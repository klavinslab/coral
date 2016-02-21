Changelog for Coral
===================

0.5.0 (2016-02-20)
* Separated `ssDNA` (single-stranded) and `DNA` (implicitly double-stranded)
classes. Convert between them using `.to_ds()` and `.to_ss()` methods,
respectively.
* `coral.Primer` instances now contain `ssDNA` instances, not `DNA`.
* `feature_type` is now an optional argument when initializing new `Feature`
instances. The default is 'misc_feature'.
* Added `to_feature()` convenience method to `coral.DNA` for generating a
feature from a given `DNA` instance.
* Fixed an issue where `coral.analysis.tm` returned an arcane message when the
input is a non-DNA character.
* Fixed an issue where read_sequencing returned `ssDNA`.
* Overhauled alignment and Sanger sequencing analysis modules. Added support
for using the MAFFT command-line tool as an alignment method, located at
`coral.analysis.MAFFT`. Added MAFFT as an option for the Sanger analysis class.
Added a new function, `needle_msa`, that generates a reference-aligned MSA
representation of a set of Needleman-Wunsch pairwise alignments. Added
`coral.analysis.substitution_matrices` module, which adds a SubstitutionMatrix
class for easily specifying customized substitution matrices for
Needleman-Wunsch alignment as well as built-in matrices such as BLOSUM62, DNA,
and DNA_SIMPLE.
* dev note: started using zest.releaser to automate releases.


0.4.1-0.4.3
* Bugfixes to 0.4.0 (brown-bag style bugs, fixed immediately after elease).

0.4.0
* Renamed `coral.DNA.rotate()` method to `coral.DNA.rotate_to()`.
* Created new `coral.DNA.rotate()` method that rotates a sequence
'counter-clockwise', acting as a deque.
* Created new `coral.DNA.rotate_to_feature`, which rotates a sequence to a
given feature's start location.
* Created new `coral.DNA.excise` feature, which removes a feature's sequence
from a circular DNA object, generating a linear product (useful for swapping
out features).
* Improved `coral.DNA` `__getitem__` behavior.
* Added .material property to all sequence types for pseudo-type checking.
* Made `coral.DNA.top()` and `coral.DNA.bottom()` methods into properties
(`coral.DNA.top` and `coral.DNA.bottom`) that can be overwritten and accessed
directly.
* Replaced `coral.DNA.topology` and `coral.DNA.stranded` properties (which were
strings) with boolean-valued `.circular` and `.ds` values, respectively.
* Re-wrote (and renamed) the `coral.analysis.NUPACK` and
`coral.analysis.ViennaRNA` packages to be more feature-complete.
* Fixed an issue where re-running `coral.DNA.display()` in a Jupyter notebook
resulted in non-updated text labels.

0.3.3
* Added Python 2 version check to prevent installation on Python 3.

0.3.2
* Removed cython dependencies entirely to ease installation.

0.3.1
* Fixed issues with PCR simulation and annealing, can now handle all cases of
primer directionality and overlaps, linear and circular templates.
* Added pyx to package manifest for case where user already has cython
installed.

0.3.0
* Separated out annealing behavior into analysis function (`analysis.anneal`).
* Functions that use annealing (e.g. `reaction.pcr`) can now accept partial
annealing + overhang matches due to annealing overhaul.
* Gibson reactions now retain features of the inputs.
* Added `strip()` method to `coral.DNA`.
* Installation now works on Mac OS X.
* Made most dependencies optional.
* Fixed issue where features were not being copied, resulting in unexpected
behavior (assign by reference vs. value).
* Fixed an issue where slicing the last N bases of a sequence (e.g.
`y =x[-4:]`) would modify the feature locations of the parent (`x`).

0.2.1
* Added HISTORY.md (this file) changelog.
* Fixed version bump issue, added javascript to manifest, added
dev-requirements.txt.

0.2.0
* plasmid visualizations for iPython notebooks using `coral.DNA.display`.
* features are now searchable using `coral.DNA.select_features`.
* `seqio.read_dna` now keeps all feature qualifiers when reading genbank files
(thanks @eyu-bolthreads!).

0.1.0

Initial Release.
