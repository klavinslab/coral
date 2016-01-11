Changelog for Coral
===================

### 0.3.0
* Separated out annealing behavior into analysis function (`analysis.anneal`).
* Functions that use annealing (e.g. `reaction.pcr`) can now accept partial
annealing + overhang matches due to annealing overhaul.
* Gibson reactions now retain features of the inputs.
* Added `strip()` method to `coral.DNA`
* Installation now works on Mac OS X
* Made most dependencies optional
* Fixed issue where features were not being copied, resulting in unexpected
behavior (assign by reference vs. value)
* Fixed an issue where slicing the last N bases of a sequence (e.g.
`y =x[-4:]`) would modify the feature locations of the parent (`x`).

### 0.2.1
* Added HISTORY.md (this file) changelog
* Fixed version bump issue, added javascript to manifest, added dev-requirements.txt

### 0.2.0
* plasmid visualizations for iPython notebooks using `coral.DNA.display`
* features are now searchable using `coral.DNA.select_features`.
* `seqio.read_dna` now keeps all feature qualifiers when reading genbank files (thanks @eyu-bolthreads!)

### 0.1.0

Initial Release
