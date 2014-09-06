PatchMatch
================

PatchMatch algorithm for MATLAB. I referenced mostly [1].
(not "generalized" version.)

It is only for grayscale image! Modification of color image is future work #3.

.mex file (C++ codes) are not included. This is implemented for MATLAB only.
Though .mex or C++ is much faster then MATLAB and PatchMatch algorithm is a little hard to parallelize,
we have no intention of implementing in .mex or C++.



References
-------

* [1] Barnes, Connelly, et al. "PatchMatch: A randomized correspondence algorithm for structural image editing." ACM Transactions on Graphics-TOG 28.3 (2009): 24.
* [2] Barnes, Connelly, et al. "The generalized patchmatch correspondence algorithm." Computer Vision–ECCV 2010. Springer Berlin Heidelberg, 2010. 29-43.
