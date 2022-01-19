[![R-CMD-check-bioc](https://github.com/sneumann/xcms/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/sneumann/xcms/actions?query=workflow%3AR-CMD-check-bioc)
[![codecov.io](https://codecov.io/github/sneumann/xcms/coverage.svg?branch=master)](https://codecov.io/github/sneumann/xcms?branch=master)
[![Years in Bioconductor](http://www.bioconductor.org/shields/years-in-bioc/xcms.svg)](http://www.bioconductor.org/packages/release/bioc/html/xcms.html)
[![Bioconductor release build status](http://www.bioconductor.org/shields/build/release/bioc/xcms.svg)](http://www.bioconductor.org/packages/release/bioc/html/xcms.html)
[![Bioconductor devel build status](http://www.bioconductor.org/shields/build/devel/bioc/xcms.svg)](http://www.bioconductor.org/checkResults/devel/bioc-LATEST/xcms.html)

# The `xcms` package (version >= 3)

<img align = "right" src="https://raw.githubusercontent.com/Bioconductor/BiocStickers/master/xcms/xcms.png" height="200">


Version >= 3 of the `xcms` package are updated and partially re-written versions
of the original `xcms` package. The version number *3* was selected to avoid
confusions with the `xcms2` (http://pubs.acs.org/doi/abs/10.1021/ac800795f)
software. While providing all of the original software's functionality, `xcms`
version >= 3 aims at:

1) Better integration into the Bioconductor framework:
  - Make use and extend classes defined in the `MSnbase` package.
  - Implement class versioning (Biobase's `Versioned` class).
  - Use `BiocParallel` for parallel processing.
2) Implementation of validation methods for all classes to ensure data
   integrity.
3) Easier and faster access to raw spectra data.
4) Cleanup of the source code:
  - Remove obsolete and redundant functionality (`getEIC`, `rawEIC` etc).
  - Unify interfaces, i.e. implement a layer of base functions accessing all
    analysis methods (which are implemented in C, C++ or R).
5) Using a more consistent naming scheme of methods that follows established
   naming conventions (e.g. `correspondence` instead of `grouping`).
6) Update, improve and extend the documentation.
7) Establishing a layer of base R-functions that interface all analysis
   methods. These should take M/Z, retention time (or scan index) and intensity
   values as input along with optional arguments for the downstream functions
   (implemented in C, C++ or R). The input arguments should be basic R objects
   (numeric vectors) thus enabling easy integration of analysis methods in other
   R packages.
8) The user interface's analysis methods should take the (raw) data object and a
   parameter class, that is used for dispatching to the corresponding analysis
   algorithm.

Discussions and suggestions are welcome:
https://github.com/sneumann/xcms/issues

For more information see the package [vignette](vignettes/xcms.Rmd).
