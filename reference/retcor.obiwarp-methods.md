# Align retention times across samples with Obiwarp

Calculate retention time deviations for each sample. It is based on the
code at <http://obi-warp.sourceforge.net/>. However, this function is
able to align multiple samples, by a center-star strategy.

For the original publication see

Chromatographic Alignment of ESI-LC-MS Proteomics Data Sets by Ordered
Bijective Interpolated Warping John T. Prince and, Edward M. Marcotte
Analytical Chemistry 2006 78 (17), 6140-6152

## Methods

- object = "xcmsSet":

  retcor(object, method="obiwarp", plottype = c("none", "deviation"),
  profStep=1, center=NULL, col = NULL, ty = NULL, response=1,
  distFunc="cor_opt", gapInit=NULL, gapExtend=NULL, factorDiag=2,
  factorGap=1, localAlignment=0, initPenalty=0)

## Arguments

- object:

  the `xcmsSet` object

- plottype:

  if `deviation` plot retention time deviation

- profStep:

  step size (in m/z) to use for profile generation from the raw data
  files

- center:

  the index of the sample all others will be aligned to. If
  center==NULL, the sample with the most peaks is chosen as default.

- col:

  vector of colors for plotting each sample

- ty:

  vector of line and point types for plotting each sample

- response:

  Responsiveness of warping. 0 will give a linear warp based on start
  and end points. 100 will use all bijective anchors

- distFunc:

  DistFunc function: cor (Pearson's R) or cor_opt (default, calculate
  only 10% diagonal band of distance matrix, better runtime), cov
  (covariance), prd (product), euc (Euclidean distance)

- gapInit:

  Penalty for Gap opening, see below

- gapExtend:

  Penalty for Gap enlargement, see below

- factorDiag:

  Local weighting applied to diagonal moves in alignment.

- factorGap:

  Local weighting applied to gap moves in alignment.

- localAlignment:

  Local rather than global alignment

- initPenalty:

  Penalty for initiating alignment (for local alignment only) Default: 0

Default gap penalties: (gapInit, gapExtend) \[by distFunc type\]: 'cor'
= '0.3,2.4' 'cov' = '0,11.7' 'prd' = '0,7.8' 'euc' = '0.9,1.8'

## Value

An `xcmsSet` object

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
