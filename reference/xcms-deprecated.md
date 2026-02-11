# Deprecated functions in package ‘xcms’

These functions are provided for compatibility with older versions of
‘xcms’ only, and will be defunct at the next release.

## Details

The following functions/methods are deprecated.

- `profBin`, `profBinM`, `profBinLin`, `profBinLinM`, `profBinLinBase`,
  `profBinLinBaseM` have been deprecated and
  [`binYonX`](https://sneumann.github.io/xcms/reference/binYonX.md) in
  combination with
  [`imputeLinInterpol`](https://sneumann.github.io/xcms/reference/imputeLinInterpol.md)
  should be used instead.

- `extractMsData`: replaced by `as(x, "data.frame")`.

- `plotMsData`: replaced by `plot(x, type = "XIC")`.
