# Calibrant mass based calibration of chromatgraphic peaks

Calibrate peaks using mz values of known masses/calibrants. mz values of
identified peaks are adjusted based on peaks that are close to the
provided mz values. See details below for more information.

The `isCalibrated` function returns `TRUE` if chromatographic peaks of
the
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object `x` were calibrated and `FALSE` otherwise.

## Usage

``` r
CalibrantMassParam(
  mz = list(),
  mzabs = 1e-04,
  mzppm = 5,
  neighbors = 3,
  method = "linear"
)

isCalibrated(object)

# S4 method for class 'XCMSnExp'
calibrate(object, param)
```

## Arguments

- mz:

  a `numeric` or `list` of `numeric` vectors with reference mz values.
  If a `numeric` vector is provided, this is used for each sample in the
  `XCMSnExp` object. If a `list` is provided, it's length has to be
  equal to the number of samples in the experiment.

- mzabs:

  `numeric(1)` the absolute error/deviation for matching peaks to
  calibrants (in Da).

- mzppm:

  `numeric(1)` the relative error for matching peaks to calibrants in
  ppm (parts per million).

- neighbors:

  `integer(1)` with the maximal number of peaks within the permitted
  distance to the calibrants that are considered. Among these the mz
  value of the peak with the largest intensity is used in the
  calibration function estimation.

- method:

  `character(1)` defining the method that should be used to estimate the
  calibration function. Can be `"shift"`, `"linear"` (default) or
  `"edgeshift"`.

- object:

  An
  [XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object.

- param:

  The `CalibrantMassParam` object with the calibration settings.

## Value

For `CalibrantMassParam`: a `CalibrantMassParam` instance. For
`calibrate`: an
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object with chromatographic peaks being calibrated. **Be aware** that
the actual raw mz values are not (yet) calibrated, but **only** the
identified chromatographic peaks.

The `CalibrantMassParam()` function returns an instance of the
`CalibrantMassParam` class with all settings and properties set.

The `calibrate` method returns an
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object with the chromatographic peaks being calibrated. Note that
**only** the detected peaks are calibrated, but not the individual mz
values in each spectrum.

## Details

The method does first identify peaks that are close to the provided mz
values and, given that there difference to the calibrants is smaller
than the user provided cut off (based on arguments `mzabs` and `mzppm`),
their mz values are replaced with the provided mz values. The mz values
of all other peaks are either globally shifted (for `method = "shift"`
or estimated by a linear model through all calibrants. Peaks are
considered close to a calibrant mz if the difference between the
calibrant and its mz is `<= mzabs + mz * mzppm /1e6`.

**Adjustment methods**: adjustment function/factor is estimated using
the difference between calibrant and peak mz values only for peaks that
are close enough to the calibrants. The availabel methods are:

- `shift`: shifts the m/z of each peak by a global factor which
  corresponds to the average difference between peak mz and calibrant
  mz.

- `linear`: fits a linear model throught the differences between
  calibrant and peak mz values and adjusts the mz values of all peaks
  using this.

- `edgeshift`: performs same adjustment as `linear` for peaks that are
  within the mz range of the calibrants and shift outside of it.

For more information, details and examples refer to the
*xcms-direct-injection* vignette.

## Note

`CalibrantMassParam` classes don't have exported getter or setter
methods.

## Author

Joachim Bargsten, Johannes Rainer
