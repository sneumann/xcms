# Extract isolation window target m/z definition

`isolationWindowTargetMz` extracts the isolation window target m/z
definition for each spectrum in `object`.

## Usage

``` r
# S4 method for class 'OnDiskMSnExp'
isolationWindowTargetMz(object)
```

## Arguments

- object:

  [MSnbase::OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
  object.

## Value

a `numeric` of length equal to the number of spectra in `object` with
the isolation window target m/z or `NA` if not specified/available.

## Author

Johannes Rainer
