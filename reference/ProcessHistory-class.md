# Tracking data processing

Objects of the type ProcessHistory allow to keep track of any data
processing step in an metabolomics experiment. They are created by the
data processing methods, such as
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
and added to the corresponding results objects. Thus, usually, users
don't need to create them.

The `XProcessHistory` extends the `ProcessHistory` by adding a slot
`param` that allows to store the actual parameter class of the
processing step.

`processParam()`, `processParam<-`: get or set the parameter class from
an `XProcessHistory` object.

[`msLevel()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
returns the MS level on which a certain analysis has been performed, or
`NA` if not defined.

The `processType()` method returns a character specifying the processing
step *type*.

The `processDate()` extracts the start date of the processing step.

The `processInfo()` extracts optional additional information on the
processing step.

The `fileIndex()` extracts the indices of the files on which the
processing step was applied.

## Usage

``` r
# S4 method for class 'XProcessHistory'
processParam(object)

# S4 method for class 'XProcessHistory'
msLevel(object)

# S4 method for class 'ProcessHistory'
processType(object)

# S4 method for class 'ProcessHistory'
processDate(object)

# S4 method for class 'ProcessHistory'
processInfo(object)

# S4 method for class 'ProcessHistory'
fileIndex(object)
```

## Arguments

- object:

  A `ProcessHistory` or `XProcessHistory` object.

## Value

For `processParam`: a parameter object extending the `Param` class.

The `processType()` method returns a character string with the
processing step type.

The `processDate()` method returns a character string with the time
stamp of the processing step start.

The `processInfo()` method returns a character string with optional
additional informations.

The `fileIndex()` method returns a integer vector with the index of the
files/samples on which the processing step was applied.

## Slots

- `type`:

  `character(1)`: string defining the type of the processing step. This
  string has to match predefined values. Use
  [`processHistoryTypes()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  to list them.

- `date`:

  `character(1)`: date time stamp when the processing step was started.

- `info`:

  `character(1)`: optional additional information.

- `fileIndex`:

  integer of length 1 or \> 1 to specify on which samples of the object
  the processing was performed.

- `error`:

  (ANY): used to store eventual calculation errors.

- `param`:

  (Param): an object of type `Param` (e.g.
  [`CentWaveParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md))
  specifying the settings of the processing step.

- `msLevel:`:

  `integer` definining the MS level(s) on which the analysis was
  performed.

## Author

Johannes Rainer
