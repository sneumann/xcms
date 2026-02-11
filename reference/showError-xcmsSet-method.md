# Extract processing errors

If peak detection is performed with `findPeaks` setting argument
`stopOnError = FALSE` eventual errors during the process do not cause to
stop the processing but are recorded inside of the resulting `xcmsSet`
object. These errors can be accessed with the `showError` method.

## Usage

``` r
# S4 method for class 'xcmsSet'
showError(object, message. = TRUE, ...)
```

## Arguments

- object:

  An `xcmsSet` object.

- message.:

  Logical indicating whether only the error message, or the error itself
  should be returned.

- ...:

  Additional arguments.

## Value

A list of error messages (if `message. = TRUE`) or errors or an empty
list if no errors are present.

## Author

Johannes Rainer
