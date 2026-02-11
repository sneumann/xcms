# Generate unique group (feature) names based on mass and retention time

`groupnames` generates names for the identified features from the
correspondence analysis based in their mass and retention time. This
generates feature names that are equivalent to the group names of the
*old* user interface (aka xcms1).

## Usage

``` r
# S4 method for class 'XCMSnExp'
groupnames(object, mzdec = 0, rtdec = 0, template = NULL)
```

## Arguments

- object:

  `XCMSnExp` object containing correspondence results.

- mzdec:

  `integer(1)` with the number of decimal places to use for m/z (
  defaults to `0`).

- rtdec:

  `integer(1)` with the number of decimal places to use for the
  retention time (defaults to `0`).

- template:

  `character` with existing group names whose format should be emulated.

## Value

`character` with unique names for each feature in `object`. The format
is `M(m/z)T(time in seconds)`.

## See also

[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md).
