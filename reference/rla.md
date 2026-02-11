# Calculate relative log abundances

`rla` calculates the relative log abundances (RLA, see reference) on a
`numeric` vector.

## Usage

``` r
rla(x, group, log.transform = TRUE)

rowRla(x, group, log.transform = TRUE)
```

## Arguments

- x:

  `numeric` (for `rla`) or `matrix` (for `rowRla`) with the abundances
  (in natural scale) on which the RLA should be calculated.

- group:

  `factor`, `numeric` or `character` with the same length than `x` that
  groups values in `x`. If omitted all values are considered to be from
  the same group.

- log.transform:

  `logical(1)` whether `x` should be log2 transformed. Set to
  `log.transform = FALSE` if `x` is already in log scale.

## Value

`numeric` of the same length than `x` (for `rla`) or `matrix` with the
same dimensions than `x` (for `rowRla`).

## Details

The RLA is defines as the (log) abundance of an analyte relative to the
median across all abundances of the same group.

## References

De Livera AM, Dias DA, De Souza D, Rupasinghe T, Pyke J, Tull D,
Roessner U, McConville M, Speed TP. Normalizing and integrating
metabolomics data. *Anal Chem* 2012 Dec 18;84(24):10768-76. doi:
[10.1021/ac302748b](https://doi.org/10.1021/ac302748b)

## Author

Johannes Rainer

## Examples

``` r

x <- c(3, 4, 5, 1, 2, 3, 7, 8, 9)

grp <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)

rla(x, grp)
#>          1          1          1          2          2          2          3 
#> -0.4150375  0.0000000  0.3219281 -1.0000000  0.0000000  0.5849625 -0.1926451 
#>          3          3 
#>  0.0000000  0.1699250 
```
