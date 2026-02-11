# Derive experimental design from file paths

The `phenoDataFromPaths` function builds a `data.frame` representing the
experimental design from the folder structure in which the files of the
experiment are located.

## Usage

``` r
phenoDataFromPaths(paths)
```

## Arguments

- paths:

  `character` representing the file names (including the full path) of
  the experiment's files.

## Note

This function is used by the *old* `xcmsSet` function to guess the
experimental design (i.e. group assignment of the files) from the
folders in which the files of the experiment can be found.

## Examples

``` r
## List the files available in the faahKO package
base_dir <- system.file("cdf", package = "faahKO")
cdf_files <- list.files(base_dir, recursive = TRUE, full.names = TRUE)
```
