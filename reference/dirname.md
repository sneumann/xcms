# Change the file path of an `OnDiskMSnExp` object

`dirname` allows to get and set the path to the directory containing the
source files of the `OnDiskMSnExp` (or
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md))
object.

## Usage

``` r
# S4 method for class 'OnDiskMSnExp'
dirname(path)

# S4 method for class 'OnDiskMSnExp'
dirname(path) <- value
```

## Arguments

- path:

  `OnDiskMSnExp`.

- value:

  `character` of length 1 or length equal to the number of files
  defining the new path to the files.

## Author

Johannes Rainer
