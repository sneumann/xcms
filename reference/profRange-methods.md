# Specify a subset of profile mode data

Specify a subset of the profile mode matrix given a mass, time, or scan
range. Allow flexible user entry for other functions.

## Methods

- object = "xcmsRaw":

  `profRange(object, mzrange = numeric(), rtrange = numeric(), scanrange = numeric(), ...)`

## Arguments

- object:

  the `xcmsRaw` object

- mzrange:

  single numeric mass or vector of masses

- rtrange:

  single numeric time (in seconds) or vector of times

- scanrange:

  single integer scan index or vector of indecies

- ...:

  arguments to other functions

## Details

This function handles selection of mass/time subsets of the profile
matrix for other functions. It allows the user to specify such subsets
in a variety of flexible ways with minimal typing.

Because R does partial argument matching, `mzrange`, `scanrange`, and
`rtrange` can be specified in short form using `m=`, `s=`, and `t=`,
respectively. If both a `scanrange` and `rtrange` are specified, then
the `rtrange` specification takes precedence.

When specifying ranges, you may either enter a single number or a
numeric vector. If a single number is entered, then the closest single
scan or mass value is selected. If a vector is entered, then the range
is set to the [`range()`](https://rdrr.io/r/base/range.html) of the
values entered. That allows specification of ranges using shortened,
slightly non-standard syntax. For example, one could specify 400 to 500
seconds using any of the following: `t=c(400,500)`, `t=c(500,400)`, or
`t=400:500`. Use of the sequence operator (`:`) can save several
keystrokes when specifying ranges. However, while the sequence operator
works well for specifying integer ranges, fractional ranges do not
always work as well.

## Value

A list with the folloing items:

- mzrange:

  numeric vector with start and end mass

- masslab:

  textual label of mass range

- massidx:

  integer vector of mass indecies

- scanrange:

  integer vector with stat ane end scans

- scanlab:

  textual label of scan range

- scanidx:

  integer vector of scan range

- rtrange:

  numeric vector of start and end times

- timelab:

  textual label of time range

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
