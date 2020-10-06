raytracing
================

[![Travis-CI Build
Status](https://img.shields.io/travis/com/salvatirehbein/raytracing?style=for-the-badge)](https://travis-ci.com/github/salvatirehbein/raytracing)
[![Codecov test
coverage](https://codecov.io/gh/salvatirehbein/raytracing/branch/master/graph/badge.svg)](https://codecov.io/gh/salvatirehbein/raytracing?branch=master)

## Raytracing Documentation

The identification of the atmospheric Rossby wave ray paths is of
paramount importance for scientists, climatologists, meteorologists and
students seeking for a better understanding of the dynamics of the
atmosphere. `raytracing` is designed to detects the triggering regions
of Rossby waves, their characteristics, and their paths, in a completely
automated way.

It relies on the theory described in:

  - [Hoskins and Karoly
    (1981)](https://journals.ametsoc.org/jas/article/38/6/1179/20507/The-Steady-Linear-Response-of-a-Spherical)
  - [Karoly
    (1983)](https://www.sciencedirect.com/science/article/abs/pii/0377026583900131)
  - [Hoskins and Ambrizzi
    (1993)](https://journals.ametsoc.org/jas/article/50/12/1661/23207/Rossby-Wave-Propagation-on-a-Realistic)
  - [Yang and Hoskins
    (1996)](https://journals.ametsoc.org/jas/article/53/16/2365/24038/Propagation-of-Rossby-Waves-of-Nonzero-Frequency)

For a brief review: Rehbein et al. (2020) \<\>

## Including Code

The `raytracing` include the below functions:

| Functions | Arguments                                                        | Description                                                                                                                         |
| :-------- | :--------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------- |
| calcUg    | betamz, umz, y, lat, K, w, a                                     | Calculates the zonal group velocity                                                                                                 |
| calcVg    | betamz, umz, y, lat, K, direction, tl, a                         | Calculates the meridional group velocity                                                                                            |
| betaks    | ifile, varname, latname, lonname, ofile, a, plots, show.warnings | Calculates the stationary Rossby wave number, meridional gradient of the absolute vorticity, and zonal wind in mercator coordinates |
| ray       | betamz, umz, lat, x0, y0, K, dt, itime, direction, tl, a         | Calculates the Rossby wave ray paths                                                                                                |
| ypos      | y, lat                                                           | Obtain the closest position of a latitude in a vector                                                                               |

**Table** **1** **-** `raytracing` functions, arguments, and its
description.

## Example

Simple usage example, reproduced from [Coelho et
al. (2015)](https://link.springer.com/article/10.1007/s00382-015-2800-1).

![](README_files/figure-gfm/onda-1.png)<!-- -->

The `ray` or the `ray_source` functions return a `data.table` exampled
below:

    head(a)

## Installing

Install development versions from github with

    library(devtools)
    install_github("salvatirehbein/raytracing")
