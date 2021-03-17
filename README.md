
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRSea

The latest version is **1.3.1** (15/03/2021)

The MRSea packages allows the fitting of **spatially adaptive regression
splines using SALSA**.

It was developed to examine animal survey data for signs of changes in
animal abundance and distribution following marine renewables
development. However the methods are suitable for a wide range of
applications.

The functions of this package can be used to analyse segmented line
transect (alongside the `mrds` package) or digital aerial data. The
package includes functions for fitting spatially adaptive one and 2D
smoothers using SALSA and CReSS. Euclidean or Geodesic distances can be
used to underpin the smoothed 2D surface and a choice of Gaussian or
exponential radial basis functions are available. Non-parametric
bootstrapping is available to estimate uncertainty. Several model
assessment tools are also available. Recent updates include the direct
estimation of robust standard errors, given a panel structure.

## Installation

You can install the latest bugfix release of MRSea from
[GitHub](https://github.com/lindesaysh/MRSea) with:

``` r
# install.packages("devtools")
devtools::install_github("lindesaysh/MRSea", ref="stable")
```

You can install the development version of inlabru from
[GitHub](https://github.com/lindesaysh/MRSea) with:

``` r
devtools::install_github("lindesaysh/MRSea", ref="master")
```

The package may also be downloaded as a `.zip` or `.tar.gz` from the
latest release

## Documentation

There are two vignettes available with the package:

-   Statistical Modelling of bird and cetacean distributions in offshore
    renewables development areas
-   MRSea: 2D Interaction Example

These are available here;

-   PDF/html versions on
    [Github](https://github.com/lindesaysh/MRSea/tree/master/inst/docs)
-   or `browseVignettes(package='MRSea')` if you installed from the
    `.tar.gz` in the latest release.
