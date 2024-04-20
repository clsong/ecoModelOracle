
# ecoModelOracle

<!-- badges: start -->
<!-- badges: end -->

The goal of ecoModelOracle is to validate or invalidate ecological models with ecological time series.

## Installation

You can install the development version of ecoModelOracle from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("clsong/ecoModelOracle")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(ecoModelOracle)
set.seed(1010)
validation_statistics(
  example_data,
  species_variable = "prey",
  growth_function = "prey",
  death_function = "prey*predator + prey^2"
)
```

