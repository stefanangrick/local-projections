lp.r
====

R code to calculate and plot impulse responses from local projections (work in progress)

Overview
--------

`lp.r` contains R code I wrote in 2017 to calculate and plot impulse responses from local projections. Conceptually this is based on work by Jordà (2005) and draws on Matlab code made available by Ramey & Zubairy (2018) and Jordà (2005).

Models are specified using the `lpirf()` function. Plots are created using [ggplot2](https://ggplot2.tidyverse.org/). `lp.r` is fully documented in a [roxygen2](https://cran.r-project.org/web/packages/roxygen2/)-compatible format. That said, please consider the code work in progress.

Alternative implementations
---------------------------

For a more complete, ready-to-use implementation of local projections, have a look at [lpirfs](https://cran.r-project.org/web/packages/lpirfs/) by [Philipp Adämmer](https://github.com/AdaemmerP/) (available on [CRAN](https://cran.r-project.org/web/packages/lpirfs/) and [GitHub](https://github.com/AdaemmerP/lpirfs)). Also make sure to check out the work of Barnichon & Brownlees (2018) on [smooth local projections](http://www.econ.upf.edu/~cbrownlees/), R code for which can be found on [GitHub](https://github.com/ctbrownlees/R-Package-lproj/).

References
----------

Barnichon, R., & Brownlees, C. T. (2018). *Impulse Response Estimation by Smooth Local Projections* \[SSRN Working Paper\]. <https://doi.org/10.2139/ssrn.2892508>

Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local Projections. *The American Economic Review*, *95*(1), 161–182. <https://doi.org/10.1257/0002828053828518>

Ramey, V. A., & Zubairy, S. (2018). Government Spending Multipliers in Good Times and in Bad: Evidence from US Historical Data. *Journal of Political Economy*, *126*(2), 850–901. <https://doi.org/10.1086/696277>
