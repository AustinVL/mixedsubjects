# mixedsubjects

An open-source R package for designing and analyzing mixed-subjects design randomized control trials (MSD-RCTs).

## Installation

```r
# install.packages("remotes")
remotes::install_github("your-org/mixedsubjects")
```

## Usage

```r
library(mixedsubjects)

# data frame with treatment assignment D, observed outcome Y (NA if unlabeled),
# and model predictions S0 and S1 for each unit
est <- msd_estimate(Y ~ D | S0 + S1, data = my_data, estimator = "dip_pp")

est$estimate
summary(est)
```
