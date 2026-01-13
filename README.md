# mixedsubjects

An open-source R package for designing and analyzing mixed-subjects design randomized control trials (MSD-RCTs).

## Who is this for?

Applied social scientists who want to combine a small human sample with inexpensive model-based predictions
(e.g., LLM outputs) to estimate treatment effects with higher precision while maintaining credible identification.

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
est <- msd_estimate(Y ~ D | S1 - S0, data = my_data, estimator = "dip_pp")

est$estimate
summary(est)
```

You can also use `msd_data()` to auto-detect column names by convention and still
override the formula when needed.

```r
auto_data <- msd_data(my_data)
msd_estimate(data = auto_data, estimator = "dip_pp")
```

## What estimator should I use?

- **DIM**: standard difference-in-means using only labeled outcomes.
- **GREG**: calibrates predictions using labeled data.
- **PPI++ / D–T**: tunes the calibration strength via cross-fitting.
- **DiP / DiP++ / D–T DiP**: exploit both-arm predictions to reduce noise further.

For more details, see the introductory vignette:

```r
vignette("mixedsubjects-intro")
```
