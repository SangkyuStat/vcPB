# vcPB: Varying Coefficient Peters-Belson Method for Longitudinal Data
A package for for estimating the disparity between a majority group and minority group based on the extended model of the Peters-Belson method. Our model is the first extension of Peters-Belson method to the longitudinal data. 

Furthermore, our method can set a modifiable variable to find the complicated association between other variables and the modifiable variable.

### Installation

The current version can be installed from source using the package `devtools`
```
devtools::install_github("SangkyuStat/vcPB")
```

### Usage Examples

`vc.pb` function provides three different types of models based on the different input arguments: `modifier` and time varying coefficients. The user needs to define `group` properly to measure the disparity between two groups in `group` variable, there should be 2 levels for this variable. 

If `modifier` is `NULL` (the default setting is `NULL`) and at least a time-varying variable exists, then the simple varying-coefficient Peters-Belson method using a gaussian kernel regression can be performed as below:
```
vc.pb(formula = response ~ (time varying variable | time variable) + variable, data = input data, group = disparity_group)
```
If `modifier` is not `NULL` and is a discrete variable, and at least a time-varying variable exists, then the modifiable varying-coefficient Peters-Belson method using a gaussian kernel regression can be performed as below:
```
vc.pb(formula = response ~ (time varying variable | time variable) + variable + discrete variable, data = input_ _data, group = disparity_group, modifier = "discrete variable")
```
If `modifier` is not `NULL` and is a continuous variable, and at least a time-varying variable exists, then the simple varying-coefficient Peters-Belson method using a gaussian kernel regression can be performed as below:
```
vc.pb(formula = response ~ (time varying variable | time variable) + variable + continuous variable, data = input_data, group = disparity_group, modifier = "continuous variable")
```
The type of modifier returns the different results.

The selection of bandwidths is essential and important for the kernel regression. If there is nothing given as initial values, we get and use the default marginal bandwidth from the function `KernSmooth::dpill`. For all models, `bandwidth_M`, `bandwidth_m`, `bandwidth_xM` and `bandwidth_xm` are essential. If `modifier` is not `NULL` and is a continuous variable, then `bandwidth_Z_M`, `bandwidth_Z_m`, `bandwidth_Z_xM` and `bandwidth_Z_xm` are needed more.

### Developing

- The cross-validation function for choosing the bandwidths will be developed.
- The original PB method will be included in the package.
- We are trying to develop other methods as well.
