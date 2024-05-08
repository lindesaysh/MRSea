# MRSea 1.6

## Notes
* SALSA1D: 

- addition of parameter (`logfile`) to remove dependence on log file sink if `suppress.printout = TRUE`. Default is set to `FALSE` so if `suppress.printout = TRUE` then no log file is produced unless `logfile = TRUE`
- seed added to initialisation of candidate knotsites. The selection of knotsites is only stochastic if you have more than 800 unique values for a variable. With the seed set this allows repeatability of analyses. 
 
* SALSA2D:

- addition of parameter (`logfile`) to remove dependence on log file sink if `suppress.printout = TRUE`. Default is set to `FALSE` so if `suppress.printout = TRUE` then no log file is produced unless `logfile = TRUE`
 

* Vignettes:

- updates to the images

* Other: 

- `plotmeanVar`: the labelling for the Gaussian and Gamma alternative lines has been changed from "Poisson" to "1:1 line". 
- warning suppression has been added to `do.boostrap.cress.robust`.  There is already a try catch for "svd" warnings with an alternative used so not necessary to see the warnings. 

  
## Bug Fixes

* SALSA 1D:

* SALSA 2D:
  

* Other:

- `plotMeanVar` line labelling was the wrong way round for the Gaussian case. This has been remedied. 
- `cv.gamMRSea` occasional error relating to not finding CV0.  Code edited so that if there is an issue with CV calculation, `inf` is returned rather than an error. 


# MRSea 1.5


## Notes
* SALSA1D: 

  - update function to allow the inclusion of a variable that is a `Date` class.
  
* SALSA2D:
  
  - the `getRadiiChoices` and `getRadiiChoices.vario` functions have been superseded and merged into one function: `getRadiiSequence()`.
  - functions have been added to help create a factor based knot grid (`selectFctrKnots`) and get an associated infinity-block distance matrix (updates to the `makeDists` function).  There is also a function to select space-filled starting knots for each factor level (`selectFctrStartk`)

* Vignettes:
  
  - Introduction to MRSea: added section on using `make.gamMRSea`
  - Interactions: edited to use the new functions for factor based knots.
  - Renewables Case Study: vignette re-instated and completely updated for all new functions. 
  
  
* Other: 

  - The seed used for CV score calculation is appended as an attribute to the output of cv.gamMRSea. Additionally it is also an optional input parameter. 
  - The seed used for generating the CV folds (`getCVids`) is now appended as an attribute fo the returned vector of folds. 
  - the method for shortening the coefficient names in the summary output has been updated to be more efficient and to work in more situations. 
  
## Bug Fixes

* SALSA 1D:

* SALSA 2D:

  - The dispersion parameter estimate using `getDisperson` function was incorrectly specified as 1 for the Tweedie distribution.  This affected the dropping of "bad" knots in the initialise step resulting in more knots being dropped than should have been. This did not affect knots being subsequently added in the exchange step. The correct calculation is now included in the `getDispersion` function.
  - When using `initialise = FALSE` and specifying initial knot locations via `initialKnPos` code has been added to snap the initial locations to the nearest candidate knot locations in the knot grid. Previously it was assumed these initial locations were a subset of the candidate set. Now any locations may be used. 
  - When using an odd number of radii and with `chooseradii = FALSE` there may sometimes have been a mismatch from the output of SALSA2D and the output given from `runSALSA2D`. The default is an even number so mostly this will not have been an issue. 
  

* Other:

  - When using the Tweedie distribution, the power parameter is not stored in the model call and so issues arise when using a stored model object where this parameter is no longer available.  The parameter is stored within the model but updates have been made to `runSALSA1D` and `runSALSA2D` to ensure that the variance power and variance link are specified as numbers in the model call. 
  - The CV score was unchanged when changing K.  This occurred owing to the folds stored within the model object.  If `K` specified in `cv.gamMRSea` is different to that in the model object new folds are now created with a warning message.  The folds in the model object remain unchanged. Note that there is now an error message if the number of folds specified exceeds the number of unique panels. 
  - bug in the variance calculation for the gamma distribution in the `plotMeanVar` function. Additionally, the legend mismatched the data.
  - bug in the `runInfluence` and `timeInfluenceCheck`.  It worked fine for models without a 2D smooth but a specification error caused issues for models with one. Fixed and working as it should now. 
  


# MRSea 1.4


## Notes
* SALSA1D: 

  - Tweedie distribution and associated information criterion added (`AICtweedie` and `BICtweedie`)
  
* SALSA2D:

  - Tweedie distribution and associated information criterion added (`AICtweedie` and `BICtweedie`)

* Vignettes:

  - Added a vignette for using the Tweedie distribution in gamMRSea models
  - Added a vignette for model diagnostics with general modelling information and MRSea specific functions
  - Added a vignette for user specification of the range of radii choices using a variogram
  - Added a vignette for a brief introduction to distance sampling using the `mrds` package
  - Added an introduction to MRSea vignette.
  - The case study vignette has been removed for updating and will return soon. 
  
  
* Other: 

  - Mean variance plot function added. This allows the user to assess the suitability of the mean-variance distributional assumption for a number of distributions (Gaussian, Poisson, Quasi-Poisson, Gamma and Tweedie)
  - Minor edits to all diagnostic functions (ACF, influence measures, diagnostic plots, cumulative residuals)
  
## Bug Fixes

* SALSA 1D:

* SALSA 2D:

* Other:

  - ensure data is a `data.frame` (not `tibble`) in `create.bootstrap.data` and `plotCumRes` functions
  - minor bug fixes to diagnostic functions when no spline parameter object present in model


# MRSea 1.3.3


## Notes
* SALSA1D:

  - Update documentation for `runSALSA1D` to include `gap` parameter and `splines` parameter. Default settings are `gap=0` and `splines="bs"`.
  - Change specification of cyclic smooths (parameter now specified in `salsa1dlist`) and include additional option for a natural cubic spline.
  
* SALSA2D:

* Other:

  - Website created
  - Update to the vignettes included within the package and additional web only tutorials created.  Still under construction but included now is a basic usage vignette and information on the different splines. 
  - Add option for user specified label to plots in `runDiagnostics` function
  
## Bug Fixes

* SALSA 1D

 - Issue with using Q*IC and the update function for multiple 1D covariates. Problem resolved.

* SALSA 2D
  - Fix issue in initialise drop step when initial model has NA coefficients


# MRSea 1.3.2

## Notes
* Gap check added to initialise step to ensure gap parameter is obeyed during initialisation
* Remove storage of `models` object. Hang over from the early days of `MRSea` and taking up workspace memory.
* Update to the "Interaction" vignette to make it clearer and the two options of ways to include one. 
* Added an option to allow the user to specify their own sequence of radii. Use by adding `r_seq = ...` to `salsa2dlist`
 

## Bug Fixes

* CV fold bug fixed for SALSA1D. If correlated blocks specified, these were not maintained in the base model for use later in the algorithm. 
* reinstate `require(parallel)` in `do.bootstrap.cress.robust`


# MRSea 1.3.1

## Notes
* Package compiled using version of R >=4.0
* Removal of various now obsolete functions (predict.cress, getCVcress, getPvalues).  These have been superceded by (`predict.gamMRSea`, `cv.gamMRSea` and `anova.gamMRSea`)
* Updates to vignettes
 
## Bug Fixes

* Various fixes to comply with updated version of R

# MRSea 1.02


## Notes
* option added to runACF/plotacf functions to specify the the maximum lag to be plotted (maxlag=NULL).  This helps when there may be one long panel and you wish to see the detail in the smaller ones. 
* vignette updated
* drop step added to the initialise 2D step.  If the initial starting knots do not give a model that has converged then knots are dropped until convergence is achieved. Then the algorithm proceeds as normal (looping over exchange, improve, drop).
 
 
## Bug Fixes

* Fix to exchange step - error finding largest residual
* Fix issue clsandcov function which occurred when old factor levels remained in the panel variable (droplevels used to counter this)
* Fix issue in plotCumRes and runInfluence functions.  These work with gamMRSea type model now. 


# MRSea 1.01


## Notes

* Addition of an interaction term where the choice of knot locations for each level of the interaction are different.  Requires specialist set up of the knot grid and distance matrix.  These will be provided in a vignette shortly. 
* The user may also now specifiy which knots from the knotgrid they would like to use for initialising. See the help file for runSALSA2D for more. 

## Bug Fixes
  
* Fixed major bug in the improve step - issue with finding the nearest neighbours for a given knot.  


# MRSea 1.0-beta


## Notes

* Addition of GEODESIC distance calculation using the `makeDists` function. If you supply a polygon defining an exclusion area, geodesic distances are calculated. 
* Addition of cross-validation (cv.gamMRSea) as fitness measure
* Addition of option for gaussian (default) or exponential basis function (see runSALSA2D)
* Update to runPartialPlots function to allow inclusion (or not) of intercept uncertainty.
  
## Bug Fixes

* Fixed partial plot bug - issue with factor covariates which are characters.  This class now allowed.
  

# MRSea 0.99-beta


## Notes

* Major overhaul of package to include a new model class `gamMRSea`.
* Update to summary function to allow robust standard errors to be presented and used for hypothesis testing
* knotgrid for 2D smooth no longer required to be a regular grid
* calculation of basis radii absorbed into `runSALSA2D` function.
* predict.gamMRSea function call updated to have same names as predict.glm
* new cv.gamMRSea function.  This is the cv.glm function from the boot library edited to allow for use with gamMRSea models.

## Bug Fixes

* Fixed bug in summary and anova functions; term names now correct.
* `tol` option in runSALSA2D re-instated
