#' Nearshore data with redistribution post-impact
#'
#' A simulated dataset containing the observed counts, the effort data and other variables of 
#' grid data. The variables are as follows:
#' 
#' \itemize{
#' \item \code{x.pos} spatial location in the horizontal axis in UTMs
#' \item \code{y.pos} spatial location in the vertical axis in UTMs
#' \item \code{area} area surveyed in the gridcell in km squared
#' \item \code{floodebb} 3 level factor covariate for tides
#' \item \code{observationhour} hour of observation
#' \item \code{GridCode} identifier for the different grids that were surveyed
#' \item \code{Year} Year of the survey
#' \item \code{DavOfMonth} Day of the survey
#' \item \code{MonthOfYear} Month of the survey 
#' \item \code{impact} numerical indicator for before (0) and after (1) impact
#' \item \code{birds} observed number of birds
#' \item \code{cellid} identifier for the individual records 
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 27798 rows and 12 variables
#' 
#' @name ns.data.re
NULL
