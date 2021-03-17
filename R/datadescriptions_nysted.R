# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Polygon of the Nysted coastline
#'
#' The variables are as follows:
#' 
#' \itemize{
#' \item \code{x} spatial location in the horizontal axis in UTMs. Unit is kilometers.
#' \item \code{y} spatial location in the vertical axis in UTMs. Unit is kilometers.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 46 rows and 2 variables
#' 
#' @name nysted.coast
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Polygon of the outline of the Nysted study area
#'
#' The variables are as follows:
#' 
#' \itemize{
#' \item \code{x.pos} spatial location in the horizontal axis in UTMs. Unit is meters.
#' \item \code{y.pos} spatial location in the vertical axis in UTMs. Unit is meters.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 2210 rows and 2 variables
#' 
#' @name nysted.studybnd
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Prediction grid data for Nysted Data
#'
#' A simulated prediction dataset containing the true counts, the effort data and other variables of
#' grid data. The variables are as follows:
#' 
#' \itemize{
#' \item \code{area} area surveyed in the gridcell in km squared
#' \item \code{x.pos} spatial location in the horizontal axis in UTMs (km)
#' \item \code{y.pos} spatial location in the vertical axis in UTMs (km)
#' \item \code{depth} depth in m
#' \item \code{segment.id} Identifier for individual visits to the segment
#' \item \code{season} Numerical indicator for the four different seasons
#' \item \code{impact} Numerical indicator for before (0) and after (1) impact
#' \item \code{truth.re} number of birds in the redistribution simulation
#' \item \code{truth.de} number of birds in the site-wide decrease simulation
#' \item \code{truth.no} number of birds in the no change simulation
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 37928 rows and 10 variables
#' 
#' @name nysted.predictdata
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Nysted Data
#'
#' A simulated dataset containing the distance corrected `dis.data.re`. The variables are as follows:
#' 
#' \itemize{
#' \item \code{transect.id} Identifier for the individual visits to the transects
#' \item \code{transect.label} Labels for transects
#' \item \code{season} Numerical indicator for the four different seasons
#' \item \code{impact} Numerical indicator for before (0) and after (1) impact
#' \item \code{segment.id} Identifier for individual visits to the segment
#' \item \code{segment.label} Label for segments
#' \item \code{length} Length of segment in km
#' \item \code{x.pos} spatial location in the horizontal axis in UTMs
#' \item \code{y.pos} spatial location in the vertical axis in UTMs
#' \item \code{depth} Depth in m
#' \item \code{area} area of each segment in kmsq
#' \item \code{NHAT} Bird counts in each segment
#' \item \code{response} Same as for `NHAT` 
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 9232 rows and 12 variables
#' 
#' @name nysted.analysisdata
NULL