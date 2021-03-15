#' Line transect data with no post-impact consequence
#'
#' A simulated dataset containing the observed perpendicular distances, the effort data and other variables of 
#' segmented line transect data. The variables are as follows:
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
#' \item \code{object} Id for detected object      
#' \item \code{distance} Perpendicular distance from the line
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 10771 rows and 12 variables
#' 
#' @name dis.data.no
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Line transect data with decrease post-impact
#'
#' A simulated dataset containing the observed perpendicular distances, the effort data and other variables of 
#' segmented line transect data. The variables are as follows:
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
#' \item \code{object} Id for detected object      
#' \item \code{distance} Perpendicular distance from the line
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 10759 rows and 12 variables
#' 
#' @name dis.data.de
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Line transect data with redistribution post-impact
#'
#' A simulated dataset containing the observed perpendicular distances, the effort data and other variables of 
#' segmented line transect data. The variables are as follows:
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
#' \item \code{object} Id for detected object      
#' \item \code{distance} Perpendicular distance from the line
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 10951 rows and 12 variables
#' 
#' @name dis.data.re
NULL


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Knot grid data for nearshore example
#'
#' @name knotgrid.ns
#' @docType data
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#' Knot grid data for offshore example
#'
#' @name knotgrid.off
#' @docType data
NULL


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Prediction grid data for no post-impact consequence
#'
#' A simulated dataset containing the true number of birds, the effort data and other variables of 
#' prediction grid data. The variables are as follows:
#' 
#' \itemize{
#' \item \code{area} area surveyed in the gridcell in km squared
#' \item \code{x.pos} spatial location in the horizontal axis in UTMs
#' \item \code{y.pos} spatial location in the vertical axis in UTMs
#' \item \code{depth} depth in m
#' \item \code{segment.id} Identifier for individual visits to the segment
#' \item \code{season} Numerical indicator for the four different seasons
#' \item \code{impact} Numerical indicator for before (0) and after (1) impact
#' \item \code{truth} number of birds
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 37928 rows and 8 variables
#' 
#' @name predict.data.no
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Prediction grid data for post-impact decrease
#'
#' A simulated dataset containing the true number of birds, the effort data and other variables of 
#' prediction grid data. The variables are as follows:
#' 
#' \itemize{
#' \item \code{area} area surveyed in the gridcell in km squared
#' \item \code{x.pos} spatial location in the horizontal axis in UTMs
#' \item \code{y.pos} spatial location in the vertical axis in UTMs
#' \item \code{depth} depth in m
#' \item \code{segment.id} Identifier for individual visits to the segment
#' \item \code{season} Numerical indicator for the four different seasons
#' \item \code{impact} Numerical indicator for before (0) and after (1) impact
#' \item \code{truth} number of birds
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 37928 rows and 8 variables
#' 
#' @name predict.data.de
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Prediction grid data for post-impact redistribution
#'
#' A simulated dataset containing the true number of birds, the effort data and other variables of 
#' prediction grid data. The variables are as follows:
#' 
#' \itemize{
#' \item \code{area} area surveyed in the gridcell in km squared
#' \item \code{x.pos} spatial location in the horizontal axis in UTMs
#' \item \code{y.pos} spatial location in the vertical axis in UTMs
#' \item \code{depth} depth in m
#' \item \code{segment.id} Identifier for individual visits to the segment
#' \item \code{season} Numerical indicator for the four different seasons
#' \item \code{impact} Numerical indicator for before (0) and after (1) impact
#' \item \code{truth} number of birds
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 37928 rows and 8 variables
#' 
#' @name predict.data.re
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Nearshore data with no effect of impact
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
#' @name ns.data.no
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Nearshore data with decrease post-impact
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
#' @name ns.data.de
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Prediction grid data for nearshore post-impact decrease
#'
#' A simulated prediction dataset containing the true counts, the effort data and other variables of 
#' grid data. The variables are as follows:
#' 
#' \itemize{
#' \item \code{x.pos} spatial location in the horizontal axis in UTMs
#' \item \code{y.pos} spatial location in the vertical axis in UTMs
#' \item \code{area} Area surveyed in the gridcell in km squared
#' \item \code{floodebb} 3 level factor covariate for tide state
#' \item \code{observationhour} hour of observation
#' \item \code{GridCode} identifier for the different grids that were surveyed
#' \item \code{Year} Year of the survey
#' \item \code{DavOfMonth} Day of the survey
#' \item \code{MonthOfYear} Month of the survey 
#' \item \code{impact} numerical indicator for before (0) and after (1) impact
#' \item \code{birds} true density of birds
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 27798 rows and 11 variables
#' 
#' @name ns.predict.data.de
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Prediction grid data for nearshore post-impact redistribution
#'
#' A simulated prediction dataset containing the true counts, the effort data and other variables of 
#' grid data. The variables are as follows:
#' 
#' \itemize{
#' \item \code{x.pos} spatial location in the horizontal axis in UTMs
#' \item \code{y.pos} spatial location in the vertical axis in UTMs
#' \item \code{area} Area surveyed in the gridcell in km squared
#' \item \code{floodebb} 3 level factor covariate for tide state
#' \item \code{observationhour} hour of observation
#' \item \code{GridCode} identifier for the different grids that were surveyed
#' \item \code{Year} Year of the survey
#' \item \code{DavOfMonth} Day of the survey
#' \item \code{MonthOfYear} Month of the survey 
#' \item \code{impact} numerical indicator for before (0) and after (1) impact
#' \item \code{birds} true density of birds
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A data frame with 27798 rows and 11 variables
#' 
#' @name ns.predict.data.re
NULL

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
