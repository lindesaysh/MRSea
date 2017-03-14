#' Runs Test for Randomness
#'
#' @description This function performs the runs test for randomness. Users can choose whether to plot the correlation graph or not, and whether to test against two-sided, negative or positive correlation. NAs from the data are omitted. An empirical distribution may be used for the distribution of the test statistic under the null hypothesis of independence.
#' 
#' @param y a numeric vector of data values
#' @param plot.it logical flag.  If 'TRUE' then the graph will be plotted. If 'FALSE', then it is not plotted.
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "negative.correlated" or "positive.correlated"
#' @param emp.distribution vector containing the empirical distribution of test statistics under the nyll hypothesis. Generated using \code{\link{getEmpDistribution}}.
#' 
#' @details On the graph observations which are less than the sample median are represented by letter "A" in red color, and observations which are greater or equal to the sample median are represented by letter "B" in blue color.
#' 
#' @return 
#' A list with the following components.
#' 
#' \code{statistic}	the value of the standardized Runs statistic.
#' \code{p.value} the p-value for the test.
#' \code{data.name} a character string giving the names of the data.
#' \code{alternative} a character string describing the alternative hypothesis.
#' 
#' @references 
#' 
#' Mendenhall, W (1982), Statistics for Management and Economics, 4th Ed., 801-807, Duxbury Press, Boston.
#' 
#' J.L. Gastwirth; Y.R. Gel, W. L. Hui, V. Lyubchich, W. Miao and K. Noguchi (2015). lawstat: Tools for Biostatistics, Public Policy, and Law. R package version 3.0
#' 
#' @examples 
#' 
#' data(ns.data.re)
#' 
#' model<-gamMRSea(birds ~ observationhour + as.factor(floodebb) + as.factor(impact),  
#'               family='poisson', data=ns.data.re)
#'
#' runsTest(residuals(model))
#' 
#' ## Empirical distribution:
#' 
#' simData<-generateNoise(n=500, response=fitted(model), family='poisson')
#' 
#' empdist<-getEmpDistribution(500, simData, model, data=ns.data.re, plot=FALSE, 
#'                             returnDist=TRUE,dots=FALSE)
#' runsTest(residuals(model), emp.distribution=empdist)
#' 
#' @export
#' 
runsTest<-function (y, plot.it = FALSE, alternative = c("two.sided", "positive.correlated",
                                                        "negative.correlated"), emp.distribution = NULL)
{
  alternative <- match.arg(alternative)
  DNAME = deparse(substitute(y))
  y <- na.omit(y)
  med <- median(y, na.rm = TRUE)
  for (k in 2:length(y)) {
    if ((y[k] == med) & (y[k - 1] < med)) {
      y[k] = y[k - 1]
    }
    else if ((y[k] == med) & (y[k - 1] > med)) {
      y[k] = y[k - 1]
    }
  }
  q <- rep(0.05, length(y))
  p <- rep(-0.05, length(y))
  d <- y
  q[I(d < med) | I(d == med)] <- NA
  p[I(d >= med)] <- NA
  if (plot.it) {
    plot(q, type = "p", pch = "A", col = "red", ylim = c(-0.5,
                                                         0.5), xlim = c(1, length(y)), xlab = "", ylab = "")
    points(p, pch = "B", col = "blue")
    abline(h = 0)
  }
  m <- length(na.omit(q))
  n <- length(na.omit(p))
  R <- 1
  s <- sign(y - med)
  for (k in 1:(length(y) - 1)) {
    if (s[k] != s[k + 1]) {
      R <- R + 1
    }
  }
  E <- 1 + 2 * n * m/(n + m)
  s2 <- (2 * n * m * (2 * n * m - n - m))/((n + m)^2 * (n +
                                                          m - 1))
  statistic <- (R - E)/sqrt(s2)
  if (alternative == "positive.correlated") {
    if (is.null(emp.distribution)) {
      p.value = pnorm(statistic)
      METHOD = "Runs Test - Positive Correlated"
    }else{
      p.value=length(which(emp.distribution<statistic))/length(emp.distribution)
      METHOD = "Runs Test - Positive Correlated; Empirical Distribution"
    }
  }
  else if (alternative == "negative.correlated") {
    if (is.null(emp.distribution)) {
      p.value = 1 - pnorm(statistic)
      METHOD = "Runs Test - Negative Correlated"
    }else{
      p.value=length(which(emp.distribution>statistic))/length(emp.distribution)
      METHOD = "Runs Test - Negative Correlated; Empirical Distribution"
    }
  }
  else {
    if (is.null(emp.distribution)) {
      p.value = 2 * min(pnorm(statistic), 1 - pnorm(statistic))
      alternative = "two.sided"
      METHOD = "Runs Test - Two sided"
    }
    else {
      if(alternative == "two.sided"){
        p.value=2 * min(length(which(emp.distribution>statistic))/length(emp.distribution), length(which(emp.distribution<statistic))/length(emp.distribution))
        alternative = "two.sided"
        METHOD = "Runs Test - Two sided; Empirical Distribution"
      }
    }
  }
  STATISTIC = statistic
  names(STATISTIC) = "Standardized Runs Statistic"
  structure(list(statistic = STATISTIC, p.value = p.value,
                 method = METHOD, data.name = DNAME), class = "htest")
}
