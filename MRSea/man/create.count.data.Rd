\name{create.count.data}
\alias{create.count.data}
\title{Aggregate distance data into count data}
\usage{
  create.count.data(dis.data, column.numbers = NULL)
}
\arguments{
  \item{dis.data}{Data frame containing distance data (one
  row for each detection). Expects a column \code{NHAT},
  i.e. size of detection divided by its probability of
  detection (see \code{create.NHAT}) and that ids in
  \code{segment.id} are unique regardless of what transect
  they belong to}

  \item{column.numbers}{Optional argument: vector of
  integers indicating which columns other than \code{NHAT}
  from \code{dis.data} should be retained in the returned
  data. Generally all columns from the segment and higher
  levels should be kept while those from the observation
  level should be discarded. If the default is used,
  \code{column.numbers=NULL}, the columns \code{distance},
  \code{object}, \code{size}, \code{distbegin} and
  \code{distend} from the observation level are
  automatically discarded. Note that for those columns from
  the observation layer that are kept, only the first
  recorded value will be transferred.}
}
\value{
  This function returns count data that is suited for
  second stage count modelling of distance sampling data.
  The data includes the columns \code{NHAT} and \code{area}
  which are the response and the offset required by
  functions concerned with second stage modelling from this
  package.
}
\description{
  This function creates a new data set where
  \code{dis.data} is aggregated for each visit to a segment
  (\code{segment.id}). The sum of the estimated number of
  individuals for each segment from \code{dis.data$NHAT} is
  given in the column \code{NHAT} in the new data. Only
  columns from the segment or higher layers should be
  carried over into \code{count.data} from \code{dis.data}.
  Use argument \code{column.numbers} to identify these.
}
\examples{
data(dis.data.re)
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
           meta.data=list(width=250))
dis.data.re<-create.NHAT(dis.data.re,result)
count.data<-create.count.data(dis.data.re)
}

