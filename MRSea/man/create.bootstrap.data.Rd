\name{create.bootstrap.data}
\alias{create.bootstrap.data}
\title{Create bootstrap data for non-parametric bootstrapping}
\usage{
create.bootstrap.data(dis.data, resample = "transect.id",
  rename = "segment.id", stratum = NULL)
}
\arguments{
  \item{dis.data}{Original data to be bootstrapped.
  Requires a column that matches argument \code{resample}
  exactly.}

  \item{resample}{Specifies the resampling unit for
  bootstrapping, default is \code{transect.id}. Must match
  a column name in \code{dis.data} exactly}

  \item{rename}{A vector of column names for which a new
  column needs to be created for the bootstrapped data.
  This defaults to \code{segment.id} for line transects,
  however others might be added A new column with new ids
  will automatically be created for the column listed in
  \code{resample}}

  \item{stratum}{The column name in \code{dis.data} that
  identifies the different strata. The default \code{NULL}
  returns un-stratified bootstrap data. If stratum is
  specified, this requires a column in \code{dis.data} that
  matches argument \code{stratum} exactly}
}
\value{
Returns one realisation of bootstrapped distance data. Note
that a new column (in addition to those listed under
argument \code{rename}) is created. If the default for
\code{resample} is used, a column with new unique ids
called \code{transect.id2}. Note that a new column is
created with renamed bootstrap resamples to preserve the
number of unique bootstrap resamples. If the default for
\code{resample} is used, i.e. \code{transect.id}, this new
column is called \code{transect.id2}. In addition, a new
column \code{segment.id2} is created which is required for
other bootstrap functions.
}
\description{
This function creates one realisation of bootstrapped data
based on \code{dis.data}. The default resampling unit is
\code{transect.id} which may be modified using the argument
\code{resample}.
}
\examples{
data(dis.data.re)
# run distance analysis to create NHATS
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
result<-ddf(dsmodel=~mcds(key="hn", formula=~1), data=dis.data.re, method="ds",
             meta.data=list(width=250))
dis.data.re<-create.NHAT(dis.data.re,result)

# bootstrap data without stratification
bootstrap.data<-create.bootstrap.data(dis.data.re)
# boostrap data with stratification (here by survey which is composed of
# season and impact)
dis.data.re$survey.id<-paste(dis.data.re$season,dis.data.re$impact,sep="")
bootstrap.data.str<-create.bootstrap.data(dis.data.re, stratum = "survey.id")
}

