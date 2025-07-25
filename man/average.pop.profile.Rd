% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate.R
\name{average.pop.profile}
\alias{average.pop.profile}
\title{Calculate Average Population Profile}
\usage{
average.pop.profile(
  profile,
  avg.fn = base::mean,
  min.var.fn = std.dev.min,
  max.var.fn = std.dev.max
)
}
\arguments{
\item{profile}{A data structure containing the profile information (population simulation data).}

\item{avg.fn}{A function to compute the average. Defaults to `base::mean`.}

\item{min.var.fn}{A function to compute the minimum variance. Defaults to `std.dev.min`.}

\item{max.var.fn}{A function to compute the maximum variance. Defaults to `std.dev.max`.}
}
\value{
The average population profile as a result. The return type and structure depend
on the implementation details of the function.
}
\description{
This function calculates the average population profile based on the provided profile data.
}
\examples{
\dontrun{
# Example usage:
result <- average.pop.profile(profile_data)
print(result)
}

}
