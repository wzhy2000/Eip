\name{TIAN.snp1_test}
\alias{TIAN.snp1_test}

\title{1 SNP analysis}

\description{ 
  Perform 1 SNP analysis.
} 

\usage{ TIAN.snp1_test(sCsvFile, options=NULL, output_file=NULL) }

\arguments{
\item{sCsvFile}{CSV file as a data set.}
\item{options}{ Parameters for the data file, see \code{\link{options}}.}
\item{output_file}{File name for the summary report. }
}

\value{
A list that includes the following results.

\item{ret1}{The result without correction for 1 SNP analysis.}
\item{ret_b1}{The results with Bonferroni correction for 1 SNP analysis.}
\item{ret_h1}{The results with HOLM correction for 1 SNP analysis.}
\item{ret_f1}{The results with FDR correction for 1 SNP analysis.}

}

\details{

1 SNP analysis is performed to get addititive and dominance effector by chi-square test. 
No figures are exported.

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{

library(Eip);

options  <- list(
	type        = "auto", 
	correction  = c("FDR, Holm"),
	ignoreRow  = NULL, 
	bHead      = TRUE);
	
FM.snp1_step( "data.csv", options=options );

}

\keyword{functions}
