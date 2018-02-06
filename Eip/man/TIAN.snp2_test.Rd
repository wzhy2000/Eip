\name{TIAN.snp2_test}
\alias{TIAN.snp2_test}

\title{2 SNPs analysis}

\description{ 
  Perform 2 SNPs analysis.
} 

\usage{ TIAN.snp2_test(sCsvFile, options=NULL, output_file=NULL) }

\arguments{
\item{sCsvFile}{CSV file as a data set.}
\item{options}{ parameters for the data file, see \code{\link{options}}.}
\item{output_file}{file name for the summary report. }
}

\value{
A list that includes the following results.

\item{ret2}{The results without correction for 2 SNP analysis.}
\item{ret_b2}{The results with Bonferroni correction for 2 SNP analysis.}
\item{ret_h2}{The results with HOLM correction for 2 SNP analysis.}
\item{ret_f2}{The results with FDR correction for 2 SNP analysis.}
}

\details{

2 SNP analysis is performed to get the following effector by chi-square test. 

\tabular{lll}{
No.\tab	Effectors\tab	Description for 2 SNPs A and B\cr
1\tab	a1\tab		Additive of SNP A\cr
2\tab	a2\tab		Additive of SNP B\cr
3\tab	d1\tab		Dominance of SNP A\cr
4\tab	d2\tab		Dominance of SNP B\cr
5\tab	a1a2\tab	Additive of SNP A Å~ Additive of SNP B\cr
6\tab	a1d2\tab	Additive of SNP A Å~ Dominance of SNP B\cr
7\tab	d1a2\tab	Dominance of SNP A Å~ Additive of SNP B\cr
8\tab	d1d2\tab	Dominance of SNP A Å~ Dominance of SNP B\cr
}

In the summary, the significant SNP combinations with epistatic effector are selected and 
exported into the log file. Besides the significant values,

a). Correlation figures for all epistatic effectors are exported into the PDF files.
b). Full data set are exported into the CSV files.


}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
library(Eip);

options  <- list(
	type        = "auto", 
	correction  = c("FDR, Holm"),
	ignoreRow  = NULL, 
	bHead      = TRUE);
	
FM.snp2_step( "data.csv", options=options );

}

\keyword{functions}
