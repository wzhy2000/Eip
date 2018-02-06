\name{TIAN.snp3_test}
\alias{TIAN.snp3_test}

\title{3 SNPs analysis}

\description{ 
  Perform 3 SNPs analysis.
} 

\usage{ TIAN.snp3_test( sCsvFile, options=NULL, output_file=NULL ) }

\arguments{
\item{sCsvFile}{CSV file as a data set.}
\item{options}{parameters for the data file, see \code{\link{options}}.}
\item{output_file}{file name for the summary report. }
}

\value{
A list that includes the following results.

\item{ret3}{The results without correction for 3 SNP analysis.}
\item{ret_b3}{The results with Bonferroni correction for 3 SNP analysis.}
\item{ret_h3}{The results with HOLM correction for 3 SNP analysis.}
\item{ret_f3}{The results with FDR correction for 3 SNP analysis.}

}

\details{

3 SNP analysis is performed to get the following effector by chi-square test. 

\tabular{lll}{
No.\tabEffectors\tab Description for 3 SNPs A , B and C\cr
1\tab	a1\tab		Additive of SNP A\cr
2\tab	a2\tab		Additive of SNP B\cr
3\tab	a3\tab		Additive of SNP C\cr
4\tab	d1\tab		Dominance of SNP A\cr
5\tab	d2\tab		Dominance of SNP B\cr
6\tab	d3\tab		Dominance of SNP C\cr
7\tab	a1a2\tab	Additive of SNP A Å~ Additive of SNP B\cr
8\tab	a1a3\tab	Additive of SNP A Å~ Additive of SNP C\cr
9\tab	a2a3\tab	Additive of SNP B Å~ Additive of SNP C\cr
10\tab	d1d2\tab	Dominance of SNP A Å~ Dominance of SNP B\cr
11\tab	d2d3\tab	Dominance of SNP B Å~ Dominance of SNP C\cr
13\tab	d1d3\tab	Dominance of SNP A Å~ Dominance of SNP C\cr
12\tab	a1d2\tab	Additive of SNP A Å~ Dominance of SNP B\cr
14\tab	a1d3\tab	Additive of SNP A Å~ Dominance of SNP C\cr
15\tab	d1a2\tab	Dominance of SNP A Å~ Additive of SNP B\cr
16\tab	d1a3\tab	Dominance of SNP A Å~ Additive of SNP C\cr
17\tab	a2d3\tab	Additive of SNP B Å~ Dominance of SNP C\cr
18\tab	d2a3\tab	Dominance of SNP B Å~ Additive of SNP C\cr
19\tab	a1a2a3\tab	Additive of SNP A Å~ Additive of SNP B Å~ Additive of SNP C\cr
20\tab	a1a2d3\tab	Additive of SNP A Å~ Additive of SNP B Å~ Dominance of SNP C\cr
21\tab	a1d2a3\tab	Additive of SNP A Å~ Dominance of SNP B Å~ Additive of SNP C\cr
22\tab	a1d2d3\tab	Additive of SNP A Å~ Dominance of SNP B Å~ Dominance of SNP C\cr
23\tab	d1a2a3\tab	Dominance of SNP A Å~ Additive of SNP B Å~ Additive of SNP C\cr
24\tab	d1a2d3\tab	Dominance of SNP A Å~ Additive of SNP B Å~ Dominance of SNP C\cr
25\tab	d1d2a3\tab	Dominance of SNP A Å~ Dominance of SNP B Å~ Additive of SNP C\cr
26\tab	d1d2d3\tab	Dominance of SNP A Å~ Dominance of SNP B Å~ Dominance of SNP C\cr
}

The singificant SNP combinations with with epistatic effector, correclation figures and full data set are exported into the log file in the summary. 

}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{

library(Eip);

options  <- list(
	type        = "auto", 
	correction  = c("FDR, Holm"),
	ignoreRow  = NULL, 
	bHead      = TRUE);
	
FM.snp3_step( "data.csv", options=options );

}

\keyword{functions}
