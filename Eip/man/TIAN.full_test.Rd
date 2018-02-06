\name{TIAN.full_test}
\alias{TIAN.full_test}

\title{Main function}

\description{ 
  Perform all analyses for the real data as a main function.
} 

\usage{ TIAN.full_test(sCsvFile, options=NULL, output_file=NULL, model=NULL) }

\arguments{
\item{sCsvFile}{CSV file as a data set.}
\item{options}{ parameters for the data file, see \code{\link{options}}.}
\item{output_file}{file name for the summary report. }
\item{model}{1,2,3 are available that indicate the interaction of 1 SNP, 2 SNPs, 3 SNPs.}
}

\value{
A list that includes the following results.

\item{ret1}{The result without correction for 1 SNP analysis.}
\item{ret_b1}{The results with Bonferroni correction for 1 SNP analysis.}
\item{ret_h1}{The results with Holm correction for 1 SNP analysis.}
\item{ret_f1}{The results with FDR correction for 1 SNP analysis.}

\item{ret2}{The results without correction for 2 SNP analysis.}
\item{ret_b2}{The results with Bonferroni correction for 2 SNP analysis.}
\item{ret_h2}{The results with Holm correction for 2 SNP analysis.}
\item{ret_f2}{The results with FDR correction for 2 SNP analysis.}

\item{ret3}{The results without correction for 3 SNP analysis.}
\item{ret_b3}{The results with Bonferroni correction for 3 SNP analysis.}
\item{ret_h3}{The results with Holm correction for 3 SNP analysis.}
\item{ret_f3}{The results with FDR correction for 3 SNP analysis.}

}

\details{

\strong{TIAN.full_test} is main function of the Eip package, it 
encapsulates 3 analysis procedures, including:

\bold{1) 1 SNP analysis.}

If parameter model is undefined or includes number 1, 1 SNP analysis is performed 
to get addititive and dominance effector by chi-square test. No figures are exported.


\bold{2) 2 SNPs analysis.}

If parameter model is undefined or includes number 2, 2 SNP analysis is performed 
to get the following effector by chi-square test. 

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

a). Correlation figures for all epistatic effectors are exported into the PDF files.\cr
b). Full data set are exported into the CSV files.


\bold{3) 3 SNPs analysis.}


If parameter model is undefined or includes number 3, 3 SNP analysis is performed 
to get the following effector by chi-square test. 

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

sys$set_value("legend.peak" , 0.000001);		

sCsvFile <- "test.csv";
options  <- list(
	type       = "auto", 
	case_cols  = "Disease",
	case_values= c("CD"),
	ctrl_values= c("Healthy"),
	snp_cols   = c("X1007FS","R702W","G908R"),
	snp_labels = c("1007FS","R702W","G908R"),
	correction = c("none"),
	ignoreRow  = NULL, 
	bHead      = TRUE );

r <- TIAN.full_test( sCsvFile, options=options, output="grp1", model=2 );

r$ret2;

}

\keyword{functions}
