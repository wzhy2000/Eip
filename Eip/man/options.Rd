\name{options}
\alias{options}

\title{Options}

\description{ 
   Parameters for the data file.
} 

\usage{ options  <- list(
		arguments   = "value");
 }

\arguments{
\item{type}{The format of data set. \strong{TA}, \strong{TB}, \strong{TC} , \strong{TD} and auto recognition(\strong{auto}) are available, see the following section for details.\cr
\strong{(Default:auto)}.}
\item{correction}{The correction method of multiple tests. \strong{bonf}, \strong{holm}, \strong{fdr} or combination of these methods are available.
bonf, holm and fdr indicate Bonferroni, Holm and FDR methods respectively.\cr
\strong{(Default:all)}.}
\item{case_col}{This optional argument shows which column provides the case-control information. Default value is the first colum provides case-control values.\cr
\strong{(Default:1)}.}
\item{case_values}{The case_values indicates which items in case-control column will be regarded as case value.\cr
\strong{(defualt 1)}.}
\item{ctrl_values}{The ctrl_values indicates which items in case-control column will be regarded as control value.\cr
\strong{(defualt 0)}.}
\item{snp_cols}{The snp_cols indicated which columns provide the SNP information. The analysis starts from the column which is next to case_col and stops at the last column if the value is unspecified.\cr
\strong{(Default:NA)}.}
\item{snp_labels}{The snp_labels gives the snp labels shown in the log files and figures. The labels will be assigned from the column names of SNP if the value is unspecified.\cr
\strong{(Default:NA)}.}
\item{filters}{This optional argument defines some filter ceriations.The individual data will be removed if some columns meet these ceriations.\cr
\strong{(Default:NA)}.}
\item{ignoreRow}{This optional argument indicates whether the n row below the head will be removed.\cr
\strong{(Default:NULL)}.}
\item{bHead}{This optional argument indicates whether the CSV includes a data head.\cr
\strong{(Default:TRUE)}.}
}

\value{
}

\details{

The goal of the Eip package is to find the significant SNP combination with epistatic effector from the candidate-gene case-control data set. The data set should be a CSV file, and at least includes one case-control column and some SNP columns. Now the current package supports 4 SNP formats. Hereby, we describle 4 types of data format in this guide.


\bold{Type A(TA)}

One SNP column includes one allele. In following table three SNP columns are listed with one case-control column.

\preformatted{
Individual ID, \dots, Case, SNP 1, SNP 2, SNP 3, \dots\dots
1            , \dots,    1,    GG,    CC,    CC, \dots\dots
2            , \dots,    1,    AA,    CT,    CA, \dots\dots
3            , \dots,    0,    AG,    CC,    CC, \dots\dots
\dots\dots
}

\bold{Type B(TB)}

One SNP data is splited into three columns that includs two homczygous alleles and one heterozygote allele. Above data set can be converted into the following table.

\preformatted{
Individual ID, \dots, Case, SNP 1, SNP 1-1, SNP 1-2, SNP 2, SNP 2-1, SNP 2-2,\dots\dots
1,             \dots,     1,   GG,        ,        ,    CC,        ,        ,\dots\dots
2,             \dots,     1,     ,        ,      AA,      ,      CT,        ,\dots\dots
3,             \dots,     0,     ,      AG,        ,    CC,        ,        ,\dots\dots
\dots\dots
}

\bold{Type C(TC)}

One SNP data is splited into two columns that every column includes one nucleotide. The following table is still same data set as the above.

\preformatted{
Individual ID, \dots, Case, SNP 1, SNP 1-1, SNP 2, SNP 2-1, \dots\dots
1,             \dots,    1,     G,       G,     C,       C, \dots\dots
2,             \dots,    1,     A,       A,     C,       T, \dots\dots
3,             \dots,    0,     A,       G,     T,       T, \dots\dots
\dots\dots
}

\bold{Type D(TD)}

The type D is kind of same as the type A. Every column is an allele, but in Type D, the allele is mapped to 0, 1, 2, which means qq, Qq, QQ.

\preformatted{
Individual ID, \dots, Case, SNP 1, SNP 2, SNP 3, SNP 4, \dots\dots
1,             \dots,     1,    1,     2,     0,     1, \dots\dots
2,             \dots,     1,    0,     1,     0,     1, \dots\dots
3,             \dots,     0,    1,     2,     1,     0, \dots\dots
\dots\dots
}


}

\author{Zhong Wang, \email{zzw2@psu.edu} }

\examples{
library(Eip);

options  <- list(
	type        = "auto", 
	correction  = c("FDR, Holm"),
	case_col    = "Disease",
	case_values = c("CD"),
	ctrl_values = c("Healthy"),
	snp_cols    = c("X1007S","X380Mut","Arg30Gln","L503F"),
	snp_labels  = c("1007S","X380Mut","Arg30Gln","L503F"),
	filters     = c("X1007S==\"TT\""),
	ignoreRow   = NULL, 
	bHead       = TRUE );		

TIAN.full_test( sCsvFile, options );

}

\keyword{functions}
