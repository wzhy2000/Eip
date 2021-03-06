# Eip

Epistatic Interaction Package for R (Version 1.0)

The package has been moved from http://statgen.psu.edu/software/eip.html to here.

## Abstract

The Epistatic Interaction Package (Eip) is developed to estimate the epistatic effector in a candidate-gene case-control data set based on the Tian's model [1, 2]. In general, the case-control data has two groups, one group includes m cases who display a disease, and another group includes n control with no disease. All the cases and controls are genotyped for some SNPs. This package can calculate the p-value for all the combinations of 2 SNPs and 3 SNPs of case-control data. Based on the p-values results, some significant SNP combinations with epistatic effector are selected to export in the summary. For more details, please refer to the document.

[1] Wang, Zhong, Tian Liu, Zhenwu Lin, John Hegarty, Walter A. Koltun, and Rongling Wu. "A general model for multilocus epistatic interactions in case-control studies." PLoS One 5, no. 8 (2010): e11384. (http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0011384)

[2] Liu, Tian, A. Thalamuthu, J. J. Liu, C. Chen, Zhong Wang, and Rongling Wu. "Asymptotic distribution for epistatic tests in case–control studies." Genomics 98, no. 2 (2011): 145-151.


## Installation

The Eip package does not depend on any package, so it is very easy to install. After you download this package file, please type the following command or click the menu item "Install packages from local zip files".

Please use the install_github to install this package directly from GitHub. 

Windows/Linux/Mac OS:

>library(devtools);

>install_github("wzhy2000/Eip/Eip")

Before the package is used in R, the package importation is necessary by the following command:
> library(Eip)

After it is loaded, all functions within Funmap will be readily available to the user.

## Document

https://github.com/wzhy2000/Eip/blob/master/eip-intro.pdf


## Sample Script

```
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
```
## Sample Figures

![Image of 2 SNP test ](https://github.com/wzhy2000/Eip/blob/master/img/correlation-sample-snp2.jpg)

![Image of 3 SNP test ](https://github.com/wzhy2000/Eip/blob/master/img/correlation-sample-snp3.jpg)


