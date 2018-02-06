# Eip: Epistatic Interaction Package for R (Version 1.0)

## Description:
The Epistatic Interaction Package (Eip) is developed to estimate the epistatic effector in a candidate-gene case-control data set based on the Tian's model [1, 2]. In general, the case-control data has two groups, one group includes m cases who display a disease, and another group includes n control with no disease. All the cases and controls are genotyped for some SNPs. This package can calculate the p-value for all the combinations of 2 SNPs and 3 SNPs of case-control data. Based on the p-values results, some significant SNP combinations with epistatic effector are selected to export in the summary. For more details, please refer to the document.

[1] Tian Liu, A,Thalamuthu, J.J.Liu, C.Chen,Yao Li,and Rongling Wu. A Model for testing epistatic interactions of complex diseases in Case-Control studies.

[2] Zhong Wang, Zhenwu Lin, Arthur Berg, John Hegarty, Walter A. Koltun and Rongling Wu. A general model for multilocus epistatic interactions in Case-Control studies.

## Install
The Eip package does not depend on any package, so it is very easy to install. After you download this package file, please type the following command or click the menu item "Install packages from local zip files".

Windows OS:
>install.packages("x:/fullpath/Eip_1.0-1.zip", repos=NULL)

Linux/Mac OS:
>install.packages("/fullpath/ Eip_1.0-1.tar.gz", repos=NULL)

Before the package is used in R, the package importation is necessary by the following command:
> library(Eip)

After it is loaded, all functions within Funmap will be readily available to the user.

## Document


## Sample Script:

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
## Sample Figures:


