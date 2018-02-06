###############################################################
# References:
#1. A Model for Testing Epistatic Interactions of Complex Diseases
#   in Case-Control Studies.
#
#   Tian Liu, etc.
#
#2. A General Model for Multilocus Epistatic Interactions in
#   Case-Control Studies
#
#   Zhong Wang, etc.
###############################################################

N2_A1	  <- 1;
N2_A2	  <- 2;
N2_D1	  <- 3;
N2_D2	  <- 4;
N2_A1A2	  <- 5;
N2_A1D2	  <- 6;
N2_D1A2	  <- 7;
N2_D1D2	  <- 8;

N3_A1	  <- 1;
N3_A2	  <- 2;
N3_A3	  <- 3;
N3_D1	  <- 4;
N3_D2	  <- 5;
N3_D3	  <- 6;
N3_A1A2	  <- 7;
N3_A2A3	  <- 8;
N3_A1A3	  <- 9;
N3_D1D2	  <- 10;
N3_D2D3	  <- 11;
N3_D1D3	  <- 12;
N3_A1D2	  <- 13;
N3_A1D3	  <- 14;
N3_D1A2	  <- 15;
N3_D1A3	  <- 16;
N3_A2D3	  <- 17;
N3_D2A3	  <- 18;
N3_A1A2A3 <- 19;
N3_A1A2D3 <- 20;
N3_A1D2A3 <- 21;
N3_A1D2D3 <- 22;
N3_D1A2A3 <- 23;
N3_D1A2D3 <- 24;
N3_D1D2A3 <- 25;
N3_D1D2D3 <- 26;

EFACTOR1 <- c(	"AA:Aa:aa ",
				"Additive ",
				"Dominance");

EFACTOR2 <- c(	"a1  ",
				"a2  ",
				"d1  ",
				"d2  ",
				"a1a2" ,
				"a1d2" ,
				"d1a2" ,
				"d1d2" );

EFACTOR3 <- c(	"a1   ",
				"a2   ",
				"a3   ",
				"d1   ",
				"d2   ",
				"d3   ",
				"a1a2 " ,
				"a2a3 " ,
				"a1a3 " ,
				"d1d2 " ,
				"d2d3 " ,
				"d1d3 " ,
				"a1d2 " ,
				"a1d3 " ,
				"a2d3 " ,
				"d1a2 " ,
				"d2a3 " ,
				"d1a3 " ,
				"a*a*a" ,
				"a*a*d" ,
				"a*d*a" ,
				"a*d*d" ,
				"d*a*a" ,
				"d*a*d" ,
				"d*d*a" ,
				"d*d*d");

sys<-NULL;

.First.lib<-function(lib, pkg)
{
	cat(lib, "Thank you!\r\n");

	#var:  snp2_pv_table
	#format: SampleSize, E_factor, a1,a2,0,par_a, par_b,0, var(estimate-simulation), mean(estimate-simulation);
	#load("cc_snp2_l100000_s.RData");
	data(cc_snp2_l100000_s);

	#var:  snp3_pv_table
	#format: SampleSize, E_factor, a1,a2,a3,par_a, par_b,0, var(estimate-simulation), mean(estimate-simulation);
	#load("cc_snp3_l100000_s.RData");
	data(cc_snp3_l100000_s);

	sys <<- sys_mgr.base("TIAN");
	sys$reset();
}

#--------------------------------------------------------------
# Table Structur TA
#       case/ctrl  snp1 snp2 snp3
#       1          AG   CT   AA
#       0          GG   CC   AA
#       1          -    CC   AG
#       ......
#
#
# Table Structur TB
#       case/ctrl  snp1         snp2        snp3
#       1             , AG,   ,   , CT,   , AA,   ,
#       0             ,   , GG, CC,   ,   , AA,   ,
#       1           --, --, --, CC,   ,   ,   ,   ,AG,
#       ......
#
# Table Structur TC
#       case/ctrl  snp1a,b snp2a,b snp3a,b
#       1              A,G,    C,T,    A,A
#       0              G,G,    C,C,    A,A,
#       1              -,-,    C,C,    A,G,
#       ......
#
# Table Structur TD
#       case/ctrl  snp1 snp2 snp3
#       1          1    1    0
#       0          2    0    0
#       1          -    0    1
#
#
# --- or---
# data.frame structure DE
#       case/ctrl  snp1 snp2 snp3
#       1          AG   CT   AA
#       0          GG   CC   AA
#
# data.frame structure DF
#       case/ctrl  snp1 snp2 snp3
#       1          1    1    0
#       0          2    0    0
#       1          -    0    1
#--------------------------------------------------------------



#--------------------------------------------------------------
# Public:TIAN.readtable
#    for TA: sSnpCols=c(i1,i2,i3)
#    for TB: sSnpCols=c(i11,i12,i13; i21,i22,i23; i31,i32,i33)
#    for TC: sSnpCols=c(i11,i12;i21,i22;i31,i32)
#    for TD: sSnpCols=c(i1,i2,i3)
#--------------------------------------------------------------
TIAN.readtable<-function( sCsvFile, options=NULL )
{
	#options<-list(
	#	type       = "TA",
	#	case_col   = "Disease",
	#	case_values= c("CD"),
	#	ctrl_values= c("Healthy"),
	#	snp_cols   = c(3,4,5,6,7,8),
	#	#snp_cols   = c("C.Ins","T.mut","X380Mut","Arg30Gln","L503F","R381Q"),
	#	filters    = c("C.Ins==\"TT\""),
	#	ignoreRow  = NULL,
	#	bHead      = TRUE );

	if (is.null(options) || is.null(options$bHead))
		bHead = TRUE
	else
		bHead = options$bHead;

	tb <- read.csv(sCsvFile, head = bHead);
	if ( is.null(options) )
		options <- list(
					type       = "auto",
					case_col   = 1,
					case_value = c(1),
					ctrl_value = c(0),
					snp_cols   = c(2:length(tb[1,])),
					ignoreRow  = NULL,
					bHead      = TRUE)

	if (!is.null(options$ignoreRow))
		tb <- tb[-(1:options$ignoreRow),];

	#check the case/ctrl column.
	case_col <- options$case_col;
	if (!is.numeric(options$case_col))
	{
		case_cols <- which(colnames(tb)==options$case_col);
		if (length(case_cols)!=1)
		{
			stop(paste( length(case_cols), "column(s) to be specified as case/ctrl data.", sep=""));
		}
		else
			case_col <- case_cols[1];
	}

	#check the snp columns.
	if (is.null(options$snp_cols))
	{
		#stop("Snp columns is not specified.");
		start <- options$case_col+1;
		end <- length( tb[1,] );
		snp_cols <- c(start:end);
	}
	else
		snp_cols <- options$snp_cols;

	if (!is.numeric(snp_cols))
	{
		snp_cols <- rep(-1, length(options$snp_cols) );
		for (i in 1:length(options$snp_cols))
		snp_cols[i] <- which(colnames(tb)==options$snp_cols[i]);
		if ( length( which(snp_cols < -1) ) >0 )
			stop(paste( "No column to be matched your request as snp data.", sep=""));
	}

	if (options$bHead==TRUE)
		sSnpNames <- colnames(tb)[snp_cols]
	else
		sSnpNames <- paste("SNP", c(1:length(snp_cols)),sep="");
	if (!is.null( options ) &&
	    !is.null( options$snp_labels ) &&
	    length(options$snp_labels)==length(snp_cols) )
	    sSnpNames <- options$snp_labels;

	#check the filetype
	#TA Snpcol1: GG;AG;
	#TB Snpcol1: GG,,;,AG,;
	#TB Snpcol1: G,G;A,G;
	#TA Snpcol1: 0;1;2;
	snpdata1 <- tb[,snp_cols[1]];
	snpdata2 <- try( tb[,snp_cols[1]+1], FALSE) ;
	if (class(snpdata2)=="try-error")
		snpdata2 <-	 snpdata1;
	snpdata3 <- try( tb[,snp_cols[1]+2], FALSE) ;
	if (class(snpdata3)=="try-error")
		snpdata3 <- snpdata1;

	samps <- which( !is.na(snpdata1) & snpdata1!="" &
					!is.na(snpdata2) & snpdata2!="" &
					!is.na(snpdata3) & snpdata3!="" );
	if (length(samps)<0)
		stop(paste("No data in the column ", snp_cols[1], ".", sep=""));

	if (nchar(as.character( snpdata1[samps[1]]))==2)
	{
		snpdata <- paste(snpdata1,snpdata2, snpdata3, sep="");
		if ( nchar( as.character( snpdata[1])) > 2 )
			ftype="TA"
		else
			ftype="TB";
	}
	else if (nchar(as.character( snpdata1[samps[1]]))==1)
	{
		if (is.numeric(snpdata1[samps[1]]) )
			ftype <- "TD"
		else
			ftype <- "TC";
	}
	else
		stop( paste( "Unknow snp data found, data=(", snpdata1[samps[1]], ")", sep=""));

	if (options$type=="auto" || is.null(options$type))
		ftype=ftype
	else if (ftype!=options$type)
	{
		browser();
		stop( paste( "Wrong file format is specified.", sep=""));
	}

	#create data object
	dat<-list(
		name    = "TIAN.dat",
		DataFile= sCsvFile,
		FileType= ftype,
		CaseCol = case_col,
		SnpCols = snp_cols,
		SnpNames= sSnpNames,
		CaseList= NULL,
		SnpList = NULL,
		options = options  );
	class(df) <- "TIAN.dat";

	#filt the data firstly,
	if (!is.null(options$filters))
	{
		for (k in 1:length(options$filters))
		{
			cond_str <- paste("which(tb$", options$filters[k], ")", sep="");
			remove <- try(  eval(parse( text=cond_str ) ), TRUE);
			if (class(remove)=="try-error")
			{
				warnings(paste("Failed to execute the filter ",options$filters[k], ".", sep="" ));
				next;
			}

			if (length(remove)>0)
				tb <- tb[-remove,];
		}
	}

	tb1 <- c();
	#select the case data
	if (!is.null(options$case_values))
	{
		for (k in 1:length(options$case_values))
		{
			select <- which(tb[, case_col]==options$case_values[k]);
			if (length(select)>0)
				tb1 <- rbind( tb1, tb[select,] );
		}
	}
	else
	{
		select <- which(tb[,case_col]==1);
		if (length( select ) > 0 )
			tb1 <- tb[select,];
	}

	if (length(tb1)==0)
		stop("No data for case value.")

	tb0 <- c();
	#select the control data
	if (!is.null(options$ctrl_values))
	{
		for (k in 1:length(options$ctrl_values))
		{
			select <- which(tb[, case_col]==options$ctrl_values[k]);
			if (length(select)>0)
				tb0 <- rbind( tb0, tb[select,] );
		}
	}
	else
	{
		select <- which(tb[,case_col]==0);
		if (length(select)>0)
			tb0 <- tb[select,];
	}

	if (length(tb0)==0)
		stop("No data for control value.")

	tb1[,case_col]<-1;
	tb0[,case_col]<-0;
	tb <- rbind(tb1, tb0);
	dat$CaseList <- tb[,case_col];
	n<-length(dat$CaseList);
	dat$SnpList <- list();

	if (ftype=="TA")
	{
		for( i in 1:length(snp_cols))
		{
			snplist <- sapply( tb[,snp_cols[i]], trim_space );
			snplist[snplist==""] <- NA;
			if ( !check_file( snplist ) )
				stop("The genotype should be less than three.");

			dat$SnpList[[ i ]] <- snplist;
		};
	}

	if (ftype=="TB")
	{
		for( i in 1:length(snp_cols))
		{
			snplist <- sapply( paste( tb[,snp_cols[i]], tb[,snp_cols[i]+1], tb[,snp_cols[i]+2], sep="") , trim_space );
			snplist[snplist==""] <- NA;
			if ( !check_file( snplist ) )
				stop("The genotype should be less than three.");

			dat$SnpList[[ i ]] <- snplist;
		}
	}

	if (ftype=="TC")
	{
		for( i in 1:length(snp_cols))
		{
			snplist <- sapply( paste( tb[,snp_cols[i]], tb[,snp_cols[i]+1], sep=""), trim_space);
			snplist[snplist==""] <- NA;
			if ( !check_file( snplist ) )
				stop("The genotype should be less than three.");

			dat$SnpList[[ i ]] <- snplist;
		}
	}

	if (ftype=="TD")
	{
		for( i in 1:length(snp_cols))
		{
			snp1 <- as.character(tb[,snp_cols[i]]);
			snp1[snp1==0]<-"AA";
			snp1[snp1==1]<-"AB";
			snp1[snp1==2]<-"BB";

			if ( !check_file( snp1 ) )
				stop("The genotype should be less than three.");

			dat$SnpList[[ i ]] <-snp1;
		}
	}

	class(dat)<-dat$name;
	return(dat);
}

check_file<-function( snplist )
{
	ls <- union(snplist, snplist);
	if (NA %in% ls)
		ls <- setdiff(ls, c(NA));

	if (length(ls)>3)
		return(FALSE);

	gens <- c();
	lens <- c();
	for (i in 1:length(ls))
	{
		if ( nchar(as.character(ls[i]) )>2 )
			return(FALSE);
		lens <- union(lens, nchar(as.character(ls[i]) ) );
		gens<- union(gens, substr(ls[i], 1,1) );
		gens<- union(gens, substr(ls[i], 2,2) );
	}

	if (length(gens)>2)
		return(FALSE);

	return(TRUE);
}

TIAN.summary_dat<-function(  dat )
{
	if ( class(dat) != "TIAN.dat" )
	{
		sErrMsg<- "Error: Not a data set for Epistatic Interactions Model.";
		stop( sErrMsg );
	}

	#df<-list(
	#	name    = "TIAN.dat",
	#	DataFile= sCsvFile,
	#	FileType= type,
	#	CaseCol = nCaseCol,
	#	SnpCols = sSnpCols,
	#	SnpNames= sSnpNames,
	#	CaseList= NULL,
	#	SnpList = NULL  );

	nCaseCnt <- sum(dat$CaseList);
	nCtrlCnt <- length(dat$CaseList)-sum(dat$CaseList);

	stru <- sprintf("------------------------------------\n");
	str0 <- sprintf("%10s: %s\n", 	"Data File",	dat$DataFile );
	str1 <- sprintf("%10s: %s\n", 	"Date", 		Sys.time() );
	str2 <- sprintf("%10s: %s\n", 	"File type", 	dat$FileType );
	str3 <- sprintf("%10s: %-10.0f\n", "SNPs Count", 	length( dat$SnpList ) );
	str4 <- sprintf("%10s: %-10.0f\n", "Cases", 		nCaseCnt );
	str5 <- sprintf("%10s: %-10.0f\n", "Controls", 	nCtrlCnt );
	strd <- sprintf("------------------------------------\n");

	str<-paste(stru, str0, str1, str2, str3, str4, str5, strd, sep="" );

	str<-paste(str, "\nSNPs List:", sep="");
	str<-paste(str, "\nName\tPos.\tAA(case:ctrl)\tAa\taa\n", sep="");

	for (i in  1:length(dat$SnpList))
	{
		r <-pFuncCalcStandard( dat$CaseList, dat$SnpList[[i]], 1, list(c(0), c(1),c(2))  );
		genAA<-paste(r$GenA, r$GenA, sep="");
		genAa<-paste(r$GenA, r$Gena, sep="");
		genaa<-paste(r$Gena, r$Gena, sep="");

		##SNP1\t  2\t GG(89:90)\t  AG(89:90)\t  AA(89:90)\t
		str <- paste(str, 	dat$SnpNames[[i]],"\t",
		 					dat$SnpCols[[i]], "\t",
		 					genAA, "(", r$M[1,1], ":", r$M[2,1], ")\t",
		 					genAa, "(", r$M[1,2], ":", r$M[2,2], ")\t",
		 					genaa, "(", r$M[1,3], ":", r$M[2,3], ")\n", sep="" );
	}

	return(str);
}

#--------------------------------------------------------------
# Public:TIAN.snp1
#
#--------------------------------------------------------------
TIAN.snp1<-function(df, correct=NULL)
{
	ret<-list(
		type     = 1,
		name     = "TIAN.ret.1",
		dataFile = df$DataFile,
		Exception = "",
		correction = "",
		nCaseCnt = sum(df$CaseList),
		nCtrlCnt = length(df$CaseList)-sum(df$CaseList),
		SnpNames = df$SnpNames,
		model    = list() );

	nSnpCnt<-length(df$SnpCols);
	if (nSnpCnt<1)
	{
		ret$Exception<-"Snp count is less than 1.";
		return(ret);
	}

	class(ret)<- "TIAN.ret.1";

	sys$task_start("\nExecute the model(1 SNP) task, data file:", df$DataFile,", SNP count:", length(df$SnpCols),".\n");

	nTaskIdx <- 1;
	nTaskCnt <- 0;
	for (i in 1:nSnpCnt)
		nTaskCnt <- nTaskCnt+1;

	for (i in 1:nSnpCnt)
	{
		snp1_list <- df$SnpList[[i]];
		ret1 <- pInteractionSnps1( df$CaseList, snp1_list )
		ret$model <- rbind( ret$model, c(i, ret1 ));

		sys$task_elapsed( nTaskIdx/nTaskCnt, "Task Index:", nTaskIdx, "/",nTaskCnt, ", ", "$SYS_PROMPT$", "\n");
		nTaskIdx <- nTaskIdx + 1;
	}

	sys$task_stop("The model(1 SNP) task is done.\n\n");

	if (!is.null(correct))
	{
		if (correct=="bonf")
			return ( TIAN.adjust_by_bonferroni( ret ) )
		else if (correct=="holm")
			return ( TIAN.adjust_by_holm( ret ) )
		else if (correct=="fdr")
			return ( TIAN.adjust_by_fdr( ret ) )
   	}

	return(ret);
}

#--------------------------------------------------------------
# Public:TIAN.snp2
#
#--------------------------------------------------------------
TIAN.snp2<-function(df, correct=NULL)
{
	ret<-list(
		type     = 2,
		name     = "TIAN.ret.2",
		Exception = "",
		correction = "",
		dataFile = df$DataFile,
		nCaseCnt = length(which(df$CaseList==1)),
		nCtrlCnt = length(df$CaseList) - length(which(df$CaseList==1)),
		SnpNames = df$SnpNames,
		model    = c() );
	class(ret)<- "TIAN.ret.2";
	nSnpCnt<-length(df$SnpCols);

	if (nSnpCnt<2)
	{
		ret$Exception<-"Snp count is less than 2.";
		return(ret);
	}

	sys$task_start("\nExecute the model(2 SNPs) task, data file:", df$DataFile,", SNP count:", length(df$SnpCols),".\n");


	nTaskIdx <- 0;
	nTaskCnt <- 0;
	for (i in 1:(nSnpCnt-1))
	for (j in (i+1):nSnpCnt)
		nTaskCnt <- nTaskCnt+1;


	for (i in 1:(nSnpCnt-1))
	{
		for (j in (i+1):nSnpCnt)
		{
			snp1_list <- df$SnpList[[i]];
			snp2_list <- df$SnpList[[j]];

			ret2 <- pInteractionSnps2( df$CaseList, snp1_list, snp2_list,i, j );
			ret$model<- rbind( ret$model, c(ret2));

			nTaskIdx <- nTaskIdx + 1;
		}

		sys$task_elapsed( nTaskIdx/nTaskCnt, "Task Index:", nTaskIdx, "/", nTaskCnt, ",", "$SYS_PROMPT$", "\n");
	}

	sys$task_stop("The model(2 SNPs) task is done.\n\n");

	if (!is.null(correct))
	{
		if (correct=="bonf")
			return ( TIAN.adjust_by_bonferroni( ret ) )
		else if (correct=="holm")
			return ( TIAN.adjust_by_holm( ret ) )
		else if (correct=="fdr")
			return ( TIAN.adjust_by_fdr( ret ) );
   	}

	return(ret);
}

#--------------------------------------------------------------
# Public:TIAN.snp3
#
#--------------------------------------------------------------
TIAN.snp3<-function(df, correct=NULL)
{
	ret<-list(
		type     = 3,
		name     = "TIAN.ret.3",
		Exception = "",
		correction = "",
		dataFile = df$DataFile,
		nCaseCnt = length(which(df$CaseList==1)),
		nCtrlCnt = length(df$CaseList) - length(which(df$CaseList==1)),
		SnpNames = df$SnpNames,
		model    = c() )
	class(ret)<- "TIAN.ret.3";

	nSnpCnt<-length(df$SnpCols);
	if (nSnpCnt<3)
	{
		ret$Exception<-"Snp count is less than 3.";
		return(ret);
	}

	sys$task_start("\nExecute the model(3 SNPs) task, data file:", df$DataFile,", SNP count:", length(df$SnpCols),".\n");


	nTaskIdx <- 0;
	nTaskCnt <- 0;
	for (i in 1:(nSnpCnt-2))
	for (j in (i+1):(nSnpCnt-1))
	for (k in (j+1):(nSnpCnt-0))
		nTaskCnt <- nTaskCnt+1;

	for (i in 1:(nSnpCnt-2))
	for (j in (i+1):(nSnpCnt-1))
	{
		for (k in (j+1):(nSnpCnt-0))
		{
			if(i==j || i==k || k==j)
				next;

			snp1_list <- df$SnpList[[i]];
			snp2_list <- df$SnpList[[j]];
			snp3_list <- df$SnpList[[k]];

			ret3 <- pInteractionSnps3( df$CaseList, snp1_list, snp2_list, snp3_list, i, j, k );
			ret$model <- rbind( ret$model, c(ret3 ));

			nTaskIdx <- nTaskIdx + 1;
		}

		sys$task_elapsed( nTaskIdx/nTaskCnt, "Task Index:", nTaskIdx, "/", nTaskCnt, ",", "$SYS_PROMPT$", "\n");
	}

	sys$task_stop("The model(3 SNPs) task is done.\n\n");

	if (!is.null(correct))
	{
		if (correct=="bonf")
			return ( TIAN.adjust_by_bonferroni( ret ) )
		else if (correct=="holm")
			return ( TIAN.adjust_by_holm( ret ) )
		else if (correct=="fdr")
			return ( TIAN.adjust_by_fdr( ret ) );
   	}

	return(ret);
}

#--------------------------------------------------------------
# Public:TIAN.summary1
#
#------
# 1 SNPs interactions
#
#         AA:Aa:aa AA:aa AA+aa:2Aa
# snpx1   0.10     0.09  0.06
#--------------------------------------------------------------
TIAN.summary1<-function(res, outputFile=NULL)
{
	if (res$type != 1)
	{
		sErrMsg<- "Error: Not a result for 1 SNPs interaction.";
		stop( sErrMsg );
	}

	if ( res$Exception!= "")
	{
		sErrMsg<- paste("Exception: ", res$Exception, sep="");
		stop( sErrMsg );
	}

	str = pSummaryHeader(res);

	sig_marker <- c(1e-8, 1e-6, 1e-4, 1e-2, 0.05);
	for (k in 1:length(sig_marker))
	{
		pvs.len = (length(res$model[1,])-5+1);
		pvs <- matrix( res$model[, 5:length(res$model[1,]) ], ncol=pvs.len);
		#pvs <- res$model[, 5:length(res$model[1,]) ];
		if (k== 1 )
			sigs <- which( pvs <= sig_marker[k] )
		else
			sigs <- which( (pvs <= sig_marker[k]) & (pvs > sig_marker[k-1]) );

		if (length(sigs)>0)
		{
			sigs_col <- (sigs-1) %/% length(pvs[,1]) +1 ;
			sigs_row <- (sigs-1) %% length(pvs[,1]) +1;

			pvalue <- c();
			chi2_v <- c();
			for (i in 1:length(sigs_row))
			{
				pvalue <- c(pvalue, unlist(res$model[sigs_row[i],sigs_col[i]+4]) );
				chi2_v <- c(chi2_v, unlist(res$model[sigs_row[i],sigs_col[i]+1]) );
			}

			sigs_list <- cbind( snp 	= res$model[sigs_row,1],
							    efactor = sigs_col,
							    pvalue	= pvalue,
							    chi2	= chi2_v);
			str_sigs <-  pSummaryInteraction( res$SnpNames, sigs_list );

			str<-paste(str, "p-Value <= ", sig_marker[k], ":\n", str_sigs, "\n", sep="");

		}
	}

	if ( !is.null(outputFile) )
	{
		correct <- ifelse( res$correction=="", "none", res$correction);
		csvFile1 <- paste(outputFile, ".s1.",correct,".pvalue.csv", sep="");
		TIAN.export_csv( res, csvFile1, FALSE  );
		csvFile2 <- paste(outputFile, ".s1.",correct, ".chi2.csv", sep="");
		TIAN.export_csv( res, csvFile2, TRUE  );

		str<- paste(str, "#1: Full p-value data is exported to the CSV file(", csvFile1,").\n", sep="" );
		str<- paste(str, "#2: Full chi2-value data is exported to the CSV file(", csvFile2,").\n\n", sep="" );
	}


	return(str);
}

#--------------------------------------------------------------
# Public:TIAN.summary2
#
#------
# Case:     xx
# Control:  xx
# SNP count:xx
#
# p-Value <=1e-8
#   snpx1, snpx2, a*d   x.xxxxx
#   snpx1, snpx2, a*d   x.xxxxx
#   snpx1, snpx2, a*d   x.xxxxx
#   snpx1, snpx2, a*d   x.xxxxx
#
# p-Value <=1e-6
#   snpx1, snpx2, a*d   x.xxxxx
#   snpx1, snpx2, a*d   x.xxxxx
#   snpx1, snpx2, a*d   x.xxxxx
#   snpx1, snpx2, a*d   x.xxxxx
#--------------------------------------------------------------
TIAN.summary2<-function(res, outputFile=NULL)
{
	if (res$type != 2)
	{
		sErrMsg<- "Error: Not a result for 2 SNPs interaction.";
		stop( sErrMsg );
	}

	if ( res$Exception!= "")
	{
		sErrMsg<- paste("Exception: ", res$Exception, sep="");
		stop( sErrMsg );
	}

	str = pSummaryHeader(res);

	sig_marker <- c(1e-8, 1e-6, 1e-4, 1e-2, 0.05);
	for (k in 1:length(sig_marker))
	{
		pvs.len = (length(res$model[1,])-11+1);
		pvs <- matrix( res$model[, 11:length(res$model[1,]) ], ncol=pvs.len);
		#pvs <- res$model[, 11:length(res$model[1,]) ];
		if (k== 1 )
			sigs <- which( pvs <= sig_marker[k] )
		else
			sigs <- which( (pvs <= sig_marker[k]) & (pvs > sig_marker[k-1]) );

		if (length(sigs)>0)
		{
			sigs_col <- (sigs-1) %/% length(pvs[,1]) +1 ;
			sigs_row <- (sigs-1) %% length(pvs[,1]) +1;

			pvalue <- c();
			chi2_v <- c();
			for (i in 1:length(sigs_row))
			{
				pvalue <- c(pvalue, unlist(res$model[sigs_row[i],sigs_col[i]+10]) );
				chi2_v <- c(chi2_v, unlist(res$model[sigs_row[i],sigs_col[i]+2]) );
			}

			sigs_list <- cbind(  snp1 	 = res$model[sigs_row,1],
								 snp2 	 = res$model[sigs_row,2],
								 efactor = sigs_col,
								 pvalue  = pvalue,
								 chi2	= chi2_v);

			sigs_list <- matrix(unlist(sigs_list), ncol=5 );
			sigs_list <- sigs_list[ order(sigs_list[,1], sigs_list[,2]),];
			sigs_list <- matrix( sigs_list, ncol=5  );
			str_sigs  <- pSummaryInteraction( res$SnpNames, sigs_list );
			str<-paste(str, "p-Value <= ", sig_marker[k], ":\n", str_sigs, "\n", sep="");
		}
	}

	if ( !is.null(outputFile) )
	{
		correct <- ifelse( res$correction=="", "none", res$correction);
		csvFile1 <- paste(outputFile, ".s2.",correct,".pvalue.csv", sep="");
		TIAN.export_csv( res, csvFile1, FALSE  );
		csvFile2 <- paste(outputFile, ".s2.",correct,".chi2.csv", sep="");
		TIAN.export_csv( res, csvFile2, TRUE  );
		str<- paste(str, "#1: Full p-value data is exported to the CSV file(", csvFile1,").\n", sep="" );
		str<- paste(str, "#2: Full chi2-value data is exported to the CSV file(", csvFile2,").\n\n", sep="" );

		pdfFile <- paste(outputFile, ".s2.",correct,".pdf", sep="");
		TIAN.draw_correlate2( res, pdfFile  );
		str <-	paste( str, "#3: The figure is ", pdfFile, "\n", sep="");
	}


	return(str);
}

#--------------------------------------------------------------
# Public:TIAN.summary3
#
#------
# Case:     xx
# Control:  xx
# SNP count:xx
#
# p-Value <=1e-8
#   snpx1, snpx2, snpx3, a*d   x.xxxxx
#   snpx1, snpx2, snpx3, a*d   x.xxxxx
#   snpx1, snpx2, snpx3, a*d   x.xxxxx
#   snpx1, snpx2, snpx3, a*d   x.xxxxx
#
# p-Value <=1e-6
#   snpx1, snpx2, snpx3, a*d   x.xxxxx
#   snpx1, snpx2, snpx3, a*d   x.xxxxx
#   snpx1, snpx2, snpx3, a*d   x.xxxxx
#   snpx1, snpx2, snpx3, a*d   x.xxxxx
#--------------------------------------------------------------
TIAN.summary3<-function(res, outputFile=NULL)
{
	if (res$type != 3 )
	{
		sErrMsg<- "Error: Not a result for 3 SNPs interaction.";
		stop( sErrMsg );
	}

	if ( res$Exception!= "")
	{
		sErrMsg<- paste("Exception: ", res$Exception, sep="");
		stop( sErrMsg );
	}

	str = pSummaryHeader(res);

	sig_marker <- c(1e-8, 1e-6, 1e-4, 1e-2, 0.05);

	for (k in 1:length(sig_marker))
	{
		pvs.len = (length(res$model[1,])-30+1);
		pvs <- matrix( res$model[, 30:length(res$model[1,]) ], ncol=pvs.len);
		if (k== 1 )
			sigs <- which( pvs <= sig_marker[k] )
		else
			sigs <- which( (pvs <= sig_marker[k]) & (pvs > sig_marker[k-1]) );

		if (length(sigs)>0)
		{
			sigs_col <- (sigs-1) %/% length(pvs[,1]) +1 ;
			sigs_row <- (sigs-1) %% length(pvs[,1]) +1;
			pvalue <- c();
			chi2_v <- c();
			for (i in 1:length(sigs_row))
			{
				pvalue <- c(pvalue, unlist(res$model[sigs_row[i],sigs_col[i]+3+26]) );
				chi2_v <- c(chi2_v, unlist(res$model[sigs_row[i],sigs_col[i]+3]) );
			}

			sigs_list <- cbind(  snp1 	 = res$model[sigs_row,1],
					     snp2 	 = res$model[sigs_row,2],
					     snp3 	 = res$model[sigs_row,3],
					     efactor	 = sigs_col,
					     pvalue	 = pvalue,
					     chi2	 = chi2_v);

			sigs_list <- matrix(unlist(sigs_list), ncol=6 );
			sigs_list <- sigs_list[ order( sigs_list[,1], sigs_list[,2], sigs_list[,3] ), ];
			sigs_list <- matrix( sigs_list, ncol=6  );
			str_sigs <-  pSummaryInteraction( res$SnpNames, sigs_list )
			str<-paste(str, "p-Value <= ", sig_marker[k], ":\n", str_sigs, "\n", sep="");
		}
	}


	if ( !is.null(outputFile) )
	{
		correct <- ifelse( res$correction=="", "none", res$correction);
		csvFile1 <- paste(outputFile, ".s3.", correct, ".pvalue.csv", sep="");
		TIAN.export_csv( res, csvFile1, FALSE  );
		csvFile2 <- paste(outputFile, ".s3.", correct, ".chi2.csv", sep="");
		TIAN.export_csv( res, csvFile2, TRUE  );

		str<- paste(str, "#1: Full p-value data is exported to the CSV file(", csvFile1,").\n", sep="" );
		str<- paste(str, "#2: Full chi2-value data is exported to the CSV file(", csvFile2,").\n\n", sep="" );

		pdfFile <- paste(outputFile, ".s3.", correct, ".pdf", sep="");
		TIAN.draw_correlate3( res, pdfFile  );
		str <-	paste( str, "#3: The figure is ", pdfFile, "\n", sep="");
	}

	return(str);
}

select_flank3<-function(fulldat, snpNames, sigs)
{
	sigs    <- unlist(sigs);
	snp1    <- sigs[1];
	snp2    <- sigs[2];
	snp3    <- sigs[3];
	efactor <- sigs[4];

	s1    <- which( fulldat[,1]==snp1);
	dat1  <- fulldat[ s1, c( 2,3,efactor + 3) ]
	title1<- paste( EFACTOR3[efactor],"(", snpNames[snp1], ",*,*)" , sep="");
	label1<- paste( snpNames[dat1[,1]], ",", snpNames[dat1[,2]], sep="");

	s2    <- which( fulldat[,2]==snp2);
	dat2  <- fulldat[ s2, c( 1,3,efactor + 3) ];
	title2<- paste( EFACTOR3[efactor],"(*,", snpNames[snp2], ",*)" , sep="");
	label2<- paste( snpNames[dat1[,1]], ",", snpNames[dat1[,2]], sep="");

	s3    <- which( fulldat[,3]==snp3);
	dat3  <- fulldat[ s3, c( 1,2,efactor + 3) ];
	title3<- paste( EFACTOR3[efactor],"(*,*,", snpNames[snp3], ")" , sep="");
	label3<- paste( snpNames[dat1[,1]], ",", snpNames[dat1[,2]], sep="");

	return (list(data   = list( dat1,  dat2,  dat3),
	             title  = list( title1,title2, title3),
	             label  = list( label1,label2,label3) ) );
}

#--------------------------------------------------------------
# Public:TIAN.export_csv
#
#--------------------------------------------------------------
TIAN.export_csv<-function( res, outputFile, bChi2=FALSE )
{
	if (res$type == 1 )
	{
		export_csv1(res, outputFile, bChi2);
	}
	else if (res$type == 2 )
	{
		export_csv2(res, outputFile, bChi2);
	}
	else if (res$type == 3 )
	{
		export_csv3(res, outputFile, bChi2);
	}
}

#--------------------------------------------------------------
# Public:export_csv1
#
#--------------------------------------------------------------
export_csv1<-function( res, outputFile, bChi2=FALSE )
{
	if (res$type != 1 )
	{
		sErrMsg<- "Error: Not a result for 1 SNPs interaction.";
		stop( sErrMsg );
	}

	if (!bChi2)
	{
		items <- data.frame( res$SnpNames[ unlist(res$model[,1]) ],
							 res$model[, c( 5:7 ) ] );
		colnames(items)<-  c(  "SNP1", paste("pv-", EFACTOR1, sep="") );
	}
	else
	{
		items <- data.frame( res$SnpNames[ unlist(res$model[,1]) ],
							 res$model[, c( 2:4 ) ] );
		colnames(items)<-  c(  "SNP1", paste("chi2-", EFACTOR1, sep="") );
	}

	write.table(as.matrix(items), file=outputFile, sep=",");

	return( read.csv(outputFile) );
}

#--------------------------------------------------------------
# Public:export_csv2
#
#--------------------------------------------------------------
export_csv2<-function( res, outputFile, bChi2=FALSE  )
{
	if (res$type != 2 )
	{
		sErrMsg<- "Error: Not a result for 2 SNPs interaction.";
		stop( sErrMsg );
	}

	if (!bChi2)
	{
		items <- data.frame( res$SnpNames[ unlist(res$model[,1]) ],
							 res$SnpNames[ unlist(res$model[,2]) ],
							 res$model[, c( 11:18 ) ] );
		colnames(items)<-  c(  "SNP1", "SNP2", paste("pv-", EFACTOR2, sep="") );
	}
	else
	{
		items <- data.frame( res$SnpNames[ unlist(res$model[,1]) ],
							 res$SnpNames[ unlist(res$model[,2]) ],
							 res$model[, c( 3:10 ) ] );
		colnames(items)<-  c(  "SNP1", "SNP2", paste("chi2-", EFACTOR2, sep="") );
	}

	write.table(as.matrix(items), file=outputFile, sep=",");

	return( read.csv(outputFile) );
}

#--------------------------------------------------------------
# Public:export_csv3
#
#--------------------------------------------------------------
export_csv3<-function( res, outputFile, bChi2=FALSE  )
{
	if (res$type != 3 )
	{
		sErrMsg<- "Error: Not a result for 3 SNPs interaction.";
		stop( sErrMsg );
	}

	if (!bChi2)
	{
		items <- data.frame( res$SnpNames[ unlist(res$model[,1]) ],
							 res$SnpNames[ unlist(res$model[,2]) ],
							 res$SnpNames[ unlist(res$model[,3]) ],
							 res$model[, c( 30:55 ) ] );
		colnames(items)<-  c( "SNP1", "SNP2", "SNP3", paste("pv-", EFACTOR3, sep=""));
	}
	else
	{
		items <- data.frame( res$SnpNames[ unlist(res$model[,1]) ],
							 res$SnpNames[ unlist(res$model[,2]) ],
							 res$SnpNames[ unlist(res$model[,3]) ],
							 res$model[, c( 4:29 ) ] );
		colnames(items)<-  c( "SNP1", "SNP2", "SNP3", paste("chi2-", EFACTOR3, sep=""));
	}

	write.table(as.matrix(items), file=outputFile, sep=",");

	return( read.csv(outputFile) );
}

#--------------------------------------------------------------
# private:pSummaryTable1
#
#------
# 1 SNPs interactions( SNP)
#
#         AA:Aa:aa AA:aa AA+aa:2Aa
# chi2  2.76   2.89  3.01  5.10  12.11    2.31     2.33     2.11
# p-v   0.10   0.09  0.06  0.03  <=0.001  ...      ...      ...
#--------------------------------------------------------------
pSummaryTable1<-function(ret, outputFile)
{
	str_list <- list();
	index <- 1;
	model <- ret$model;

	for (i in 1:length(model))
	{
		ret1 <- model[[i]];
		str0 <- paste( "#---------------------------------------------------------------------\n",
			  "Snp1=", ret1$snp1_id,"\n", sep="");

		line1 <- pFormatSummary( "%10s", list("AA:Aa:aa", "Additive", "Dominance" ) );
		line2 <- pFormatSummary( "%10s", list(ret1$chi2_gen3, ret1$chi2_AA_to_aa, ret1$chi2_AAaa_to_Aa ) );
		line3 <- pFormatSummary( "%10s", list(ret1$p_gen3, ret1$p_AA_to_aa, ret1$p_AAaa_to_Aa ) );
		str1  <- paste( "      ",line1, "chi2  ", line2, "p-v   ", line3, "Bonfer", line4, "Holm  ", line5, "FDR   ", line6, sep="");

		str_list[index] <- paste( str0, str1, sep="");
		index<-index+1;

	}

	str<-"";
	for ( i in 1:length(str_list) )
		str <- paste( str, str_list[[i]],sep="");

	if (outputFile!="")
	{
		zz <- file( outputFile, "w");
		cat(str, file=zz, sep="");
		close(zz);
	}

	return(str);
}


#--------------------------------------------------------------
# private:pSummaryTable2
#
#------
# 2 SNPs interactions( SNP1 X SNP2 )
#
# SNP: snpx1, snpx2
#       a1     a2    d1    d2    i_a1a2   i_a1d2   i_d1a2   i_d1d2
# chi2  2.76   2.89  3.01  5.10  12.11    2.31     2.33     2.11
# p-v   0.10   0.09  0.06  0.03  <=0.001  ...      ...      ...
#--------------------------------------------------------------
pSummaryTable2<-function(ret, outputFile)
{
	str_list <- list();
	index <- 1;
	model <- ret$model;

	for (i in 1:length(model))
	{
		ret2 <- model[[i]];
		str0 <- paste( "#---------------------------------------------------------------------\n",
			  "Snp1=", ret2$snp1_id, ",","Snp2=", ret2$snp2_id, "\n", sep="");

		line1 <- pFormatSummary( "%10s", list("a1","a2","d1","d2","a1a2","a1d2","d1a2","d1d2") );
		line2 <- pFormatSummary( "%10s", list(ret2$chi2_a1, ret2$chi2_a2, ret2$chi2_d1, ret2$chi2_d2, ret2$chi2_aa, ret2$chi2_ad, ret2$chi2_da, ret2$chi2_dd ) );
		line3 <- pFormatSummary( "%10s", list(ret2$p_a1, ret2$p_a2, ret2$p_d1, ret2$p_d2, ret2$p_aa, ret2$p_ad, ret2$p_da, ret2$p_dd ) );
		str1  <- paste( "      ",line1, "chi2  ", line2, "p-v   ", line3, sep="");

		str_list[index] <- paste( str0, str1, sep="");
		index <- index+1;

	}

	str<-"";
	for ( i in 1:length(str_list) )
		str <- paste( str, str_list[[i]],sep="");

	if (outputFile!="")
	{
		zz <- file( outputFile, "w");
		cat(str, file=zz, sep="");
		close(zz);
	}


	return(str);
}

#--------------------------------------------------------------
# private:pSummaryTable3
#
#------
# SNP: snpx1, snpx2, snpx3
#       a1          a1a2       a1d2       a1a3       a1d3       a1a2a3     a1a2d3     a1d2a3     a1d2d3
# chi2  +00.12345  -00.12345  -00.12345  -00.12345  -00.12345  -00.12345  -00.12345  -00.12345  -00.12345
# p-v   *<0.0001   *<0.0001    ...        ...        ...        ...        ...        ...        ...
#       d1          d1a2       d1d2       d1a3       d1d3       d1a2a3     d1a2d3     d1d2a3     d1d2d3
# chi2  +xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx
# p-v   *x.xxxxx   *x.xxxxx    ...        ...        ...        ...        ...        ...        ...
#
#       a2          (a1a2)     (d1a2)     a2a3       a2d3       (a1a2a3)   (a1a2d3)   (d1a2a3)   (d1a2d3)
# chi2  +xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx
# p-v   *x.xxxxx   *x.xxxxx    ...        ...        ...        ...        ...        ...        ...
#       d2          (a1d2)     (d1d2)     d2a3       d2d3       (a1d2a3)   (a1d2d3)   (d1d2a3)   (d1d2d3)
# chi2  +xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx
# p-v   *x.xxxxx   *x.xxxxx    ...        ...        ...        ...        ...        ...        ...
#
#       a3          (a1a3)     (d1a3)     (a2a3)     (d2a3)     (a1a2a3)   (a1d2d3)   (d1a2a3)   (d1d2a3)
# chi2  +xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx
# p-v   *x.xxxxx   *x.xxxxx    ...        ...        ...        ...        ...        ...        ...
#       d3          (a1d3)     (d1d3)     (a2d3)     (d2d3)     (a1a2d3)   (a1d2d3)   (d1a2d3)   (d1d2d3)
# chi2  +xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx  -xx.xxxxx
# p-v   *x.xxxxx   *x.xxxxx    ...        ...        ...        ...        ...        ...        ...
#--------------------------------------------------------------
pSummaryTable3<-function(ret, outputFile="")
{
	str_list<-list();
	index <- 1;
	model <- ret$model;

	for (i in 1:length(model))
	{
		ret3<-model[[i]];
		str0 <- paste( "#-----------------------------------------------------------------------------\n",
					  "Snp1=", ret3$snp1_id, ",","Snp2=", ret3$snp2_id,  ",", "Snp3=", ret3$snp3_id, "\n", sep="");

		line1 <- pFormatSummary( "%10s", list("a1","a1a2","a1d2","a1a3","a1d3","a1a2a3","a1a2d3","a1d2a3","a1d2d3") );
		line2 <- pFormatSummary( "%10s", list(ret3$chi2_a1, ret3$chi2_a1a2, ret3$chi2_a1d2, ret3$chi2_a1a3, ret3$chi2_a1d3, ret3$chi2_aaa, ret3$chi2_aad, ret3$chi2_ada, ret3$chi2_add) );
		line3 <- pFormatSummary( "%10s", list(ret3$p_a1, ret3$p_a1a2, ret3$p_a1d2, ret3$p_a1a3, ret3$p_a1d3, ret3$p_aaa, ret3$p_aad, ret3$p_ada, ret3$p_add) );
		str1  <- paste( "      ",line1, "chi2  ", line2, "p-v   ", line3, sep="");

		line1 <- pFormatSummary( "%10s", list("d1","d1a2","d1d2","d1a3","d1d3","d1a2a3","d1a2d3","d1d2a3","d1d2d3") );
		line2 <- pFormatSummary( "%10s", list(ret3$chi2_d1, ret3$chi2_d1a2, ret3$chi2_d1d2, ret3$chi2_d1a3, ret3$chi2_d1d3, ret3$chi2_daa, ret3$chi2_dad, ret3$chi2_dda, ret3$chi2_ddd) );
		line3 <- pFormatSummary( "%10s", list(ret3$p_d1, ret3$p_d1a2, ret3$p_d1d2, ret3$p_d1a3, ret3$p_d1d3, ret3$p_daa, ret3$p_dad, ret3$p_dda, ret3$p_ddd) );
		str2  <- paste( "      ",line1, "chi2  ", line2, "p-v   ", line3, sep="");

		line1 <- pFormatSummary( "%10s", list("a2","(a1a2)","(d1a2)","a2a3","a2d3","(a1a2a3)","(a1a2d3)","(d1a2a3)","(d1a2d3)" ) );
		line2 <- pFormatSummary( "%10s", list(ret3$chi2_a2, ret3$chi2_a1a2,ret3$chi2_d1a2,ret3$chi2_a2a3,ret3$chi2_a2d3,ret3$chi2_aaa,ret3$chi2_aad,ret3$chi2_daa,ret3$chi2_dad ) );
		line3 <- pFormatSummary( "%10s", list(ret3$p_a2, ret3$p_a1a2, ret3$p_d1a2, ret3$p_a2a3, ret3$p_a2d3, ret3$p_aaa, ret3$p_aad, ret3$p_daa, ret3$p_dad ) );
		str3  <- paste( "      ",line1, "chi2  ", line2, "p-v   ", line3, sep="");

		line1 <- pFormatSummary( "%10s", list("d2","(a1d2)","(d1d2)","d2a3","d2d3","(a1d2a3)","(a1d2d3)","(d1d2a3)","(d1d2d3)" ) );
		line2 <- pFormatSummary( "%10s", list(ret3$chi2_d2, ret3$chi2_a1d2, ret3$chi2_d1d2, ret3$chi2_d2a3, ret3$chi2_d2d3, ret3$chi2_ada, ret3$chi2_add, ret3$chi2_dda, ret3$chi2_ddd ) );
		line3 <- pFormatSummary( "%10s", list(ret3$p_d2, ret3$p_a1d2, ret3$p_d1d2, ret3$p_d2a3, ret3$p_d2d3, ret3$p_ada, ret3$p_add, ret3$p_dda, ret3$p_ddd) );
		str4  <- paste( "      ",line1, "chi2  ", line2, "p-v   ", line3, sep="");

		line1 <- pFormatSummary( "%10s", list("a3","(a1a3)","(d1a3)","(a2a3)","(d2a3)","(a1a2a3)","(a1d2d3)","(d1a2a3)","(d1d2a3)" ) );
		line2 <- pFormatSummary( "%10s", list(ret3$chi2_a3, ret3$chi2_a1a3, ret3$chi2_d1a3, ret3$chi2_a2a3, ret3$chi2_d2a3, ret3$chi2_aaa, ret3$chi2_add, ret3$chi2_daa, ret3$chi2_dda ) );
		line3 <- pFormatSummary( "%10s", list(ret3$p_a3, ret3$p_a1a3, ret3$p_d1a3, ret3$p_a2a3, ret3$p_d2a3, ret3$p_aaa, ret3$p_add, ret3$p_daa, ret3$p_dda ) );
		str5  <- paste( "      ",line1, "chi2  ", line2, "p-v   ", line3, sep="");

		line1 <- pFormatSummary( "%10s", list("d3","(a1d3)","(d1d3)","(a2d3)","(d2d3)","(a1a2d3)","(a1d2d3)","(d1a2d3)","(d1d2d3)" ) );
		line2 <- pFormatSummary( "%10s", list(ret3$chi2_d3, ret3$chi2_a1d3, ret3$chi2_d1d3, ret3$chi2_a2d3, ret3$chi2_d2d3, ret3$chi2_aad, ret3$chi2_add,ret3$chi2_dad, ret3$chi2_ddd ) );
		line3 <- pFormatSummary( "%10s", list(ret3$p_d3, ret3$p_a1d3, ret3$p_d1d3, ret3$p_a2d3, ret3$p_d2d3, ret3$p_aad, ret3$p_add,ret3$p_dad, ret3$p_ddd) );
		str6  <- paste( "      ",line1, "chi2  ", line2, "p-v   ", line3, sep="");

		str_list[index] <- paste( str0, str1, str2, str3, str4, str5, str6, sep="");
		index<-index+1;
	}

	str<-"";
	for ( i in 1:length(str_list) )
		str<-paste( str, str_list[[i]],sep="");

	if (outputFile!="")
	{
		zz <- file( outputFile, "w");
		cat(str, file=zz, sep="");
		close(zz);
	}

	return(str);
}

#--------------------------------------------------------------
# private:pFormatSummary
#
#--------------------------------------------------------------
pFormatSummary<-function(fmt, list1)
{
	str<-"";
	for (i in 1:length(list1))
	{
		if (class(list1[[i]])=="numeric")
			str0<-sprintf( "%2.3f", list1[[i]] )
		else
			str0<-list1[[i]];
		str<-paste( str, sprintf(fmt, str0),sep="" );
	}

	return( paste(str, "\n",sep="") );
}


#--------------------------------------------------------------
# private: pChiseq
#
#--------------------------------------------------------------
pChiseq<- function( case_ctrl_matrix )
{
	dm<-dim(case_ctrl_matrix);
	if (dm[1]!=2)
	{
		stop("Error: The row of Matrix is greater than 2.");
		return(-999);
	}

	n <- dm[1];
	m <- dm[2];

	M <- array(0, dim=c(n+1, m+1));
	M[1:n,1:m]<-case_ctrl_matrix;
	M[ M=0 ] <- 0.5;

	N <- sum(M[1,])+sum(M[2,]);
	M[1,m+1] <-  sum(M[1,1:m])/N;
	M[2,m+1] <-  sum(M[2,1:m])/N;
	for (i in 1:m )
		M[3,i] <- (M[1,i]+M[2,i])/N;

	x2 <- 0;
	for (j in 1:m )
		for (i in 1:2 )
		{
			if ( (M[n+1,j]!=0) && (M[i,m+1]!=0) && (M[ i, j]!=0) )
				x2 <- x2 + (M[ i, j]-M[n+1,j]*M[i,m+1]*N)^2/( M[n+1,j]*M[i,m+1]*N )
		}

	M[n+1,m+1] <- x2;
	return (x2);
}

#--------------------------------------------------------------
# private: pInteractionSnps1
#
#           1 SNPs
#--------------------------------------------------------------
pInteractionSnps1<-function( case_f_list, snp1_list, nHypoTests=1 )
{
	na_set <-which(is.na(case_f_list));
	na_set <- union(na_set, which(is.na(snp1_list)) )
	if (length(na_set)>0)
	{
		snp1_list   <- snp1_list[ -na_set ];
		case_f_list <- case_f_list[ -na_set ];
	}

	ret1<-list(
		chi2_gen3     = 0,
		chi2_AA_to_aa = 0,
		chi2_AAaa_to_Aa = 0,
		p_gen3        = 0,
		p_AA_to_aa    = 0,
		p_AAaa_to_Aa  = 0
	);

	#AA:Aa:aa
	r <-pFuncCalcStandard( case_f_list, snp1_list, 1.00, list(c(0), c(1),c(2))  );
	ret1$chi2_gen3 <- r$chi2;
	ret1$p_gen3    <- 1 - pchisq( ret1$chi2_gen3, (3-1)*(2-1) ) ;

	#AA:aa
	ret1$chi2_AA_to_aa <- pFuncCalc1( case_f_list, snp1_list, 1.00, c(0), c(2)  );
	ret1$p_AA_to_aa    <- 1-pchisq( ret1$chi2_AA_to_aa, 1 );

	#AA+aa:Aa
	ret1$chi2_AAaa_to_Aa <- pFuncCalc1( case_f_list, snp1_list, 1.00, c(0,2), c(1,1) );
	ret1$p_AAaa_to_Aa    <- 1-pchisq( ret1$chi2_AAaa_to_Aa, 1 );

	return( ret1 )
}


#--------------------------------------------------------------
# private: 2 SNPs interactions
#               a1, a2, d1, d2, a*a a*d d*a d*d
#--------------------------------------------------------------
pInteractionSnps2<-function( case_f_list, snp1_list, snp2_list, snp1_col, snp2_col )
{
	na_set <-which(is.na(case_f_list));
	na_set <- union(na_set, which(is.na(snp1_list)) )
	na_set <- union(na_set, which(is.na(snp2_list)) )
	if (length(na_set)>0)
	{
		snp1_list   <- snp1_list[ -na_set ];
		snp2_list   <- snp2_list[ -na_set ];
		case_f_list <- case_f_list[ -na_set ];
	}


	nSampSize<- length(case_f_list);
	a1 <- pGetMAF(snp1_list);
	a2 <- pGetMAF(snp2_list);

    ret2<-list(
      snp1_id = snp1_col,
      snp2_id = snp2_col,
      chi2_a1 = 0,
      chi2_a2 = 0,
      chi2_d1 = 0,
      chi2_d2 = 0,
      chi2_aa = 0,
      chi2_ad = 0,
      chi2_da = 0,
      chi2_dd = 0,
      p_a1 = 0,
      p_a2 = 0,
      p_d1 = 0,
      p_d2 = 0,
      p_aa = 0,
      p_ad = 0,
      p_da = 0,
      p_dd = 0
    );

   #a1
   ret2$chi2_a1 <- pFuncCalc2( case_f_list, snp1_list, snp2_list, 1.00, c(00,02), c(20,22)  );
   ret2$p_a1    <- pProbLookup2( ret2$chi2_a1, N2_A1, nSampSize, a1, a2 );

   #a2
   ret2$chi2_a2 <- pFuncCalc2( case_f_list, snp1_list, snp2_list, 1.00, c(00,20), c(02,22)  );
   ret2$p_a2    <- pProbLookup2( ret2$chi2_a2, N2_A2, nSampSize, a1, a2  );

   #d1
   ret2$chi2_d1 <- pFuncCalc2( case_f_list, snp1_list, snp2_list, 0.50, c(00,02,20,22), c(10,10,12,12) );
   ret2$p_d1    <- pProbLookup2( ret2$chi2_d1, N2_D1, nSampSize, a1, a2  );

   #d2
   ret2$chi2_d2 <- pFuncCalc2( case_f_list, snp1_list, snp2_list, 0.50, c(00,02,20,22), c(01,01,21,21) );
   ret2$p_d2    <- pProbLookup2( ret2$chi2_d2, N2_D2, nSampSize, a1, a2  );

   #a*a
   ret2$chi2_aa <- pFuncCalc2( case_f_list, snp1_list, snp2_list, 1.00, c(22,00), c(20,02)  );
   ret2$p_aa    <- pProbLookup2( ret2$chi2_aa, N2_A1A2, nSampSize, a1, a2  );

   #a*d
   ret2$chi2_ad <- pFuncCalc2( case_f_list, snp1_list, snp2_list, 0.50, c(21,21,02,00), c(01,01,20,22)  );
   ret2$p_ad    <- pProbLookup2( ret2$chi2_ad, N2_A1D2, nSampSize, a1, a2  );

   #d*a
   ret2$chi2_da <- pFuncCalc2( case_f_list, snp1_list, snp2_list, 0.50, c(12,12,20,00), c(10,10,02,22)  );
   ret2$p_da    <- pProbLookup2( ret2$chi2_da, N2_D1A2, nSampSize, a1, a2  );

   #d*d
   ret2$chi2_dd <- pFuncCalc2( case_f_list, snp1_list, snp2_list, 0.25, c(11,11,11,11,00,02,20,22), c(10,10,12,12,01,01,21,21) );
   ret2$p_dd    <- pProbLookup2( ret2$chi2_dd, N2_D1D2, nSampSize, a1, a2  );

   return( ret2 );
}

#--------------------------------------------------------------
# private: 3 SNPs interactions
#               a*a*a a*a*d a*d*a a*d*d d*a*a d*a*d d*d*a d*d*d
#--------------------------------------------------------------
pInteractionSnps3<-function( case_f_list, snp1_list, snp2_list, snp3_list, snp1_col, snp2_col, snp3_col )
{
	na_set <-which(is.na(case_f_list));
	na_set <- union(na_set, which(is.na(snp1_list)) )
	na_set <- union(na_set, which(is.na(snp2_list)) )
	na_set <- union(na_set, which(is.na(snp3_list)) )
	if (length(na_set)>0)
	{
		snp1_list   <- snp1_list[ -na_set ];
		snp2_list   <- snp2_list[ -na_set ];
		snp3_list   <- snp3_list[ -na_set ];
		case_f_list <- case_f_list[ -na_set ];
	}

	nSampSize <- length( case_f_list );
	a1 <- pGetMAF( snp1_list );
	a2 <- pGetMAF( snp2_list );
	a3 <- pGetMAF( snp3_list );

    ret3<-list(
      snp1_id = snp1_col,
      snp2_id = snp2_col,
      snp3_id = snp3_col,
      chi2_a1 = 0,
      chi2_a2 = 0,
      chi2_a3 = 0,
      chi2_d1 = 0,
      chi2_d2 = 0,
      chi2_d3 = 0,
      chi2_a1a2 = 0,
      chi2_a2a3 = 0,
      chi2_a1a3 = 0,
      chi2_d1d2 = 0,
      chi2_d2d3 = 0,
      chi2_d1d3 = 0,
      chi2_a1d2 = 0,
      chi2_a1d3 = 0,
      chi2_a2d3 = 0,
      chi2_d1a2 = 0,
      chi2_d2a3 = 0,
      chi2_d1a3 = 0,
      chi2_aaa = 0,
      chi2_aad = 0,
      chi2_ada = 0,
      chi2_add = 0,
      chi2_daa = 0,
      chi2_dad = 0,
      chi2_dda = 0,
      chi2_ddd = 0,
      p_a1 = 0,
      p_a2 = 0,
      p_a3 = 0,
      p_d1 = 0,
      p_d2 = 0,
      p_d3 = 0,
      p_a1a2 = 0,
      p_a2a3 = 0,
      p_a1a3 = 0,
      p_d1d2 = 0,
      p_d2d3 = 0,
      p_d1d3 = 0,
      p_a1d2 = 0,
      p_a1d3 = 0,
      p_a2d3 = 0,
      p_d1a2 = 0,
      p_d2a3 = 0,
      p_d1a3 = 0,
      p_aaa = 0,
      p_aad = 0,
      p_ada = 0,
      p_add = 0,
      p_daa = 0,
      p_dad = 0,
      p_dda = 0,
      p_ddd = 0
    );

	## a1
	ret3$chi2_a1 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1, c(200,202,220,222),  c(000,002,020,022)  )
	ret3$p_a1    <- pProbLookup3( ret3$chi2_a1, N3_A1, nSampSize, a1, a2, a3 );

	## a2
	ret3$chi2_a2 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1, c(020,022,220,222),  c(000,002,200,202)  )
	ret3$p_a2    <- pProbLookup3( ret3$chi2_a2, N3_A2, nSampSize, a1, a2, a3);

	## a3
	ret3$chi2_a3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1, c(002,022,202,222),  c(000,020,200,220)  )
	ret3$p_a3    <- pProbLookup3( ret3$chi2_a3, N3_A3, nSampSize, a1, a2, a3);

	## d1
	ret3$chi2_d1 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(100,100,102,102,120,120,122,122), c(000,002,020,022,200,202,220,222) )
	ret3$p_d1    <- pProbLookup3( ret3$chi2_d1, N3_D1, nSampSize, a1, a2, a3 );

	## d2
	ret3$chi2_d2 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(010,010,012,012,210,210,212,212), c(000,002,020,022,200,202,220,222) )
	ret3$p_d2    <- pProbLookup3( ret3$chi2_d2, N3_D2, nSampSize, a1, a2, a3 );

	## d3
	ret3$chi2_d3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(001,001,021,021,201,201,221,221), c(000,002,020,022,200,202,220,222) )
	ret3$p_d3    <- pProbLookup3( ret3$chi2_d3, N3_D3, nSampSize, a1, a2, a3 );

	## a1a2
	ret3$chi2_a1a2 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1, c(000,002,220,222),  c(020,022,200,202) )
	ret3$p_a1a2    <- pProbLookup3( ret3$chi2_a1a2, N3_A1A2, nSampSize, a1, a2, a3 );

	## a2a3
	ret3$chi2_a2a3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1, c(000,022,200,222),  c(002,020,202,220) )
	ret3$p_a2a3    <- pProbLookup3( ret3$chi2_a2a3, N3_A2A3, nSampSize, a1, a2, a3 );

	## a1a3
	ret3$chi2_a1a3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1, c(000,020,202,222),  c(002,022,200,220) )
	ret3$p_a1a3    <- pProbLookup3( ret3$chi2_a1a3, N3_A1A3, nSampSize, a1, a2, a3  );

	## d1d2
	ret3$chi2_d1d2 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/4, c(000,002,020,022,110,110,110,110,112,112,112,112,200,202,220,222), c(010,010,012,012,100,100,102,102,120,120,122,122,210,210,212,212) )
	ret3$p_d1d2    <- pProbLookup3( ret3$chi2_d1d2, N3_D1D2, nSampSize, a1, a2, a3 );

	## d2d3
	ret3$chi2_d2d3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/4, c(000,002,011,011,011,011,020,022,200,202,211,211,211,211,220,222), c(001,001,010,010,012,012,021,021,201,201,210,210,212,212,221,221) )
	ret3$p_d2d3    <- pProbLookup3( ret3$chi2_d2d3, N3_D2D3, nSampSize, a1, a2, a3 );

	## d1d3
	ret3$chi2_d1d3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/4, c(000,002,020,022,101,101,101,101,121,121,121,121,200,202,220,222), c(001,001,021,021,100,100,102,102,120,120,122,122,201,201,221,221) )
	ret3$p_d1d3    <- pProbLookup3( ret3$chi2_d1d3, N3_D1D3, nSampSize, a1, a2, a3  );

	## a1d2
	ret3$chi2_a1d2 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(000,002,020,022,210,210,212,212), c(010,010,012,012,200,202,220,222)  )
	ret3$p_a1d2    <- pProbLookup3( ret3$chi2_a1d2, N3_A1D2, nSampSize, a1, a2, a3 );

	## a1d3
	ret3$chi2_a1d3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(000,002,020,022,201,201,221,221), c(001,001,021,021,200,202,220,222)  )
	ret3$p_a1d3    <- pProbLookup3( ret3$chi2_a1d3, N3_A1D3, nSampSize, a1, a2, a3 );

	## a2d3
	ret3$chi2_a2d3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(000,002,021,021,200,202,221,221), c(001,001,020,022,201,201,220,222)  )
	ret3$p_a2d3    <- pProbLookup3( ret3$chi2_a2d3, N3_A2D3, nSampSize, a1, a2, a3 );

	## d1a2
	ret3$chi2_d1a2 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(000,002,120,120,122,122,200,202), c(020,022,100,100,102,102,220,222)  )
	ret3$p_d1a2    <- pProbLookup3( ret3$chi2_d1a2, N3_D1A2, nSampSize, a1, a2, a3 );

	## d1a3
	ret3$chi2_d1a3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(000,020,102,102,122,122,200,220), c(002,022,100,100,120,120,202,222)  )
	ret3$p_d1a3    <- pProbLookup3( ret3$chi2_d1a3, N3_D1A3, nSampSize, a1, a2, a3 );

	## d2a3
	ret3$chi2_d2a3 <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(000,012,012,020,200,212,212,220), c(002,010,010,022,202,210,210,222)  )
	ret3$p_d2a3    <- pProbLookup3( ret3$chi2_d2a3, N3_D2A3, nSampSize, a1, a2, a3 );

	## a1*a2*a3
	ret3$chi2_aaa <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1, c(222,200,020,002),  c(220,202,022,000)  )
	ret3$p_aaa    <- pProbLookup3( ret3$chi2_aaa, N3_A1A2A3, nSampSize, a1, a2, a3 );

	## a1*a2*d3
	ret3$chi2_aad <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(221,221,202,200,022,020,001,001), c(222,220,201,201,021,021,002,000)  )
	ret3$p_aad    <- pProbLookup3( ret3$chi2_aad, N3_A1A2D3, nSampSize, a1, a2, a3 );

	## a1*d2*a3
	ret3$chi2_ada <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(220,212,212,200,022,010,010,002), c(222,210,210,202,020,012,012,000)  )
	ret3$p_ada    <- pProbLookup3( ret3$chi2_ada, N3_A1D2A3, nSampSize, a1, a2, a3 );

	## a1*d2*d3
	ret3$chi2_add <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/4, c(222,220,211,211,211,211,202,200,201,201,012,012,010,010,001,001), c(221,221,212,212,210,210,201,201,022,020,011,011,011,011,002,000)  )
	ret3$p_add    <- pProbLookup3( ret3$chi2_add, N3_A1D2D3, nSampSize, a1, a2, a3 );

	## d1*a2*a3
	ret3$chi2_daa <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/2, c(220,202,122,122,100,100,020,002 ), c(222,200,120,120,102,102,022,000)  )
	ret3$p_daa    <- pProbLookup3( ret3$chi2_daa, N3_D1A2A3, nSampSize, a1, a2, a3  );

	## d1*a2*d3
	ret3$chi2_dad <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/4, c(222,220,201,201,121,121,121,121,102,102,100,100,022,020,001,001), c(221,221,202,200,122,122,120,120,101,101,101,101,021,021,002,000) )
	ret3$p_dad    <- pProbLookup3( ret3$chi2_dad, N3_D1A2D3, nSampSize, a1, a2, a3 );

	## d1*d2*a3
	ret3$chi2_dda <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list, 1/4, c(222,210,210,202,120,120,112,112,112,112,100,100,022,010,010,002), c(220,212,212,200,122,122,110,110,110,110,102,102,020,012,012,000)  )
	ret3$p_dda    <- pProbLookup3( ret3$chi2_dda, N3_D1D2A3, nSampSize, a1, a2, a3 );

	## d1*d2*d3
	ret3$chi2_ddd <- pFuncCalc3(case_f_list, snp1_list, snp2_list, snp3_list,  1/8,
	c(221,221,212,212,210,210,201,201,122,122,120,120,111,111,111,111,111,111,111,111,102,102,100,100,021,021,012,012,010,010,001,001),
	c(222,220,211,211,211,211,202,200,121,121,121,121,112,112,112,112,110,110,110,110,101,101,101,101,022,020,011,011,011,011,002,000) )
	ret3$p_ddd    <- pProbLookup3( ret3$chi2_ddd, N3_D1D2D3, nSampSize, a1, a2, a3 );

	return(ret3);
}

#--------------------------------------------------------------
# private: pGetGenStrs1
#
#--------------------------------------------------------------
pGetGenStrs1<- function(GenAa,mode)
{
	genA <- substring(GenAa,1,1)
	gena <- substring(GenAa,2,2)

	if (mode==0)
		x <- paste(genA,genA,sep="")
	else if (mode==2)
		x <- paste(gena,gena,sep="")
	else
		x <- paste(genA,gena,sep="");

	return(x);
}

#--------------------------------------------------------------
# Private:pFuncCalc1
#
#--------------------------------------------------------------
pFuncCalc1<- function( case_f, snp1, r, cols1, cols2 )
{
   	GenAa <- pFindHeterozygousGenes(snp1);
   	snp1  <- pReviseHeterozygousGenes(snp1, GenAa);
   	snp   <- paste(case_f, snp1, sep="");

   	M <- array(0, dim=c(2,2))
   	for (j in 1:length(cols1))
   	{
      	genStr <- pGetGenStrs1(GenAa,cols1[j])
	  	M[1,1] <- M[1,1] + r*length( snp[snp==paste("1", genStr, sep="")] )
      	M[2,1] <- M[2,1] + r*length( snp[snp==paste("0", genStr, sep="")] )
	}

	for (j in 1:length(cols2))
	{
		genStr <- pGetGenStrs1(GenAa,cols2[j])
		M[1,2] <- M[1,2] + r*length( snp[snp==paste("1", genStr, sep="")] )
        M[2,2] <- M[2,2] + r*length( snp[snp==paste("0", genStr, sep="")] )
	}

    chi2 <- pChiseq(M);
    return( chi2 );
}


#--------------------------------------------------------------
# Private:pGetGenStrs2
#
#--------------------------------------------------------------
pGetGenStrs2<- function(GenAa,GenBb,mode)
{
	genA <- substring(GenAa,1,1);
	gena <- substring(GenAa,2,2);
	genB <- substring(GenBb,1,1);
	genb <- substring(GenBb,2,2);

	x <- "";
	if (mode==22)
		x <- paste( genA,genA,genB,genB,sep="")
	else if (mode==21)
		x <- paste( genA,genA,genB,genb,sep="")
	else if (mode==20)
		x <- paste( genA,genA,genb,genb,sep="")
	else if (mode==12)
		x <- paste( genA,gena,genB,genB,sep="")
	else if (mode==11)
		x <- paste( genA,gena,genB,genb,sep="")
	else if (mode==10)
		x <- paste( genA,gena,genb,genb,sep="")
	else if (mode==2)
		x <- paste( gena,gena,genB,genB,sep="")
	else if (mode==1)
		x <- paste( gena,gena,genB,genb,sep="")
	else if (mode==0)
		x <- paste( gena,gena,genb,genb,sep="");

	return(x);
}

#--------------------------------------------------------------
# Private:pFuncCalc2
#
#--------------------------------------------------------------
pFuncCalc2<- function( case_f, snp1, snp2, r, cols1, cols2 )
{
	GenAa <- pFindHeterozygousGenes(snp1);
	GenBb <- pFindHeterozygousGenes(snp2);
	snp1  <- pReviseHeterozygousGenes(snp1, GenAa);
	snp2  <- pReviseHeterozygousGenes(snp2, GenBb);
	snp   <- paste(case_f, snp1, snp2, sep="");

	M <- array(0, dim=c(2,2));

	for (j in 1:length(cols1))
	{
		genStr <- pGetGenStrs2(GenAa,GenBb,cols1[j]);
		M[1,1] <- M[1,1] + r*length( snp[snp==paste("1", genStr, sep="")] );
		M[2,1] <- M[2,1] + r*length( snp[snp==paste("0", genStr, sep="")] );
	}

	for (j in 1:length(cols2))
	{
		genStr <- pGetGenStrs2(GenAa,GenBb,cols2[j]);
		M[1,2] <- M[1,2] + r*length( snp[snp==paste("1", genStr, sep="")] );
		M[2,2] <- M[2,2] + r*length( snp[snp==paste("0", genStr, sep="")] );
	}

	chi2<-pChiseq(M);
	return( chi2 );
}

#--------------------------------------------------------------
# Private:pGetGenStrs3
#
#--------------------------------------------------------------
pGetGenStrs3<- function(GenAa,GenBb,GenCc, mode)
{
	genA <- substring( GenAa,1,1 );
	gena <- substring( GenAa,2,2 );
	genB <- substring( GenBb,1,1 );
	genb <- substring( GenBb,2,2 );
	genC <- substring( GenCc,1,1 );
	genc <- substring( GenCc,2,2 );

	GenType <- array("", dim=c(27));
	ModType <- array(0, dim=c(27));

	AAList <- c( paste(gena,gena, sep=""), paste(genA,gena, sep=""), paste(genA,genA, sep="") );
	BBList <- c( paste(genb,genb, sep=""), paste(genB,genb, sep=""), paste(genB,genB, sep="") );
	CCList <- c( paste(genc,genc, sep=""), paste(genC,genc, sep=""), paste(genC,genC, sep="") );
	index  <- 1;

	for (i in 1:3)
	for (j in 1:3)
	for (k in 1:3)
	{
		GenType[index] <- paste(AAList[i], BBList[j], CCList[k], sep="");
		ModType[index] <- (i-1)*100+(j-1)*10+(k-1);
		index <- index+1;
	}

	r <- GenType[ModType==mode];
	return(r[1]);
}

#--------------------------------------------------------------
# Private:pFuncCalc3
#
#--------------------------------------------------------------
pFuncCalc3<- function( case_f, snp1, snp2, snp3, r, cols1, cols2 )
{
	GenAa <- pFindHeterozygousGenes(snp1);
	GenBb <- pFindHeterozygousGenes(snp2);
	GenCc <- pFindHeterozygousGenes(snp3);

    snp1 <- pReviseHeterozygousGenes(snp1, GenAa);
    snp2 <- pReviseHeterozygousGenes(snp2, GenBb);
    snp3 <- pReviseHeterozygousGenes(snp3, GenCc);
    snp  <- paste(case_f, snp1, snp2, snp3,sep="");

    M <- array(0, dim=c(2,2))
    for (j in 1:length(cols1))
    {
		genStr <- pGetGenStrs3(GenAa,GenBb,GenCc,cols1[j]);
		M[1,1] <- M[1,1] + r*length( snp[snp==paste("1", genStr, sep="")] );
        M[2,1] <- M[2,1] + r*length( snp[snp==paste("0", genStr, sep="")] );
	}

    for (j in 1:length(cols2))
    {
		genStr <- pGetGenStrs3(GenAa,GenBb,GenCc,cols2[j]);
		M[1,2] <- M[1,2] + r*length( snp[snp==paste("1", genStr, sep="")] );
        M[2,2] <- M[2,2] + r*length( snp[snp==paste("0", genStr, sep="")] );
	}

    chi2<-pChiseq(M);
    return( chi2 );
}

#--------------------------------------------------------------
# Private:pFuncCalcStandard
#
#--------------------------------------------------------------
pFuncCalcStandard<-function( case_f, snp1, r, row_cells  )
{
	genAa <- pFindHeterozygousGenes(snp1);
	genA  <- substring(genAa,1,1)
	gena  <- substring(genAa,2,2)
	snp1  <- pReviseHeterozygousGenes(snp1, genAa);
	snp   <- paste(case_f, snp1, sep="");

	M <- array(0, dim=c(2,length(row_cells)));

	for (i in 1:length(row_cells))
	{
		cols <- row_cells[[i]];
		for (j in 1:length(cols))
		{
			genStr <- pGetGenStrs1(genAa, cols[j]) ;
			M[1,i] <- M[1,i] + r*length( snp[snp==paste("1", genStr, sep="")] ) ;
			M[2,i] <- M[2,i] + r*length( snp[snp==paste("0", genStr, sep="")] ) ;
		}
	}

	names(M)<-c( paste(genA,genA,sep=""), paste(genA,gena,sep=""), paste(gena,gena,sep="") );

	chi2 <- pChiseq(M);
	return( list(chi2=chi2, M=M, GenA=genA, Gena=gena ) );
}

#--------------------------------------------------------------
# Private:pFindHeterozygousGenes
#
#--------------------------------------------------------------
pFindHeterozygousGenes<-function(snp_list)
{
    GenA  <- "";
    GenAA <- "";
    Gena  <- "";

    for (i in 1:length(snp_list))
    {
    	str <- snp_list[i];

    	if (GenA=="" && str!="")
    	{
    	   GenA  <- substr(str,1,1)
    	   GenAA <- paste(GenA, GenA, sep="");
    	}

    	if (str!=GenAA && Gena=="")
    	{
    	    if (substr(str,1,1)!=GenA)
    	        Gena <- substr(str,1,1)
    	    else
    	        Gena <- substr(str,2,2);
    	}

    }

    return(paste(GenA, Gena, sep=""));
}

#--------------------------------------------------------------
# Private:pReviseHeterozygousGenes
#
#--------------------------------------------------------------
pReviseHeterozygousGenes<-function(snp_list, GenAa)
{
	GenA  <- substr(GenAa,1,1);
	Gena  <- substr(GenAa,2,2);
	GenaA <- paste(Gena, GenA, sep="");

	snp_list <- as.character( snp_list);
	snp_list[snp_list==GenaA] <- GenAa;

	return(snp_list);
}

SAMPLE_SIZE <-c(0, 200, 400, 800, 1000, 2000, 4000, 10000, 1e+20 );

#--------------------------------------------------------------
# Private:pProbLookup3
#
#--------------------------------------------------------------
pProbLookup3<-function( chi2, eFactor, nSampSize, a1, a2, a3)
{
	if (!exists("snp3_pv_table"))
		data(cc_snp3_l100000_s);

	a1 <- round(ifelse( a1>0.5, 1-a1, a1 )*10)/10;
	a2 <- round(ifelse( a2>0.5, 1-a2, a2 )*10)/10;
	a3 <- round(ifelse( a3>0.5, 1-a3, a3 )*10)/10;

	as3 <- sort(c(a1,a2,a3));
	a1 <-  as3[1];
	a2 <-  as3[2];
	a3 <-  as3[3];

	flank_index <- -1;
	for (i in 1:(length(SAMPLE_SIZE)-1))
		if  ( ( nSampSize - SAMPLE_SIZE[i] >= 0) && ( SAMPLE_SIZE[i+1] - nSampSize>= 0) )
			flank_index <- i;

	par1<- which( snp3_pv_table[,1]== SAMPLE_SIZE[flank_index] &
				  snp3_pv_table[,2]== eFactor &
				  snp3_pv_table[,3]== round(a1*10) &
				  snp3_pv_table[,4]== round(a2*10) &
				  snp3_pv_table[,5]== round(a3*10) );
	par2<- which( snp3_pv_table[,1]== SAMPLE_SIZE[flank_index+1] &
				  snp3_pv_table[,2]== eFactor &
				  snp3_pv_table[,3]== round(a1*10) &
				  snp3_pv_table[,4]== round(a2*10) &
				  snp3_pv_table[,5]== round(a3*10) );

	if (par1==0)
	{
		ret <- 1 - pchisq( chi2/snp3_pv_table[ par2 , 6 ], snp3_pv_table[ par2 , 7 ]);
	}
	else if (par2==0)
	{
		ret <- 1 - pchisq( chi2/snp3_pv_table[ par1 , 6 ], snp3_pv_table[ par1 , 7 ]);
	}
	else
	{
		ret2 <- 1 - pchisq( chi2/snp3_pv_table[ par2 , 6 ], snp3_pv_table[ par2 , 7 ]);
		ret1 <- 1 - pchisq( chi2/snp3_pv_table[ par1 , 6 ], snp3_pv_table[ par1 , 7 ]);

		ret  <- (ret2-ret1)/( SAMPLE_SIZE[flank_index+1]- SAMPLE_SIZE[flank_index])*( nSampSize - SAMPLE_SIZE[flank_index])+ret1;
	}

	#cat(chi2, eFactor, nSampSize, a1, a2, a3, ret, "\n");
	#flush.console();

	return (ret);
}

#--------------------------------------------------------------
# Private:pProbLookup2
#
#--------------------------------------------------------------
pProbLookup2<-function(chi2, eFactor, nSampSize, a1, a2 )
{
	if (!exists("snp2_pv_table"))
		data(cc_snp2_l100000_s);

	a1 <- round(ifelse( a1>0.5, 1-a1, a1 )*10)/10;
	a2 <- round(ifelse( a2>0.5, 1-a2, a2 )*10)/10;
	if (a1>a2)
	{
		a  <- a1;
		a1 <- a2;
		a2 <- a;
	}

	flank_index <- -1;
	for (i in 1:(length(SAMPLE_SIZE)-1))
		if  ( ( nSampSize - SAMPLE_SIZE[i] >= 0) && ( SAMPLE_SIZE[i+1] - nSampSize>= 0) )
			flank_index <- i;

	par1<- which( snp2_pv_table[,1]== SAMPLE_SIZE[flank_index] &
				  snp2_pv_table[,2]== eFactor &
				  snp2_pv_table[,3]== round(a1*10) &
				  snp2_pv_table[,4]== round(a2*10));
	par2<- which( snp2_pv_table[,1]== SAMPLE_SIZE[flank_index+1] &
				  snp2_pv_table[,2]== eFactor &
				  snp2_pv_table[,3]== round(a1*10) &
				  snp2_pv_table[,4]== round(a2*10));

	if (length(par1)==0)
		browser();

	if (par1==0)
	{
		ret <- 1 - pchisq( chi2/snp2_pv_table[ par2 , 6 ], snp2_pv_table[ par2 , 7 ]);
	}
	else if (par2==0)
	{
		ret <- 1 - pchisq( chi2/snp2_pv_table[ par1 , 6 ], snp2_pv_table[ par1 , 7 ]);
	}
	else
	{
		ret2 <- 1 - pchisq( chi2/snp2_pv_table[ par2 , 6 ], snp2_pv_table[ par2 , 7 ]);
		ret1 <- 1 - pchisq( chi2/snp2_pv_table[ par1 , 6 ], snp2_pv_table[ par1 , 7 ]);
		ret  <- (ret2-ret1)/( SAMPLE_SIZE[flank_index+1]- SAMPLE_SIZE[flank_index])*( nSampSize - SAMPLE_SIZE[flank_index])+ret1;
	}

	#cat(chi2, eFactor, nSampSize, a1, a2, ret, "\n");
	#flush.console();

	return (ret);
}

#--------------------------------------------------------------
# Private:pFindTable
#
# SampSize=800
#         a*a       a*d       d*a      d*d
# 0.05    3.855868  3.091899  3.021685 2.150806
# 0.04    4.237501  3.394491  3.320417 2.359034
# 0.03    4.731543  3.788948  3.706502 2.634268
# 0.025   5.046130  4.043677  3.948806 2.810318
# 0.02    5.436830  4.363299  4.249326 3.023520
# 0.01    6.660623  5.331260  5.211649 3.705977
# 0.009   6.853475  5.481846  5.354972 3.807189
# 0.008   7.060065  5.654544  5.519499 3.928891
# 0.007   7.299579  5.841686  5.701708 4.059748
# 0.006   7.576921  6.064399  5.922590 4.214378
# 0.005   7.906502  6.334334  6.187868 4.400070
# 0.004   8.306440  6.637344  6.498593 4.631195
# 0.003   8.817538  7.055636  6.915405 4.909662
# 0.002   9.564631  7.628304  7.516572 5.329030
# 0.001  10.893624  8.607966  8.511935 6.087811
# 0.0005 12.206780  9.705727  9.519495 6.809887
# 0.0001 14.983136 11.901417 11.831287 8.390764
#
#
# SampSize=2000
#          a*a       a*d       d*a      d*d
# 0.05    3.838267  3.134887  3.132371 2.073486
# 0.04    4.214602  3.442856  3.445706 2.277714
# 0.03    4.706940  3.853665  3.838678 2.542172
# 0.025   5.016273  4.108493  4.096646 2.712266
# 0.02    5.402035  4.429127  4.415764 2.918676
# 0.01    6.622148  5.413215  5.404164 3.573128
# 0.009   6.809269  5.565340  5.553525 3.669296
# 0.008   7.015477  5.738049  5.723830 3.783566
# 0.007   7.251487  5.930714  5.907291 3.916748
# 0.006   7.530890  6.142392  6.151514 4.060562
# 0.005   7.869016  6.418774  6.410197 4.231042
# 0.004   8.281864  6.746729  6.732095 4.445650
# 0.003   8.795845  7.187406  7.124648 4.743155
# 0.002   9.513416  7.793576  7.727318 5.142597
# 0.001  10.804575  8.836399  8.786469 5.883260
# 0.0005 12.058079  9.990800  9.934400 6.585202
# 0.0001 15.126089 12.511666 12.285011 8.273149
#
# SampSize=800
#         a*a*a     a*a*d     a*d*a     a*d*d     d*a*a     d*a*d     d*d*a    d*d*d
# 0.05    3.877435  3.199245  3.135328  2.595132  3.198639  2.615121  2.587106 2.129040
# 0.04    4.249664  3.509711  3.441647  2.851589  3.508598  2.873018  2.840972 2.334768
# 0.03    4.735003  3.908701  3.838116  3.185248  3.918111  3.208925  3.178045 2.603283
# 0.025   5.046595  4.174722  4.094141  3.398342  4.174079  3.423529  3.381893 2.776175
# 0.02    5.437531  4.497771  4.411300  3.658089  4.491907  3.683341  3.645328 2.990465
# 0.01    6.662896  5.498694  5.397255  4.494223  5.486956  4.510076  4.469643 3.662616
# 0.009   6.863034  5.650404  5.556241  4.617288  5.638456  4.637061  4.593364 3.763877
# 0.008   7.072322  5.829606  5.728404  4.755732  5.813866  4.781635  4.731588 3.881165
# 0.007   7.310011  6.028351  5.920265  4.911473  5.999123  4.939334  4.896962 4.014646
# 0.006   7.574396  6.262412  6.155102  5.090695  6.220452  5.127814  5.079474 4.163860
# 0.005   7.882263  6.517103  6.416252  5.317058  6.492283  5.344820  5.291390 4.348506
# 0.004   8.283909  6.835662  6.733075  5.572280  6.833593  5.617819  5.555630 4.556122
# 0.003   8.786710  7.268163  7.174501  5.935689  7.272453  5.988361  5.890440 4.840763
# 0.002   9.528653  7.866815  7.773540  6.423090  7.857396  6.467224  6.352940 5.274100
# 0.001  10.706591  8.920975  8.784350  7.255002  8.909333  7.329702  7.218076 5.955025
# 0.0005 11.924245  9.976899  9.746644  8.146464  9.942440  8.262644  8.059227 6.617415
# 0.0001 15.008326 12.342043 12.064696 10.163803 12.514748 10.467968 10.104315 8.1672
#
# SampSize=2000
#
#         a*a*a     a*a*d     a*d*a     a*d*d     d*a*a     d*a*d     d*d*a    d*d*d
# 0.05    3.858517  3.195455  3.190595  2.642796  3.184236  2.633363  2.636684 2.151878
# 0.04    4.236325  3.505284  3.501381  2.902956  3.495572  2.892816  2.894264 2.366763
# 0.03    4.728011  3.914996  3.913851  3.239478  3.899793  3.229038  3.230673 2.644231
# 0.025   5.042871  4.178419  4.174962  3.455918  4.156996  3.443873  3.444548 2.820575
# 0.02    5.423972  4.503066  4.493133  3.720023  4.483244  3.714264  3.715074 3.041619
# 0.01    6.649821  5.528974  5.503372  4.549873  5.494257  4.550716  4.545219 3.727712
# 0.009   6.828195  5.683582  5.654279  4.676053  5.639777  4.677565  4.674846 3.834921
# 0.008   7.035489  5.848797  5.833553  4.816649  5.815752  4.822356  4.818260 3.960015
# 0.007   7.273961  6.043149  6.033462  4.984057  6.016431  4.994363  4.982673 4.096417
# 0.006   7.545545  6.267673  6.263585  5.175376  6.257115  5.187808  5.157942 4.257425
# 0.005   7.875212  6.546245  6.536951  5.401855  6.539414  5.420892  5.381921 4.441675
# 0.004   8.274572  6.893243  6.881444  5.677941  6.872947  5.710397  5.654058 4.667981
# 0.003   8.802450  7.329905  7.308081  6.032059  7.313704  6.055910  5.998642 4.958756
# 0.002   9.609165  7.938571  7.897705  6.542112  7.896997  6.574160  6.518016 5.382842
# 0.001  10.895767  8.923579  9.003315  7.403025  8.883473  7.419482  7.414816 6.088396
# 0.0005 12.169890  9.920138 10.068815  8.352486  9.945667  8.291980  8.280507 6.795205
# 0.0001 15.059864 12.196852 12.612459 10.268463 12.452363 10.308919 10.310114 8.495672
#--------------------------------------------------------------
pFindTable<-function( chi2, d)
{
	tb<-array(0, dim=c(17,4));
	#            1d        2d         3d
	tb[1,] <-c( 0.05     ,3.195455   ,2.642796   ,2.151878);
	tb[2,] <-c( 0.04     ,3.505284   ,2.902956   ,2.366763);
	tb[3,] <-c( 0.03     ,3.914996   ,3.239478   ,2.644231);
	tb[4,] <-c( 0.025    ,4.178419   ,3.455918   ,2.820575);
	tb[5,] <-c( 0.02     ,4.503066   ,3.720023   ,3.041619);
	tb[6,] <-c( 0.01     ,5.528974   ,4.549873   ,3.727712);
	tb[7,] <-c( 0.009    ,5.683582   ,4.676053   ,3.834921);
	tb[8,] <-c( 0.008    ,5.848797   ,4.816649   ,3.960015);
	tb[9,] <-c( 0.007    ,6.043149   ,4.984057   ,4.096417);
	tb[10,]<-c( 0.006    ,6.267673   ,5.175376   ,4.257425);
	tb[11,]<-c( 0.005    ,6.546245   ,5.401855   ,4.441675);
	tb[12,]<-c( 0.004    ,6.893243   ,5.677941   ,4.667981);
	tb[13,]<-c( 0.003    ,7.329905   ,6.032059   ,4.958756);
	tb[14,]<-c( 0.002    ,7.938571   ,6.542112   ,5.382842);
	tb[15,]<-c( 0.001    ,8.923579   ,7.403025   ,6.088396);
	tb[16,]<-c( 0.0005   ,9.920138   ,8.352486   ,6.795205);
	tb[17,]<-c( 0.0001  ,12.196852  ,10.268463   ,8.495672);

	for (i in 1:17)
	{
		if ( chi2 < tb[i,d+1] )
		{
			if (i==1)
				#return(">0.05")
				return("0.10")
			else
				#return ( paste("*<", tb[i-1,1], sep="" )) ;
				return ( paste("", tb[i-1,1], sep="" )) ;
		}
	}

	#return("*<=0.0001");
	return("0.0001");
}

#--------------------------------------------------------------
# Private:pGetMAF
#
#--------------------------------------------------------------
pGetMAF<-function(snp_list)
{
    genA <- "";
    cntA <- 0;
    total<- 0;

    for (i in 1:length(snp_list))
    {
    	str<-snp_list[i];
    	if (genA=="" && str!="")
    	   genA <- substr(str,1,1);

   	    if (!is.na(str) && str!="")
   	    {
   	    	if ( substr(str,1,1) ==genA)
   	    		cntA <- cntA + 1;
   	    	if ( substr(str,2,1) ==genA)
   	    		cntA <- cntA + 1;
   	    	total<-total+2;
   	    }
   }

	ret <- ifelse(cntA/total>0.5, 1-cntA/total , cntA/total );
	return(ret);
}

#--------------------------------------------------------------
# public: adjust_by_bonferroni
#
#--------------------------------------------------------------
TIAN.adjust_by_bonferroni<-function( res, coff = NA )
{
	m <- ifelse(is.na(coff), length(res$model[,1]), coff);

	if (res$type==1)
	{
		p <- matrix(data=unlist(res$model[,c(5:7)]), ncol=3, byrow=FALSE);
		res$model <- cbind( matrix(  res$model[,1:4], ncol=4),
					bon_p_gen3 = (p*m)[,1],
					bon_p_AA_to_aa = (p*m)[,2],
					bon_p_AAaa_to_Aa = (p*m)[,3] );
	}
	else if (res$type==2)
	{
		p <- matrix(data=unlist(res$model[,c((3+8):18 )]), ncol=8, byrow=FALSE);
		res$model <- cbind( matrix(  res$model[,1:10], ncol=10),
					bon_p_a1 = (p*m)[,1],
					bon_p_a2 = (p*m)[,2],
					bon_p_d1 = (p*m)[,3],
					bon_p_d2 = (p*m)[,4],
					bon_p_aa = (p*m)[,5],
					bon_p_ad = (p*m)[,6],
					bon_p_da = (p*m)[,7],
					bon_p_dd = (p*m)[,8] );
	}
	else if (res$type==3)
	{

		p <- matrix(data=unlist(res$model[,c((4+26):55 )]), ncol=26, byrow=FALSE);
		res$model <- cbind(matrix( res$model[,1:29], ncol=29),
					bon_p_a1 = (p*m)[,1],
					bon_p_a2 = (p*m)[,2],
					bon_p_a3 = (p*m)[,3],
					bon_p_d1 = (p*m)[,4],
					bon_p_d2 = (p*m)[,5],
					bon_p_d3 = (p*m)[,6],
					bon_p_a1a2 = (p*m)[,7],
					bon_p_a2a3 = (p*m)[,8],
					bon_p_a1a3 = (p*m)[,9],
					bon_p_d1d2 = (p*m)[,10],
					bon_p_d2d3 = (p*m)[,11],
					bon_p_d1d3 = (p*m)[,12],
					bon_p_a1d2 = (p*m)[,13],
					bon_p_a1d3 = (p*m)[,14],
					bon_p_a2d3 = (p*m)[,15],
					bon_p_d1a2 = (p*m)[,16],
					bon_p_d2a3 = (p*m)[,17],
					bon_p_d1a3 = (p*m)[,18],
					bon_p_aaa  = (p*m)[,19],
					bon_p_aad  = (p*m)[,20],
					bon_p_ada  = (p*m)[,21],
					bon_p_add  = (p*m)[,22],
					bon_p_daa  = (p*m)[,23],
					bon_p_dad  = (p*m)[,24],
					bon_p_dda  = (p*m)[,25],
					bon_p_ddd  = (p*m)[,26]);
	}
	else
	{
		sErrmsg <- "The model data must be for 1,2,3 SNPs.(TIAN.Adjust_by_Bonferroni)";
		stop(sErrmsg);
	}

	res$correction<-"Bonferroni";

	return( res );
}

#--------------------------------------------------------------
# public: TIAN.adjust_by_holm
#
#--------------------------------------------------------------
TIAN.adjust_by_holm<-function( res )
{
	m<-length(res$model[,1]);
	cols.name <- colnames(res$model);

	adj_model = c();
	if (res$type==1)
		adj_model <- matrix( res$model[,c((2+3):7)],   ncol=3)
	else if (res$type==2)
		adj_model <- matrix( res$model[,c((3+8):18 )], ncol=8)
	else if (res$type==3)
		adj_model <- matrix( res$model[,c((4+26):55 )],ncol=26)
	else
	{
		sErrmsg <- "The model data must be for 1,2,3 SNPs.(TIAN.adjust_by_holm)";
		stop(sErrmsg);
	}

	for (i in 1:length(adj_model[1,]) )
	{
		newp <- pAdjustByHolm( adj_model[,i] );
		adj_model[,i]<- newp;
	}


	pvalue.start = 0;
	if ( res$type == 1 )
	{
		res$model <- cbind( matrix( res$model[,c(1:4)], ncol=4), matrix(adj_model[,c(1:3)], ncol=3 ) )
		pvalue.start = 5;
	}
	else if ( res$type == 2 )
	{
		res$model <- cbind( matrix( res$model[,c(1:10)],ncol=10), matrix(adj_model[,c(1:8)], ncol=8 ) )
		pvalue.start = 11;
	}
	else if ( res$type == 3 )
	{
		res$model <- cbind( matrix( res$model[,c(1:29)],ncol=29), matrix(adj_model[,c(1:26)], ncol=26) )
		pvalue.start = 30;
	}
	else
		return(NULL):

	for (i in 1:(length(cols.name)-pvalue.start+1) )
	{
		str0 <- substr(cols.name[i],3,100);
		cols.name[i]<-paste("h_", str0, sep="");
	}
	colnames(res$model) <- cols.name;

	res$correction<-"Holm";

	return( res );
}


pAdjustByHolm<-function( p )
{
	m    <- length(p);
	new_p<- array(0, dim=c(m));
	tmp  <- sort( unlist(p), index.return=TRUE);

	for (i in 1:length(tmp$ix))
	{
		k  <- (tmp$ix)[i];
		np <- as.numeric( p[k] ) *(m-i+1);
		new_p[k] <- np;
	}

	return( new_p );
}

#--------------------------------------------------------------
# public: TIAN.adjust_by_FDR
#
#--------------------------------------------------------------
TIAN.adjust_by_fdr<-function( res )
{
	m<-length(res$model[,1]);
	cols.name <- colnames(res$model);

	adj_model = c();
	if ( res$type==1 )
		adj_model <- matrix( res$model[,c(5:7)], ncol=3)
	else if ( res$type==2 )
		adj_model <- matrix( res$model[,c((3+8):18 )], ncol=8)
	else if ( res$type==3 )
		adj_model <- matrix( res$model[,c((4+26):55 )], ncol=26)
	else
	{
		sErrmsg <- "The model data must be for 1,2,3 SNPs.(TIAN.adjust_by_FDR)";
		stop(sErrmsg);
	}

	for (i in 1:length(adj_model[1,]) )
	{
		newp <- pAdjustByFDR( adj_model[,i] );
		adj_model[,i]<- newp;
	}

	if (res$type == 1)
	{
		res$model <- cbind( matrix( res$model[,1:4],  ncol=4), matrix(adj_model[,1:3], ncol=3) )
		pvalue.start = 5;
	}
	else if (res$type == 2)
	{
		res$model <- cbind( matrix( res$model[,1:10], ncol=10), matrix(adj_model[,c(1:8)], ncol=8) )
		pvalue.start = 11;
	}
	else if (res$type == 3)
	{
		res$model <- cbind( matrix( res$model[,1:29], ncol=29), matrix(adj_model[,c(1:26)], ncol=26) )
		pvalue.start = 30;
	}
	else
		return (NULL);

	for (i in 1:(length(cols.name)-pvalue.start+1) )
	{
		str0 <- substr(cols.name[i],3,100);
		cols.name[i]<-paste("f_", str0, sep="");
	}
	colnames(res$model) <- cols.name;

	res$correction<-"FDR";
	return( res );
}

pAdjustByFDR<-function( p )
{
	m     <- length(p);
	new_p0<- array(0, dim=c(m));
	new_p <- array(0, dim=c(m));
	tmp   <- sort(unlist(p), index.return=TRUE);

	n <- length(tmp$ix);
	for (i in 1:n)
	{
		k  <- (tmp$ix)[i];
		np <- as.numeric( p[k]) *n/i;
		new_p[k] <- np;
	}

	return( new_p );
}

#--------------------------------------------------------------
# private: pSummaryInteraction
#
#--------------------------------------------------------------
pSummaryInteraction <- function(snp_list, sig_list)
{
	n        <- length(sig_list[1,]);
	#new_sig  <- sig_list[order(unlist(sig_list[, n-1])), ]
	#sig_list <- matrix(new_sig, ncol=n);
	e_factor <- sig_list[,n-2];
	pv       <- sig_list[,n-1];
	ch2v     <- sig_list[,n];
	snp1     <- sig_list[,1];
	snp2     <- sig_list[,2];
	snp3     <- sig_list[,3];

	ret <- NULL;
	#Model 1
	if (n==4)
		snp_list2 <- data.frame( SNP    = snp_list[ unlist(snp1) ],
								eFactor = EFACTOR1[unlist(e_factor)],
								pValue  = unlist(pv),
								chi2v   = unlist(ch2v) );
	#Model 2
	if (n==5)
		snp_list2 <- data.frame( SNP1   = snp_list[ unlist(snp1) ],
								SNP2    = snp_list[ unlist(snp2) ],
								eFactor = EFACTOR2[unlist(e_factor)],
								pValue  = unlist(pv),
								chi2v   = unlist(ch2v) );
	#Model 3
	if (n==6)
		snp_list2 <- data.frame( SNP1   = snp_list[ unlist(snp1) ],
								SNP2    = snp_list[ unlist(snp2) ],
								SNP3    = snp_list[ unlist(snp3) ],
								eFactor = EFACTOR3[unlist(e_factor)],
								pValue  = unlist(pv),
								chi2v   = unlist(ch2v) );
	str  <- "";
	str0 <- ""
	str1 <- "";
	for(i in 1:length(snp_list2[,1]))
	{
		str <- paste( str, i,":\t", sep="");
		str0 <- snp_list2[i,1];
		if (n>=5)
			str0 <- paste( str0, ",", snp_list2[i,2], sep="");
		if (n>=6)
			str0 <- paste( str0, ",", snp_list2[i,3], sep="");
		if (str0 != str1)
		{
			str <- paste( str, str0, "\t", snp_list2[i,n-2], sep="");
			str1 <- str0;
		}
		else
		{
			strspace <- "";
			for (k in 1:nchar(str0))
				strspace<-paste(strspace, " ", sep="");

			str <- paste( str, strspace, "\t", snp_list2[i,n-2], sep="");
		}

		str <- paste( str, "\tx2=", sprintf("%10.8f", snp_list2[i,n] ), sep="");
		str <- paste( str, "\t", sprintf("%10.8f", snp_list2[i,n-1] ), sep="");
		str <- paste( str, "\n" );
	}

	return( str );
}

#--------------------------------------------------------------
# private: pSummaryHeader
#
#--------------------------------------------------------------
pSummaryHeader <- function( res )
{
	stru <- sprintf("------------------------------------\n");
	str0 <- sprintf("%10s: %s\n", "Data File", 	res$dataFile );
	str1 <- sprintf("%10s: %s\n", "Date", 	Sys.time() );
	str2 <- sprintf("%10s: %-10.0f\n", "Cases", 	res$nCaseCnt );
	str3 <- sprintf("%10s: %-10.0f\n", "Controls", 	res$nCtrlCnt );
	str4 <- sprintf("%10s: %-10.0f\n", "SNPs Count", length( res$SnpNames ) );
	str5 <- sprintf("%10s: %d SNPs\n", "Model", res$type );
	str6 <- sprintf("%10s: %s\n", "Correction", res$correction );
	strd <- sprintf("------------------------------------\n");

	return(paste("\n\nEpistatic Interaction Analysis Report\n",
			stru,str0,str1,str2,str3,str4,str5,str6,strd, "\n", sep=""));
}

#--------------------------------------------------------------
# public: TIAN.draw_correlate2
#
#--------------------------------------------------------------
TIAN.draw_correlate2 <- function( res, pdfFile  )
{
	if (res$type != 2 )
	{
		sErrMsg<- "Error: Not a result for 2 SNPs interaction.";
		stop( sErrMsg );
	}

	pdf(pdfFile);

	for(i in 1:8)
	{
		subdata <- matrix( unlist(res$model[,c(1,2,10+i)]), ncol=3);
		title <- EFACTOR2[i];
		if (res$correction!="")
			title <- paste( EFACTOR2[i], "(", res$correction, ")", sep="");
		pDrawFigure2(  subdata , title, labels=res$SnpNames, log10=TRUE);
	}
	dev.off();

	return(pdfFile);
}

pDrawFigure2<-function( dat, sTitle, labels=NULL, log10=TRUE, significance=0.05)
{
	nLenArm <- max( max(dat[,1]), max(dat[,2]) ) ;

	len <- length(dat[,1]);
	rm <- dat[,3 ];
	rm[ rm>1 ] <- 1;
	rm <- -log10 ( dat[,3 ] );
	rm[ rm<0 ] <- 0;

	if (sys$get_value("legend.peak")!= 0)
	{
		maxv <- -log10( sys$get_value("legend.peak") );
	}
	else
	{
		maxv <- max( rm );
		if (maxv <= -log10(0.05))
			maxv <- -log10(0.04);
	}
	minv <- min( rm );
	minv <- 0;

	ramp <- colorRamp(c( sys$get_value("color.v1"), sys$get_value("color.v0") ));
	cols <- rgb ( ramp(seq(0, 1, length = 1000) ), max=255 );
	mc <- array("#FFFFFF", dim = c( nLenArm, nLenArm));
	cls <- round( (rm-minv)/(maxv-minv)*1000 ) +1;
	cls[which(cls>1000)]<-1000;

	for (n in 1:len)
		mc[ dat[n,1], dat[n,2] ] <- cols[  cls[n] ];

	par(mar=c(2,2,2,2)+0.1);
	plot(c(0, nLenArm*1.45), c(0, nLenArm*1.45), type= "n", xlab="", ylab="", xaxt="n", yaxt="n");
	for (x in 1:nLenArm)
	for (y in 1:nLenArm)
	{

		ox <- (nLenArm-x+1);
		oy <- y;
		if (oy+ox<(nLenArm+1)) next;

		x0 <- sqrt(2)/2*(ox-1) - sqrt(2)/2*(oy-1) + nLenArm*sqrt(2)/2;
		y0 <- 2*(nLenArm-1) - sqrt(2)/2*(oy-1+ox-1) - nLenArm*sqrt(2)/2*0.4;
		xs <- c(x0, x0+sqrt(2)/2, x0, x0-sqrt(2)/2, x0);
		ys <- c(y0, y0+sqrt(2)/2, y0+sqrt(2), y0+sqrt(2)/2, y0);

		if (nLenArm<=3)
		{
			xs <- xs/2 + nLenArm*1.45/2/2;
			ys <- ys/2 + nLenArm*1.45/2/2;
		}

		if ( mc[x,y] != "#FFFFFF" )
			polygon(xs, ys, col=mc[x,y], border="gray", angle=-45)
		else
			polygon(xs, ys, col=mc[x,y], border="white", angle=-45);
	}


	l <- nLenArm * sqrt(2)/2;
	x0 <- nLenArm*sqrt(2)/4;
	for (i in 0:100)
	{
		rect(x0+i*(l/100), nLenArm*sqrt(2)/2*0.1,
		     x0+(i+1)*l/100, nLenArm*sqrt(2)/2*0.2,
		     col=cols[i*10+1], border=cols[i*10+1])
	}

	text( x0+0* (l/100), nLenArm*sqrt(2)/2*0.1-strheight("1"), "1",cex=0.5 );
	x50 <- round( (-log10(0.5) - minv)/(maxv-minv)*100 ) +1;
	if (x50<100)
	{
		segments(x0+x50*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x50*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x50*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".5",cex=0.5 );
	}
	x005 <- round( (-log10(0.05) - minv)/(maxv-minv)*100 ) +1;
	if (x005<=101)
	{
		segments(x0+x005*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x005*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x005*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".05",cex=0.5 );
	}
	x001 <- round( (-log10(0.01) - minv)/(maxv-minv)*100 ) +1;
	if (x001<=101)
	{
		segments(x0+x001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".01",cex=0.5 );
	}
	x0001 <- round( (-log10(0.001) - minv)/(maxv-minv)*100 ) +1;
	if (x0001<=101)
	{
		segments(x0+x0001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x0001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x0001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".001",cex=0.5 );
	}
	x00001 <- round( (-log10(0.0001) - minv)/(maxv-minv)*100 ) +1;
	if (x00001<=101)
	{
		segments(x0+x00001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x00001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x00001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".0001",cex=0.5 );
	}

	x000001 <- round( (-log10(0.000001) - minv)/(maxv-minv)*100 ) +1;
	if (x000001<=101)
	{
		segments(x0+x000001*l/100, nLenArm*sqrt(2)/2*0.1,
		     x0+x000001*l/100, nLenArm*sqrt(2)/2*0.2, col = "white");
		text( x0+x000001*(l/100), nLenArm*sqrt(2)/2*0.1-strheight("1.0"), ".000001",cex=0.5 );
	}

	for (n in 1:len )
	{
		if ( dat[n,3] <= significance )
		{
			ox <- (nLenArm - dat[n,1] + 1);
			oy <- dat[n,2];
			if (oy+ox<(nLenArm+1)) next;

			x0 <- sqrt(2)/2*(ox-1) - sqrt(2)/2*(oy-1) + nLenArm*sqrt(2)/2;
			y0 <- 2*(nLenArm-1) - sqrt(2)/2*(oy-1+ox-1) - nLenArm*sqrt(2)/2*0.4 +sqrt(2)/2;

			if (nLenArm<=3)
			{
				x0 <- x0/2 + nLenArm*1.45/2/2;
				y0 <- y0/2 + nLenArm*1.45/2/2;
			}

			if (nLenArm<10)
			{
				str <- sprintf("%5.3f", dat[n,3]);
				text( x0, y0, str, col=sys$get_value("color.cross"), srt=-45);
			}
			else
			{
				points( x0, y0, pch=3,cex=0.8, col=sys$get_value("color.cross") );
			}
		}
	}

	if (!is.null( labels ))
	{
		for (n in 1:length(labels) )
		{
			x0 <- (n-1)*sqrt(2)+sqrt(2)/2;
			y0 <- sqrt(2)/2*nLenArm * (1.4)-0.5;
			if (nLenArm<=3)
			{
				x0 <- x0/2 + nLenArm*1.45/2/2;
				y0 <- y0/2 + nLenArm*1.45/2/2;
			}

			y0 <- y0 + strwidth(labels[n], srt=90) ;
			text(x0, y0, labels[n], adj=c(1, 0), srt=90);
		}
	}


	title( sTitle );
}

#--------------------------------------------------------------
# public: TIAN.draw_correlate3
#
#--------------------------------------------------------------
TIAN.draw_correlate3 <- function( res, pdfFile  )
{
	if (res$type != 3 )
	{
		sErrMsg<- "Error: Not a result for 3 SNPs interaction.";
		stop( sErrMsg );
	}

	pdf(pdfFile);

	fulldat <- get_fulldat3(res);
	for(i in 1:26)
	{
		subdata <- fulldat[,c(1,2,3, 3+i)];
		title <- EFACTOR3[i];
		if (res$correction!="")
			title <- paste( EFACTOR3[i], "(", res$correction, ")", sep="");
		pDrawFigure3(  subdata , title, res$SnpNames, log10=TRUE);
	}

	dev.off();

	return(pdfFile);
}


get_fulldat3 <- function( res )
{
	fulldat<-array(0, dim=c(0,29));
	nSnp <- length(res$SnpNames);
	for(i in 1:nSnp)
	for(j in 1:nSnp)
	{
		p<- get_correlate3(res, i, j, nSnp);
		fulldat<-rbind( fulldat, p)
	}
	fulldat <- fulldat[ -which(fulldat[,1]==0), ];
	fulldat <- matrix( unlist(fulldat), ncol=length(fulldat[1,]));

	return(fulldat);
}

get_correlate3 <- function( res, snp1, snp2,  snp_cnt)
{
	s3<-array(0, dim=c(6,26));
          #a1,a2,a3,d1,d2,d3,a1a2,a2a3,a1a3,d1d2,d2d3,d1d3,a1d2,a1d3,a2d3,d1a2,d2a3,d1a3,a*a*a,a*a*d,a*d*a,a*d*d,d*a*a,d*a*d,d*d*a,d*d*d
	#1,2,3
	s3[1,]<-c(1, 2, 3, 4, 5, 6, 7,   8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  19,   20,   21,   22,   23,   24,   25,   26);
	#1,3,2
	s3[2,]<-c(1, 3, 2, 4, 6, 5, 9,   8,   7,   12,  11,  10,  14,  13,  17,  18,  15,  16,  19,   21,   20,   22,   23,   25,   24,   26);
	#2,1,3
	s3[3,]<-c(2, 1, 3, 5, 4, 6, 7,   9,   8,   10,  12,  11,  16,  15,  14,  13,  18,  17,  19,   20,   23,   24,   21,   22,   25,   26);
	#2,3,1
	s3[4,]<-c(2, 3, 1, 5, 6, 4, 8,   9,   7,   11,  12,  10,  15,  16,  18,  17,  14,  13,  19,   23,   20,   24,   21,   25,   22,   26);
	#3,2,1
	s3[5,]<-c(3, 2, 1, 6, 5, 4, 8,   7,   9,   11,  10,  12,  17,  18,  16,  15,  13,  14,  19,   23,   21,   25,   20,   24,   22,   26);
	#3,1,2
	s3[6,]<-c(3, 1, 2, 6, 4, 5, 9,   7,   8,   12,  10,  11,  18,  17,  13,  14,  16,  15,  19,   21,   23,   25,   20,   22,   24,   26);

	p <- array(0, dim=c(snp_cnt, 29) );
	for (i in 1:snp_cnt)
	{
		snp3<- i;
		cs  <- sort(c( snp1, snp2, snp3 ));
		set3<- which(res$model[,1]==cs[1] & res$model[,2]==cs[2] & res$model[,3]==cs[3])
		if (length(set3)<=0)
			next;

		pi<-res$model[ set3[1], c((3+26+1):(3+26+26))];
		if (snp1<snp3 && snp3<snp2) pi<-pi[s3[2,]];
		if (snp2<snp1 && snp1<snp3) pi<-pi[s3[3,]];
		if (snp2<snp3 && snp3<snp1) pi<-pi[s3[4,]];
		if (snp3<snp2 && snp2<snp1) pi<-pi[s3[5,]];
		if (snp3<snp1 && snp1<snp2) pi<-pi[s3[6,]];
		p[i,] <- c( snp1, snp2, snp3, c(unlist(pi)));
	}

	return(p)
}

pDrawFigure3 <- function(  dat, sTitle, labels=NULL, log10=TRUE, significance=0.05)
{
	nLenArm <- max( max(dat[,1]), max(dat[,2]), max(dat[,3]) ) ;
	len <- length(dat[,1]);
	rm  <- dat[,4 ];
	rm[ rm>1 ] <- 1;
	rm <- -log10 ( dat[,4 ] );
	rm[ rm<0 ] <- 0;

	if (sys$get_value("legend.peak")!= 0)
	{
		maxv <- -log10( sys$get_value("legend.peak") );
	}
	else
	{
		maxv <- max( rm );
		if (maxv <= -log10(0.05))
			maxv <- -log10(0.04);
	}

	minv <- min( rm );
	minv <- 0;

	ramp <- colorRamp(c( sys$get_value("color.v1"), sys$get_value("color.v0") ));
	cols <- rgb ( ramp(seq(0, 1, length = 1000) ), max=255 );
	mc <- array("#FFFFFF", dim = c( nLenArm, nLenArm));
	for(i in 1:nLenArm)
	for(j in 1:nLenArm)
	{
		setx <- which(dat[,1]==i & dat[,2]==j);
		if (length(setx)>0)
		{
			maxx <- max(rm[ setx ]);
			cls <- round( (maxx-minv)/(maxv-minv)*1000 ) +1;
			cls[which(cls>1000)]<-1000;
			mc[ i, j] <- cols[  cls ];
		}
	}

	par(mar=c(2,2,2,2)+0.1);
	plot(c(0, nLenArm), c(-0.5, nLenArm), type= "n", xlab="", ylab="", xaxt="n", yaxt="n" );

	for (x in 1:nLenArm)
	for (y in 1:nLenArm)
	{
		if (x==y)
			next;

		x0 <- x;
		y0 <- (nLenArm+1)-y;
		rect(x0-1, y0-1, x0, y0, col=mc[x,y], border="black", angle=-45);
		rect(x0-1+1/8/2, y0-1+1/8/2, x0-1/8/2, y0-1/8/2, col="white", border="black", angle=-45);

		w <- 1/(nLenArm+1)*3/4;
		segments(x0-1+1/8, y0-1+1/8,x0-1+1/8, y0-1+7/8,col="black");
		segments(x0-1+1/8, y0-1+1/8,x0-1+7/8, y0-1+1/8,col="black");

		y05 <- y0-1+1/8+(-log10(0.05)-minv)/(maxv-minv)*3/4
		if ( y05 < y0-1+7/8)
			segments(x0-1+1/8, y0-1+1/8+(-log10(0.05)-minv)/(maxv-minv)*3/4,
					 x0-1+7/8, y0-1+1/8+(-log10(0.05)-minv)/(maxv-minv)*3/4,lty=3, col="gray");
		y001 <- y0-1+1/8+(-log10(0.001)-minv)/(maxv-minv)*3/4
		if ( y001 < y0-1+7/8)
			segments(x0-1+1/8, y0-1+1/8+(-log10(0.001)-minv)/(maxv-minv)*3/4,
				 x0-1+7/8, y0-1+1/8+(-log10(0.001)-minv)/(maxv-minv)*3/4,lty=3, col="gray");


		ysub<-rm[which(dat[,1]==x & dat[,2]==y) ];
		xsub<-dat[which(dat[,1]==x & dat[,2]==y),3 ];
		ysub<-c(ysub,0,0);
		xsub<-c(xsub,x,y);
		points(x0-1+1/8 +xsub*w, y0-1+1/8+ (ysub-minv)/(maxv-minv)*3/4, col=mc[x,y], cex=0.5);
	}

	if (!is.null( labels ))
		for (x in 1:length(labels) )
			text(x-1+0.5, nLenArm-x+0.5, labels[x], adj=c(0.5,0.5), srt=45, cex=0.8)

	title( sTitle );

	l <- nLenArm*3/4;
	x0 <- nLenArm/8;
	for (i in 0:100)
	{
		rect(x0+i*(l/100), -0.25,x0+(i+1)*l/100, -0.35,col=cols[i*10+1], border=cols[i*10+1])
	}

	str_maxv <- sprintf("%8.6f", 10^(-maxv));
	text( x0+100* (l/100), -0.25+2/3*strheight("1"), str_maxv, cex=0.5 );
	text( x0+0* (l/100),   -0.25+2/3*strheight("1"), "1", cex=0.5 );
	x005 <- round( (-log10(0.05) - minv)/(maxv-minv)*100 ) +1;
	if (x005<=101)
	{
		segments(x0+x005*l/100, -0.25, x0+x005*l/100, -0.35, col = "white" );
		text( x0+x005*(l/100), -0.35-strheight("1.0"), ".05", cex=0.5 );

		xs <- c(x0+x005*l/100, x0+x005*l/100-0.02, x0+x005*l/100+0.02);
		ys <- c( -0.25+0.02 , -0.25+0.07, -0.25+0.07);
		polygon(xs, ys, col="gray", border="gray");
	}

	x001 <- round( (-log10(0.01) - minv)/(maxv-minv)*100 ) +1;
	if (x001<=101)
	{
		segments(x0+x001*l/100, -0.25, x0+x001*l/100, -0.35, col = "white");


		text( x0+x001*(l/100), -0.35-strheight("1.0"), ".01", cex=0.5);
	}

	x0001 <- round( (-log10(0.001) - minv)/(maxv-minv)*100 ) +1;
	if (x0001<=101)
	{
		segments(x0+x0001*l/100, -0.25, x0+x0001*l/100, -0.35, col = "white" );
		text( x0+x0001*(l/100), -0.35-strheight("1.0"), ".001", cex=0.5 );

		xs <- c(x0+x0001*l/100, x0+x0001*l/100-0.02, x0+x0001*l/100+0.02);
		ys <- c( -0.25+0.02 , -0.25+0.07, -0.25+0.07);
		polygon(xs, ys, col="gray", border="gray");
	}

	x00001 <- round( (-log10(0.0001) - minv)/(maxv-minv)*100 ) +1;
	if (x00001<=101)
	{
		segments(x0+x00001*l/100, -0.25, x0+x00001*l/100, -0.35, col = "white");
		text( x0+x00001*(l/100), -0.35-strheight("1.0"), ".0001", cex=0.5 );
	}

	x000001 <- round( (-log10(0.00001) - minv)/(maxv-minv)*100 ) +1;
	if (x000001<=101)
	{
		segments(x0+x000001*l/100, -0.25, x0+x000001*l/100, -0.35, col = "white");
		text( x0+x000001*(l/100), -0.35-strheight("1.0"), ".00001", cex=0.5 );
	}
}

#--------------------------------------------------------------
# public: summary(TIAN.dat object)
#
#--------------------------------------------------------------
summary.TIAN.dat<-function( object, output=NULL, append=FALSE )
{
	str <- TIAN.summary_dat(  object );
	if (is.null(output))
		return ( cat( str ) )
	else
		return ( cat( str, file= output, append = append) )
}

#--------------------------------------------------------------
# public: summary(TIAN.ret.1 object)
#
#--------------------------------------------------------------
summary.TIAN.ret.1<-function( object, output=NULL, append=FALSE )
{
	str <- TIAN.summary1(  object, output )
	if (is.null(output))
		return ( cat( str ) )
	else
		return ( cat( str, file= output, append = append) )
}

#--------------------------------------------------------------
# public: summary(TIAN.ret.2 object)
#
#--------------------------------------------------------------
summary.TIAN.ret.2<-function( object, output=NULL, append=FALSE)
{
	str <- TIAN.summary2(  object, output );
	if (is.null(output))
		return ( cat( str ) )
	else
		return ( cat( str, file= output, append = append) )
}

#--------------------------------------------------------------
# public: plot(TIAN.ret.2 object)
#
#--------------------------------------------------------------
plot.TIAN.ret.2<-function( object, pdfFile )
{
		TIAN.draw_correlate2( object, pdfFile  );
}

#--------------------------------------------------------------
# public: summary(TIAN.ret.3 object)
#
#--------------------------------------------------------------
summary.TIAN.ret.3<-function( object, output=NULL, append=FALSE)
{
	str <- TIAN.summary3(  object, output );
	if (is.null(output))
		return ( cat( str ) )
	else
		return ( cat( str, file= output, append = append) )
}

#--------------------------------------------------------------
# public: plot(TIAN.ret.3 object)
#
#--------------------------------------------------------------
plot.TIAN.ret.3<-function( object, pdfFile )
{
		TIAN.draw_correlate2( object, pdfFile  );
}

#--------------------------------------------------------------
# public: full_test( sCsvFile, options )
#
#--------------------------------------------------------------
TIAN.full_test<-function( sCsvFile, options=NULL, output_file=NULL, model=NULL)
{
	dat    <- TIAN.readtable( sCsvFile, options);
	if (is.null(dat))
		stop("Failed to load the csv file.");

	if (is.null(model))
		model <- c(1,2,3);

	if (is.null(output_file))
		output_file <- sCsvFile;

	sumary_file <- paste( output_file, ".txt", sep="");
	summary( dat, output = sumary_file, append=FALSE );

	ret1   <- NULL;
	ret_b1 <- NULL;
	ret_h1 <- NULL;
	ret_f1 <- NULL;
	ret2   <- NULL;
	ret_b2 <- NULL;
	ret_h2 <- NULL;
	ret_f2 <- NULL;
	ret3   <- NULL;
	ret_b3 <- NULL;
	ret_h3 <- NULL;
	ret_f3 <- NULL;

	# 1 SNP analysis
	if (1 %in% model)
	{
		ret1     <- TIAN.snp1( dat);
		summary( ret1,   output = sumary_file, append=TRUE );

		if (include_correction(options, "bonferroni"))
		{
			ret_b1   <- TIAN.adjust_by_bonferroni( ret1 );
			summary( ret_b1, output = sumary_file, append=TRUE );
		}

		if (include_correction(options, "holm"))
		{
			ret_h1   <- TIAN.adjust_by_holm( ret1 );
			summary( ret_h1, output = sumary_file, append=TRUE );
		}

		if (include_correction(options, "fdr"))
		{
			ret_f1   <- TIAN.adjust_by_fdr( ret1 );
			summary( ret_f1, output = sumary_file, append=TRUE );
		}
	}

	# 2 SNP analysis
	if (2 %in% model)
	{
		ret2     <- TIAN.snp2( dat );
		summary( ret2, output = sumary_file, append=TRUE );

		if (include_correction(options, "bonferroni"))
		{
			ret_b2   <- TIAN.adjust_by_bonferroni( ret2 );
			summary( ret_b2, output = sumary_file, append=TRUE );
		}

		if (include_correction(options, "holm"))
		{
			ret_h2   <- TIAN.adjust_by_holm( ret2 );
			summary( ret_h2, output = sumary_file, append=TRUE );
		}

		if (include_correction(options, "fdr"))
		{
			ret_f2   <- TIAN.adjust_by_fdr( ret2 );
			summary( ret_f2, output = sumary_file, append=TRUE );
		}
	}

	# 3 SNP analysis
	if (3 %in% model)
	{
		ret3     <- TIAN.snp3( dat );
		summary( ret3, output = sumary_file, append=TRUE );

		if (include_correction(options, "bonferroni"))
		{
			ret_b3   <- TIAN.adjust_by_bonferroni( ret3 );
			summary( ret_b3, output = sumary_file, append=TRUE );
		}

		if (include_correction(options, "holm"))
		{
			ret_h3   <- TIAN.adjust_by_holm( ret3 );
			summary( ret_h3, output = sumary_file, append=TRUE );
		}

		if (include_correction(options, "fdr"))
		{
			ret_f3   <- TIAN.adjust_by_fdr( ret3 );
			summary( ret_f3, output = sumary_file, append=TRUE );
		}
	}

	return (list(ret1, ret_b1, ret_h1, ret_f1,
				 ret2, ret_b2, ret_h2, ret_f2,
				 ret3, ret_b3, ret_h3, ret_f3 ));
}

include_correction<-function(options, correction)
{
	if ( is.null(options) )
		return (TRUE);

	if ( is.null(options$correction) )
		return (TRUE);

	if ( options$correction=="all" )
		return (TRUE);

	if ( options$correction=="none" || options$correction=="" )
		return (FALSE);

	n <- length( grep(correction, options$correction) )

	return(n>=1);
}


#--------------------------------------------------------------
# public: snp1_test( sCsvFile, options )
#
#--------------------------------------------------------------
TIAN.snp1_test<-function( sCsvFile, options=NULL, output_file=NULL)
{
	r <- TIAN.full_test(sCsvFile, options=options, output_file=output_file, model=c(1));
	return (list( r$ret1, r$ret_b1, r$ret_h1, r$ret_f1 ));
}

#--------------------------------------------------------------
# public: snp2_test( sCsvFile, options )
#
#--------------------------------------------------------------
TIAN.snp2_test<-function( sCsvFile, options=NULL, output_file=NULL)
{
	r <- TIAN.full_test(sCsvFile, options=options, output_file=output_file, model=c(2));
	return (list( r$ret2, r$ret_b2, r$ret_h2, r$ret_f2 ));
}

#--------------------------------------------------------------
# public: snp3_test( sCsvFile, options )
#
#--------------------------------------------------------------
TIAN.snp3_test<-function( sCsvFile, options=NULL, output_file=NULL)
{
	r <- TIAN.full_test(sCsvFile, options=options, output_file=output_file, model=c(3));
	return (list( r$ret3, r$ret_b3, r$ret_h3, r$ret_f3 ));
}
