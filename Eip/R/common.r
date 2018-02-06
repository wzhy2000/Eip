#--------------------------------------------------------------
# COMMON.plot_to_doc
#
#--------------------------------------------------------------
COMMON.plot_to_doc<-function( filename, doctype, h=800, w=800)
{
	if (doctype=="pdf")
		pdf( paste(filename, sys$seq_id(), "pdf", sep=".") )
	else if (doctype=="jpg")
		jpeg( paste(filename, sys$seq_id(), "jpg", sep=".") )
	else if (doctype=="wmf")
		win.metafile( paste(filename, sys$seq_id(), "wmf", sep=".") )
	else if (doctype=="ps")
		postscript( paste(filename, sys$seq_id(), "ps", sep=".") )
	else 
		png( paste(filename, sys$seq_id(), "png", sep="."), h=h, w=w);
}

#--------------------------------------------------------------
# COMMON.plot_one_chromosome
#
#--------------------------------------------------------------
COMMON.plot_one_chromosome<-function(chr_logs, marker_list=NA, threshold=NA, p05=NA, p01=NA)
{
	
	par(mar=c(0.2, 0.2, 0.2, 0.2) );
	par(bty="n");
	plot( chr_logs[,2], chr_logs[,3], xlim=c( -10, 90), ylim=c(-25, 75), 
			type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i");

	#========Draw LR Curve in a range(-10, 0, 90, 75)========
	# Border(0,0,85,70);
	rect(0, 0, 85, 70, col = c(NA,"midnightblue"));
	#Draw Curve in (0,0,85,70)
	max_xlim <- max( chr_logs[,2] );
	min_xlim <- min( chr_logs[,2] );
	{
		max_ylog <- max( chr_logs[,3] )*1.1;
		min_ylog <- min( chr_logs[,3] );
		if (!any(is.na(marker_list))  )
		{
			max_xlim <- max( marker_list[,2] );
			min_xlim <- min( marker_list[,2] );
		}

		lines( chr_logs[,2]/(max_xlim - min_xlim)*85, chr_logs[,3]/(max_ylog - min_ylog)*70, lwd=2, pch=19);
		points( chr_logs[,2]/(max_xlim - min_xlim)*85, chr_logs[,3]/(max_ylog - min_ylog)*70,type = "p", pch=21);


		if ( !any(is.na(marker_list)))
			for (j in 2:length(marker_list[,1])-1 )
			{
				x0 <- marker_list[j,2]/(max_xlim - min_xlim)*85;
				segments( x0, 0, x0, 70, lty="dotted");
			}
		
		if (!is.na(p05))
			segments( 0, p05/(max_ylog - min_ylog)*70, 85, p05/(max_ylog - min_ylog)*70, lty="solid", lwd=1);
		
		if (!is.na(p01))
			segments( 0, p01/(max_ylog - min_ylog)*70, 85, p01/(max_ylog - min_ylog)*70, lty="solid", lwd=1);
			
		
		x_temp <- (max_xlim - min_xlim)/6;
		x_temp <- 10*ceiling(x_temp/10);
		
		for (j in 1:(max_xlim%/%x_temp) )
		{
			if( j*x_temp>=min_xlim && j*x_temp<=max_xlim)
			{
				segments( j*x_temp/(max_xlim - min_xlim)*85, 0, j*x_temp/(max_xlim - min_xlim)*85, 1 );
				segments( j*x_temp/(max_xlim - min_xlim)*85, 69, j*x_temp/(max_xlim - min_xlim)*85, 70 );
				text( j*x_temp/(max_xlim - min_xlim)*85, -0.3, (j*x_temp), adj=c(0.5, 1) );
			}
			
		}
		text( min_xlim/(max_xlim - min_xlim)*85, -0.3, "0", adj=c(1, 0.5) );
		
		
		y_temp <- (max_ylog - min_ylog)/12;
		y_temp <- 2*ceiling(y_temp/2);
		for (j in 1:((max_ylog-min_ylog)%/%y_temp) )
		{
			if( j*y_temp >= min_ylog && j*y_temp <= max_ylog )
			{
				segments( 0, j*y_temp/(max_ylog - min_ylog)*70, 1, j*y_temp/(max_ylog - min_ylog)*70, 1 );
				segments( 85, j*y_temp/(max_ylog - min_ylog)*70, 84, j*y_temp/(max_ylog - min_ylog)*70, 1 );
				text( 0, j*y_temp/(max_ylog - min_ylog)*70, (j*y_temp), adj=c(1, 0.5) );
			}
			
		}
	}
	

	if ( any(is.na(marker_list)))
	{
		for (j in 2:length(chr_logs[,2]) )
		{
			segments(chr_logs[j,2], min_ylog, chr_logs[j,2], min_ylog- sticker_h*3);
		}
   	}
	
	#======= Draw marker ruler in a range(-10,-25,90,25)=======
	segments( 0, -17, 85, -17, lwd=2, col="black");
	for (j in 1:length(marker_list[,1]) )
	{
		x0 <- marker_list[j,2]/(max_xlim - min_xlim)*85;
		segments( x0, -17, x0, -16);
		text( x0, -16, marker_list[j,3], adj=c(0,0.5), srt=90, cex=0.8);
	}
	
	for (j in 1:( length(marker_list[,1])-1) )
	{
		x0 <- 0.5*(marker_list[j+1,2] + marker_list[j,2])/(max_xlim - min_xlim)*85;
		segments( x0, -17, x0, -18);
		
		text( x0, -18, (marker_list[j+1,2]-marker_list[j,2]), adj=c(1,0.5), srt=90,cex=0.8);
	}

	text( -5, 35, "LR", adj=c(1,0.5), srt=90,cex=1.2);
}


#--------------------------------------------------------------
# COMMON.plot_chromosom_map
#
#--------------------------------------------------------------
COMMON.plot_chromosome_map<-function(chr_nums, level_cnt, chr_logs, marker_list=NA )
{
	max_log <- 0;
	min_log <- 0;
	ch_ev <- matrix(0, nrow=chr_nums, ncol=3);
	
	for ( i in  1:length(chr_logs[,1]) )
	{
		n <- chr_logs[i,1];
	    	if ( chr_logs[i,2]>ch_ev[n,1] )
			ch_ev[n,1] <- chr_logs[i,2];

	    	if ( chr_logs[i,3]>ch_ev[n,2] )
			ch_ev[n,2] <- chr_logs[i,3];

	    	if (chr_logs[i,3]>max_log)
			max_log <- chr_logs[i,3];

	    	if (chr_logs[i,3]<min_log)
			min_log <- chr_logs[i,3];

	    	ch_ev[n,3] <- ch_ev[n,3] + 1;
	}

	nRowChrs <- array(1, level_cnt);
	nRowChrs[1]<- chr_nums - (level_cnt-1)
	delt_sd=100000;
	delt_old_sd<-0;
	if (chr_nums>=2)
	{
		while(delt_sd>10)
		{
			for(i in 1:(level_cnt-1))
			{
				dm<-get_SmallestVar2( ch_ev, nRowChrs, i );					
				nRowChrs[i]<-dm[2];
				nRowChrs[i+1]<-dm[3];
			}
		
			dm<-get_deltvar( ch_ev, nRowChrs, 1 , level_cnt );
			delt_sd <- abs( dm[1]-delt_old_sd );
			delt_old_sd<-dm[1];
		}
	}
	else
		dm<-c(0, ch_ev[1,1], 0 );

	pos <- matrix(0, nrow=chr_nums, ncol=4);
	x0 <- 0;
	h  <- 1/level_cnt;
	y0 <- 1-h;

	nChr <- 1;
	for (i in 1:level_cnt)
	{
   		for ( j in 1: nRowChrs[i] )
    		{
    			pos[nChr,1] <- x0;
    			pos[nChr,2] <- y0;
    			pos[nChr,3] <- ch_ev[nChr,1]/dm[i+1];
    			pos[nChr,4] <- h;
    			x0 <- x0 +pos[nChr,3]; 
    			nChr <- nChr+1;
   		} 
   		x0 <- 0;
   		y0 <- y0-h;
	}


	plot.new();	
	title(xlab="", ylab='  -2 Log Likelihood Ratio');	

	for (i in 1:chr_nums )
	{
		sub_plot<-c(pos[i,1], pos[i,1]+pos[i,3], pos[i,2], pos[i,2]+pos[i,4])*0.9+0.05
		par(fig=sub_plot, new=TRUE)
		par(mai=c(0,0,0,0));
		
		ll <- matrix(0, nrow=ch_ev[i,3], ncol=2 );
		n <-1;
		for (j in 1:length(chr_logs[,1]) )
		{
			if ( chr_logs[j,1] == i )
			{
				ll[n,1] <- chr_logs[j,2];
				ll[n,2] <- chr_logs[j,3];
				n <- n+1;
			}
		}

		plot( ll[,1], ll[,2], xlim=c(0, ch_ev[i,1]), ylim=c(min_log*0.9, max_log*1.3), 
			type="n", main="", xlab="", ylab="", xaxt="n", yaxt="n",xaxs="i", yaxs="i");

		lines(ll[,1],ll[,2],lwd=1);
		
		#marker
		sticker_h<- ( max_log*1.3- min_log*0.9)/20;
		if ( any(is.na(marker_list)))
		{
			for (j in 2:length(ll[,1]) )
			{
				segments(ll[j,1], min_log*0.9, ll[j,1], min_log*0.9+ sticker_h);
			}
    		}
		else
		{
			for (j in 1:length(marker_list[,1]) )
			{
				if (marker_list[j,1]==i)
				{	
					x0<-marker_list[j,2];
					segments( x0, min_log*0.9, x0, min_log*0.9+ sticker_h);
				}
			}
		}

		#chrom no 
		text(5, max_log*1.2, paste("",i) , font=4);

	}
	
	return(pos);
}


get_SmallestVar2<-function( ch_ev, nRowChrs, iStart)
{
	nSum<-sum(nRowChrs[c(iStart,iStart+1)])
	nRowBack <- nRowChrs;
	x1<-1;
	x2<-nSum - 1
	nRowBack[iStart]<-x1;	
	nRowBack[iStart+1]<- x2;	

	delt_var <- get_deltvar(ch_ev, nRowBack, iStart , 2);
	
	for ( i in 2:nSum-1 )
    	{
		nRowBack[iStart]   <- i;	
		nRowBack[iStart+1] <- nSum-i;	
		var2 <- get_deltvar(ch_ev, nRowBack, iStart , 2);
		if (var2[1]<delt_var)
		{
			delt_var<-var2[1];
			x1<-i;
			x2<-nSum-i;
		}
	}
	
	return( c(delt_var, x1, x2))
}

get_deltvar<-function( ch_ev, nRowChrs, iStart , iCount )
{
	nWidths <- array(0, iCount );

	index<-1;
	ns<-1;

	for ( i in 1:length(nRowChrs) )
	{
		ne <- nRowChrs[i]+ns-1;

		if ( (i>=iStart) && i<(iStart+iCount) )
		{
			nWidths[index] <- sum( ch_ev[ c(ns:ne), 1] );
			index<-index+1;
		}
		
		ns <- ns + nRowChrs[i];
	}
	
	return( c( var(nWidths), nWidths ) );
}

trim_space<-function( x )
{
	if (is.na(x))
		return (NA);
		
	x <- as.character(x);
	while (substr(x,1,1)==" " )
	{
		n<-nchar(x);
		x<- substr(x,2, n);
	}
	
	n<-nchar(x);
	while (substr(x,n,n)==" " )
	{
		x<- substr(x,1, n-1);
		n<-nchar(x);
	}		
	
	return(x);
}
