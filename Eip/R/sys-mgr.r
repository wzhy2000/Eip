#--------------------------------------------------------------
# sys_mgr.base
#
#--------------------------------------------------------------

sys_mgr.base<-function( module )
{
	i_module <- module
	i_seq_id <- 0
	i_start_timer<-proc.time()
	i_hash   <- NULL
	
	get_seq_id<-function()
	{
		i_seq_id <<- i_seq_id+1;
		return (i_seq_id);
	}

	reset<-function()
	{
		i_seq_id <<- 0;
		i_start_timer<<-0;

		i_hash <<- new.env(hash=TRUE, parent=emptyenv(), size=100L);

		assign("plot_doctype",  "pdf", i_hash);
		assign("color.cross",   "yellow", i_hash);
		assign("color.v0",      "blue", i_hash);
		assign("color.v1",      "green", i_hash);
		assign("display.cross", TRUE, i_hash);
		assign("legend.peak"	, 0, i_hash);
	}
	
	set_value<-function( key, value)
	{
		assign( key, value, i_hash);
	}

	get_value<-function( key )
	{
		return ( get( key, i_hash) );
	}

	task_start<-function(...)
	{
		i_start_timer<<-proc.time();
		
		msgs <- list(...);
		if (length(msgs)==0)
			cat("The task is started...\r\n")
		else	
		{
			cat(..., sep="");
		}
		
		flush.console();
	}

	task_elapsed<-function(finished = NA, ...)
	{
		nt <- proc.time() - i_start_timer;
		
		if(is.na(finished))
		{
		}
		else
			sys_t <-  round(c(nt[3], nt[3]*( 1 - finished)/finished ));
		
		sTime1<-"";
		sTime2<-"";
		if (sys_t[1]>60)
			sTime1 <- sprintf( "%02g:%02g:%02g", sys_t[1]%/%3600, (sys_t[1]- sys_t[1]%/%3600*3600)%/%60, sys_t[1]%%60 )
		else
			sTime1 <- sprintf( "%02g seconds", sys_t[1] );

		if (sys_t[2]>60)
			sTime2 <- sprintf( "%02g:%02g:%02g", sys_t[2]%/%3600, (sys_t[2]- sys_t[2]%/%3600*3600)%/%60, sys_t[2]%%60 )
		else
			sTime2 <- sprintf( "%02g seconds", sys_t[2] );
		sTime<-paste( sTime1, " has elapsed, left time: ", sTime2, ". ");
		
		msgs <- list(...);
		if ( length(msgs) != 0 )
		{
			for (i in 1:length(msgs))
				if (msgs[[i]]=="$SYS_PROMPT$")
					msgs[[i]] <- sTime;
		}
		else
			msgs[[1]]<-paste( sTime, "\r\n");
				
		cat( unlist(msgs), sep="" );
		flush.console();
	}

	task_stop<-function(...)
	{
		msgs <- list(...);
		if (length(msgs)==0)
			cat("The task is stopped...\r\n")
		else	
			cat( ... , sep="" );

		flush.console();
	}

	return (list(reset=reset, seq_id=get_seq_id, task_start=task_start, task_stop=task_stop, task_elapsed=task_elapsed, set_value=set_value, get_value=get_value))
	
	
}
