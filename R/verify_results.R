# Base R verification and supplementary summaries of completed full runs.
# These functions never fit a model or modify the input predictions.

wrr_verify_results <- function(root, out_dir = file.path("analysis_output","verification")) {
  root <- normalizePath(root, winslash="/", mustWork=TRUE)
  f <- function(p, n) read.csv(file.path(p,"tables",n), check.names=FALSE)
  status <- readLines(file.path(root,"RUN_STATUS.txt"),warn=FALSE)
  stopifnot(any(grepl("^COMPLETED full",status)))
  cfg <- readRDS(file.path(root,"configuration.rds"))
  stopifnot(is.list(cfg), identical(cfg$mode,"full"), isTRUE(cfg$run_robustness),
            cfg$scenarios==10000L, cfg$estimation_B==199L,
            cfg$margin_window==750L, cfg$q==.10, cfg$k==3L)
  u <- f(root,"marginal_PIT.csv")
  stopifnot(nrow(u)==3273L, sum(u$split=="test")==502L)
  test <- u$date[u$split=="test"]
  eq <- function(a,b,label) {
    if(!isTRUE(all.equal(as.numeric(a),as.numeric(b),tolerance=1e-8,check.attributes=FALSE)))
      stop("Tolerance failure: ",label,call.=FALSE)
  }
  hac <- function(d,lag=10L) {
    n <- length(d); z <- d-mean(d)
    v <- sum(z*z)/n
    for(k in seq_len(min(lag,n-2L))) v <- v+2*(1-k/(lag+1))*sum(z[(k+1):n]*z[1:(n-k)])/n
    se <- sqrt(max(0,v)/n)
    c(se=se,p=if(se>0)2*pnorm(-abs(mean(d)/se)) else if(mean(d)==0)1 else 0)
  }
  specs <- c("primary","empirical_margin","train_winsor")
  universes <- c("main_no_vix","extended_with_vix")
  metrics <- c("nll","energy_common","variogram_common","brier","quantile_loss")
  ss <- pp <- local <- pit <- list()
  for(spec in specs) {
    base <- if(spec=="primary")root else file.path(root,"sensitivity",spec)
    pc <- f(base,"PIT_calibration.csv")
    deg <- pc$hit_rate %in% c(0,1)
    pc$degenerate_hit_series <- deg
    pc$hac_p_documented <- ifelse(deg,NA_real_,pc$calibration_hac_p)
    pc$holm_p_documented <- NA_real_
    for(split in unique(pc$split)) {
      i <- which(pc$split==split);p <- pc$hac_p_documented[i]
      nc <- p.adjust(ifelse(is.na(p),1,p),"holm",n=length(i))
      nc[is.na(p)] <- NA_real_;pc$holm_p_documented[i] <- nc
    }
    pit[[spec]] <- cbind(spec=spec,pc)
    for(uni in universes) {
      b <- file.path(base,uni);d <- f(b,"daily_forecasts.csv")
      s <- f(b,"model_comparison.csv");di <- f(b,"paired_inference.csv")
      stopifnot(nrow(s)==if(spec=="primary")20L else 8L,
                !anyDuplicated(d[c("date","model")]),all(table(d$model)==502L),
                all(is.finite(as.matrix(d[sapply(d,is.numeric)]))))
      for(m in unique(d$model)) {
        dm <- d[d$model==m,]
        stopifnot(identical(sort(dm$date),test))
        sr <- s[s$model==m,]
        for(metric in intersect(names(s),names(d))[sapply(s[intersect(names(s),names(d))],is.numeric)]) {
          eq(mean(dm[[metric]]),sr[[metric]],paste(spec,uni,m,metric))
        }
        eq(-mean(dm$nll),sr$mean_log_score,"mean log score")
        eq((dm$probability-dm$observed)^2,dm$brier,"daily Brier")
      }
      ref <- d[d$model=="WRR_direct",]
      for(m in setdiff(unique(d$model),"WRR_direct")) {
        other <- d[d$model==m,];j <- match(ref$date,other$date)
        for(metric in metrics) {
          loss <- ref[[metric]]-other[[metric]][j];h <- hac(loss)
          i <- which(di$benchmark==m & di$metric==metric)
          stopifnot(length(i)==3L)
          eq(rep(mean(loss),3),di$difference[i],"paired difference")
          eq(rep(h["se"],3),di$hac_se[i],"HAC standard error")
          eq(rep(h["p"],3),di$p_two_sided[i],"HAC p")
        }
      }
      for(block in unique(di$block)) {
        i <- which(di$block==block)
        eq(p.adjust(di$p_two_sided[i],"holm"),di$p_holm[i],"Holm")
      }
      ss[[paste(spec,uni)]] <- cbind(spec=spec,universe=uni,s)
      pp[[paste(spec,uni)]] <- cbind(spec=spec,universe=uni,di[di$block==20L,])
      if(spec=="primary") {
        z <- f(b,"local_estimation_bootstrap_draws.csv")
        it <- f(b,"local_estimation_intervals.csv")
        stopifnot(nrow(z)==1791L,nrow(it)==9L,all(it$B==199L),all(it$success==199L),all(is.finite(z$probability)))
        for(i in seq_len(nrow(it))) {
          zz <- z[z$date==it$date[i],]
          eq(quantile(zz$probability,c(.025,.975),names=FALSE),c(it$probability_low[i],it$probability_high[i]),"local interval")
        }
        local[[uni]] <- cbind(universe=uni,it)
      }
      cat("PASS ",spec," / ",uni," (",nrow(d)," forecast rows)\n",sep="")
    }
  }
  dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)
  write.csv(do.call(rbind,ss),file.path(out_dir,"verified_model_comparison.csv"),row.names=FALSE)
  write.csv(do.call(rbind,pp),file.path(out_dir,"verified_paired_block20.csv"),row.names=FALSE)
  write.csv(do.call(rbind,local),file.path(out_dir,"verified_local_intervals.csv"),row.names=FALSE)
  write.csv(do.call(rbind,pit),file.path(out_dir,"PIT_calibration_reporting.csv"),row.names=FALSE)
  cat("ALL SIX CASES VERIFIED. Original result files unchanged.\n")
  invisible(TRUE)
}

wrr_supplementary <- function(root, output_directory=file.path("analysis_output","supplementary_R")) {
  stopifnot(any(grepl("^COMPLETED full",readLines(file.path(root,"RUN_STATUS.txt"),warn=FALSE))))
  dir.create(output_directory,recursive=TRUE,showWarnings=FALSE)
  read_table <- function(p,f)read.csv(file.path(p,"tables",f),check.names=FALSE)
  u <- read_table(root,"marginal_PIT.csv")
  common <- c("XLE","ICLN","TAN","XLF","GLD","UUP")
  K <- rowSums(u[common]<=.10)
  A <- as.integer(K>=3L)
  wt <- cond <- loc <- state <- list()
  for(uni in c("main_no_vix","extended_with_vix")) {
    dd <- if(uni=="main_no_vix")6L else 7L
    b <- file.path(root,uni)
    rg <- read_table(b,"regimes.csv");v<-rg$lagged_vix;s<-rg$regime
    stopifnot(identical(u$date,rg$date))
    diag <- read_table(b,"refit_weight_diagnostics.csv")
    x <- diag[diag$model=="WRR_direct",]
    z <- diag[diag$model=="WRR_no_tail",]
    stopifnot(identical(x$date,z$date))
    for(h in seq_len(nrow(x))) {
      t <- x$t[h];W <- x$n_window[h];j <- seq.int(t-W,t-1L)
      w <- exp(-2*(t-1L-j)/W-1.25*abs(v[j]-v[t])/sd(v[j]))*ifelse(s[j]==s[t],3,1)
      p0 <- w/sum(w);p1 <- p0*(1+2*A[j]);p1 <- p1/sum(p1)
      p0f <- (1-z$uniform_mix[h])*p0+z$uniform_mix[h]/W
      p1f <- (1-x$uniform_mix[h])*p1+x$uniform_mix[h]/W
      e0 <- sum(p0*A[j]);e1<-sum(p1*A[j]);k0<-sum(p0*K[j])/6;k1<-sum(p1*K[j])/6
      stopifnot(abs(e1-3*e0/(1+2*e0))<1e-12,
                abs(k1-k0-2*(sum(p0*A[j]*K[j])-e0*sum(p0*K[j]))/(1+2*e0)/6)<1e-12,
                abs(x$ess[h]-1/sum(p1f^2))<1e-6,
                abs(z$ess[h]-1/sum(p0f^2))<1e-6)
      wt[[length(wt)+1L]] <- data.frame(dimension=dd,date=x$date[h],regime=s[t],window=W,
        event_mass_base=e0,event_mass_tilted=e1,mean_marginal_hit_base=k0,mean_marginal_hit_tilted=k1,
        mean_marginal_hit_increase=k1-k0,event_mass_no_tail_after_ess=sum(p0f*A[j]),
        event_mass_wrr_after_ess=sum(p1f*A[j]),mean_marginal_hit_no_tail_after_ess=sum(p0f*K[j])/6,
        mean_marginal_hit_wrr_after_ess=sum(p1f*K[j])/6,wrr_ess=x$ess[h],no_tail_ess=z$ess[h],
        wrr_uniform_mix=x$uniform_mix[h],no_tail_uniform_mix=z$uniform_mix[h])
    }
    daily <- read_table(b,"daily_forecasts.csv")
    for(m in c("WRR_direct","rolling_matched","WRR_no_tail","WRR_resampled_B20","dcc_gaussian"))
      for(ss in c("calm","medium","stress")) {
        q <- daily[daily$model==m & daily$regime==ss,]
        cond[[length(cond)+1L]]<-data.frame(dimension=dd,model=m,state=ss,n=nrow(q),events=sum(q$observed),
          observed=mean(q$observed),probability=mean(q$probability),brier=mean(q$brier),
          adjusted_brier=mean(q$brier_mc_adjusted),mean_log=-mean(q$nll))
      }
    ref <- daily[daily$model=="WRR_direct",];ref<-ref[order(ref$date),]
    ma <- daily[daily$model=="rolling_matched",];ma<-ma[match(ref$date,ma$date),]
    j <- findInterval(match(ref$date,rg$date),x$t)
    stopifnot(all(j>0L))
    last <- s[x$t[j]]
    for(ss in c("calm","medium","stress"))for(ls in c("calm","medium","stress")) {
      i <- ref$regime==ss & last==ls
      if(any(i)) state[[length(state)+1L]]<-data.frame(dimension=dd,state=ss,refit_state=ls,n=sum(i),
        mean_logscore_difference=mean(ma$nll[i]-ref$nll[i]),events=sum(ref$observed[i]))
    }
    bs <- read_table(b,"local_estimation_bootstrap_draws.csv")
    for(dt in unique(bs$date)) {
      q <- bs[bs$date==dt,];stopifnot(nrow(q)==199L,all(is.finite(q$log_score)))
      l <- quantile(q$log_score,c(.025,.975),names=FALSE);p <- quantile(q$probability,c(.025,.975),names=FALSE)
      loc[[length(loc)+1L]] <- data.frame(dimension=dd,date=dt,regime=s[match(dt,rg$date)],
        log_score_mean=mean(q$log_score),log_score_sd=sd(q$log_score),log_score_low=l[1L],log_score_high=l[2L],
        log_score_width=diff(l),probability_mean=mean(q$probability),probability_sd=sd(q$probability),
        probability_low=p[1L],probability_high=p[2L])
    }
  }
  frames<-list(tail_tilt_mass=wt,state_event_forecasts=cond,local_logscore_uncertainty=loc,forecast_refit_state_alignment=state)
  for(nm in names(frames)) {
    q <- do.call(rbind,frames[[nm]])
    write.csv(q,file.path(output_directory,paste0(nm,".csv")),row.names=FALSE)
  }
  cat("PASS: tail-event tilt identities, archived ESS and supplementary forecast diagnostics.\n")
  invisible(TRUE)
}
