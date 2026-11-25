# Weighted Rolling Regime-Adaptive Vine Copulas (WRR-TVine)
# Forecasting engine. Sourcing defines functions; use the run scripts to fit models.
# All data transformations, estimation settings and stochastic seeds are explicit.

wrr_config <- function(mode = c("full", "smoke")) {
  mode <- match.arg(mode)
  list(mode = mode, seed = 20261113L,
       start = "2010-01-01", end = "2025-12-31", train_end = "2021-12-31",
       valid_end = "2023-12-31", test_start = "2024-01-01",
       symbols = c("XLE", "ICLN", "TAN", "XLF", "GLD", "UUP", "VIX"),
       common_assets = c("XLE", "ICLN", "TAN", "XLF", "GLD", "UUP"),
       out_dir = file.path("results", mode), prices_csv = "data/prices.csv",
       download = TRUE, resume = TRUE, preprocessing = "none", winsor_p = 0.001,
       margin = "garch_select", margin_window = 750L,
       margin_candidates = c("sGARCH", "gjrGARCH", "eGARCH"),
       arma_maxit = 2000L, garch_retry = TRUE,
       progress_refits = if(mode == "full") 5L else 1L,
       regime_probs = c(.60, .85),
       regime_sensitivity = list(c(.55, .80), c(.65, .90)),
       grid = if (mode == "full") expand.grid(window = c(500L, 750L, 1000L),
                   refit = c(10L, 20L), family = c("tail", "general"),
                   stringsAsFactors = FALSE) else
         data.frame(window = 120L, refit = 20L, family = "tail"),
       recency_decay = 2, same_regime_bonus = 3, vix_distance_decay = 1.25,
       tail_bonus = 2, q = .10, k = 3L, ess_min_fraction = .15,
       min_regime_n = if (mode == "full") 150L else 40L,
       family_general = c(0L, 1L, 2L, 3L, 4L, 5L, 6L),
       family_tail = c(0L, 2L, 3L, 4L),
       trunclevel = NA_integer_, refit_on_regime_change = FALSE,
       resample_n = if (mode == "full") 750L else 100L,
       ensemble_B = if (mode == "full") c(1L, 5L, 20L) else c(1L, 2L),
       scenarios = if (mode == "full") 10000L else 300L,
       scenario_sensitivity = if (mode == "full") c(2000L, 5000L, 20000L) else c(200L, 400L),
       mc_repeats = if (mode == "full") 10L else 2L,
       audit_origins = if (mode == "full") 9L else 2L,
       estimation_B = if (mode == "full") 199L else 3L,
       estimation_block = 20L, estimation_scenarios = if (mode == "full") 20000L else 300L,
       covar_scenarios = if (mode == "full") 100000L else 2000L,
       covar_repeats = if (mode == "full") 5L else 2L,
       covar_min_joint_tail = if (mode == "full") 100L else 2L,
       inference_B = if (mode == "full") 1999L else 49L,
       block_lengths = c(10L, 20L, 40L), hac_lag = 10L,
       dcc_grid = subset(expand.grid(a = c(.01,.03,.05,.07), b = c(.85,.90,.94,.97)), a+b < .995),
       run_robustness = mode == "full", run_audits = TRUE,
       risk_alpha = .05, max_test = if (mode == "full") Inf else 8L)
}

wrr_check_packages <- function(cfg) {
  req <- c("VineCopula", "copula")
  if (cfg$download && !file.exists(cfg$prices_csv)) req <- c(req, "quantmod")
  if (cfg$margin == "garch_select") req <- c(req, "rugarch")
  missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages first: ", paste(missing, collapse = ", "))
  if (!("weights" %in% names(formals(VineCopula::RVineStructureSelect))))
    stop("Installed VineCopula lacks observation-weight support.")
  invisible(TRUE)
}

wrr_seed <- function(...) {
  bytes <- utf8ToInt(paste(..., collapse = "|")); h <- 104729
  for (v in bytes) h <- (h * 131 + v) %% 2147483629
  as.integer(h + 1)
}
wrr_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had) assign(".Random.seed", old, envir = .GlobalEnv) else
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed); force(expr)
}
wrr_clip <- function(u, eps = 1e-8) pmin(pmax(u, eps), 1-eps)
wrr_hash <- function(x) {
  # During a run, use a stable local cache directory rather than an OS-cleaned
  # temporary directory. The serialization and hash format remain unchanged.
  d <- getOption("wrr.hash_dir", file.path(tempdir(), "wrr_hash"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(d)) stop("Could not create the checkpoint hash directory: ", d)
  f <- tempfile(pattern = "hash_", tmpdir = d)
  on.exit(unlink(f), add = TRUE)
  saveRDS(x, f, version = 3)
  h <- unname(tools::md5sum(f))
  if (length(h) != 1L || is.na(h)) stop("Could not hash checkpoint input")
  h
}
wrr_cache <- function(path, key, compute, resume = TRUE) {
  if (resume && file.exists(path)) {
    z <- tryCatch(readRDS(path), error = function(e) {
      wrr_note("Unreadable checkpoint; recomputing: ", path, " (", conditionMessage(e), ")")
      NULL
    })
    if (is.list(z) && "value" %in% names(z) && identical(z$key, key)) return(z$value)
  }
  value <- compute(); dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(basename(path), "_"), tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(list(key = key, value = value), tmp)
  # On systems where rename cannot replace an existing file, use a checked copy.
  if (!suppressWarnings(file.rename(tmp, path))) {
    if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not write checkpoint: ", path)
  }
  value
}
wrr_csv <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, file, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
wrr_pd <- function(R) {
  R <- (R+t(R))/2
  if (any(!is.finite(R))) stop("Non-finite correlation matrix")
  e <- eigen(R, symmetric = TRUE)
  R <- cov2cor(e$vectors %*% diag(pmax(e$values,1e-7), nrow(R)) %*% t(e$vectors))
  diag(R) <- 1; R
}
wrr_split <- function(dates, cfg) {
  out <- ifelse(dates <= as.Date(cfg$train_end), "train",
                ifelse(dates <= as.Date(cfg$valid_end), "validation", "test"))
  if (as.Date(cfg$test_start) != as.Date(cfg$valid_end)+1)
    stop("test_start must immediately follow valid_end")
  if (any(table(factor(out, levels = c("train", "validation", "test"))) == 0))
    stop("Empty chronological split")
  out
}
wrr_prices <- function(cfg) {
  download_warnings<-list();download_checks<-list()
  if(file.exists(cfg$prices_csv)) {
    wrr_note("Using frozen prices: ",normalizePath(cfg$prices_csv,winslash="/"))
    p<-read.csv(cfg$prices_csv,check.names=FALSE)
    if(!("date"%in%names(p)))names(p)[1]<-"date"
    names(p)[names(p)=="^VIX"]<-"VIX"
  } else {
    if(!cfg$download)stop("prices.csv not found and download disabled")
    xs<-lapply(cfg$symbols,function(a) {
      sy<-if(a=="VIX")"^VIX" else a;wrr_note("Downloading ",sy)
      ans<-wrr_capture(function()quantmod::getSymbols(sy,src="yahoo",from=cfg$start,
        to=as.character(as.Date(cfg$end)+1),auto.assign=FALSE))
      if(nrow(ans$warnings))download_warnings[[a]]<<-cbind(asset=a,ans$warnings)
      if(!is.null(ans$error)) {
        wrr_csv(wrr_bind(download_warnings,data.frame(asset=character(),message=character(),call=character())),
          file.path(cfg$out_dir,"tables/download_warnings.csv"))
        stop("Download failed for ",a,": ",conditionMessage(ans$error))
      }
      x<-ans$value;y<-if(a=="VIX")quantmod::Cl(x) else quantmod::Ad(x)
      colnames(y)<-a
      download_checks[[a]]<<-data.frame(asset=a,download_rows=NROW(x),
        missing_all_fields=sum(is.na(x)),missing_selected_price=sum(is.na(y)),
        nonpositive_selected_price=sum(as.numeric(y)<=0,na.rm=TRUE))
      if(nrow(ans$warnings))wrr_note(a,": ",nrow(ans$warnings),
        " download warning(s) recorded; missing selected prices=",sum(is.na(y)))
      y
    })
    x<-do.call(merge,c(xs,all=FALSE))
    p<-data.frame(date=as.Date(zoo::index(x)),as.matrix(x),check.names=FALSE)
    wrr_csv(p,cfg$prices_csv)
  }
  # Preserve earlier download diagnostics when reading the same frozen CSV again.
  if(length(download_checks)) {
    wrr_csv(do.call(rbind,download_checks),file.path(cfg$out_dir,"tables/download_diagnostics.csv"))
    wrr_csv(wrr_bind(download_warnings,data.frame(asset=character(),message=character(),call=character())),
      file.path(cfg$out_dir,"tables/download_warnings.csv"))
  }
  if(!all(c("date",cfg$symbols)%in%names(p)))stop("Missing required price columns")
  p$date<-as.Date(p$date)
  if(anyNA(p$date)||anyDuplicated(p$date))stop("Invalid or duplicate dates in price input")
  if(!all(vapply(p[cfg$symbols],is.numeric,logical(1))))stop("Price columns must be numeric")
  p<-p[p$date>=as.Date(cfg$start)&p$date<=as.Date(cfg$end),c("date",cfg$symbols)]
  p<-p[order(p$date),];P<-as.matrix(p[cfg$symbols])
  bad<-!is.finite(P)|P<=0;ok<-rowSums(bad)==0L
  dropped<-p[!ok,,drop=FALSE]
  dropped$invalid_assets<-if(any(!ok))apply(bad[!ok,,drop=FALSE],1,function(z)paste(cfg$symbols[z],collapse=";")) else character()
  wrr_csv(dropped,file.path(cfg$out_dir,"tables/dropped_prices.csv"))
  qc<-data.frame(asset=cfg$symbols,input_rows=nrow(P),missing=sum(is.na(P)))
  qc$missing<-colSums(is.na(P));qc$nonfinite<-colSums(!is.finite(P));qc$nonpositive<-colSums(P<=0,na.rm=TRUE)
  wrr_csv(qc,file.path(cfg$out_dir,"tables/price_quality.csv"))
  wrr_note("Price alignment: ",sum(ok)," retained; ",sum(!ok)," removed; no interpolation")
  p<-p[ok,];if(nrow(p)<150)stop("Too few complete price rows")
  wrr_csv(p,file.path(cfg$out_dir,"tables/prices_snapshot.csv"))
  gaps<-data.frame(previous_date=head(p$date,-1),date=tail(p$date,-1),calendar_days=as.integer(diff(p$date)))
  wrr_csv(gaps[gaps$calendar_days>4,,drop=FALSE],file.path(cfg$out_dir,"tables/calendar_gaps_over4_days.csv"))
  p
}

wrr_data <- function(p, cfg) {
  R <- apply(as.matrix(p[cfg$symbols]), 2, function(x) diff(log(x)))
  dates <- p$date[-1]; split <- wrr_split(dates, cfg); raw <- R
  bounds <- data.frame(asset=colnames(R), lower=-Inf, upper=Inf)
  if (cfg$preprocessing == "train_winsor") {
    for (j in seq_len(ncol(R))) {
      b <- quantile(R[split=="train",j], c(cfg$winsor_p,1-cfg$winsor_p), type=8)
      bounds[j,2:3] <- b; R[,j] <- pmax(b[1],pmin(b[2],R[,j]))
    }
  } else if (cfg$preprocessing != "none") stop("Unknown preprocessing option")
  # v[t] is the observed VIX level at the previous common trading date.
  list(R=R, raw=raw, dates=dates, v=p$VIX[-nrow(p)], split=split, bounds=bounds)
}
wrr_empirical_margin <- function(x, window) {
  u <- rep(NA_real_,length(x))
  for (t in seq_along(x)) if (t>window) {
    h <- x[(t-window):(t-1)]
    u[t] <- (sum(h<x[t])+.5*sum(h==x[t])+.5)/(length(h)+1)
  }
  list(u=wrr_clip(u), z=qnorm(wrr_clip(u)), method="rolling_empirical", x=x, window=window)
}
wrr_note <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse=""))
  cat(line, "\n", sep=""); flush.console()
  con <- getOption("wrr.log_connection", NULL)
  if(!is.null(con)) try(flush(con), silent=TRUE)
  progress <- getOption("wrr.progress_file", NULL)
  if(!is.null(progress)) writeLines(line,progress)
  invisible(line)
}
wrr_capture <- function(fun) {
  rows <- list(); error <- NULL
  value <- tryCatch(withCallingHandlers(fun(), warning=function(w) {
    rows[[length(rows)+1L]] <<- data.frame(message=conditionMessage(w),
      call=paste(deparse(conditionCall(w)),collapse=" "),stringsAsFactors=FALSE)
    invokeRestart("muffleWarning")
  }),error=function(e) {error <<- e; NULL})
  warnings <- if(length(rows))do.call(rbind,rows) else
    data.frame(message=character(),call=character())
  list(value=value,error=error,warnings=warnings)
}
wrr_error_text <- function(x) if(is.null(x)) "" else conditionMessage(x)
wrr_bind <- function(rows, empty) {
  rows <- Filter(function(x)!is.null(x) && nrow(x)>0,rows)
  if(length(rows))do.call(rbind,rows) else empty
}
wrr_arma_start <- function(x,maxit=2000L) {
  # Training data only. These are optimization starts, never fixed coefficients.
  sc<-sd(x);if(!is.finite(sc)||sc<=0)stop("Degenerate training returns")
  rows<-warnings<-list();start<-NULL
  for(method in c("CSS","ML")) {
    ans<-wrr_capture(function()stats::arima(x/sc,order=c(1L,0L,1L),include.mean=TRUE,
       method=method,optim.control=list(maxit=maxit)))
    fit<-ans$value;code<-if(is.null(fit))NA_integer_ else fit$code
    pa<-if(is.null(fit))numeric() else stats::coef(fit)
    ok<-is.null(ans$error) && isTRUE(code==0L) &&
        all(c("ar1","ma1","intercept")%in%names(pa)) && all(is.finite(pa)) &&
        abs(pa["ar1"])<1 && abs(pa["ma1"])<1
    rows[[method]]<-data.frame(stage=paste0("ARMA_",method),solver="BFGS",
      initialization="training_only",convergence=code,training_loglik=if(is.null(fit))NA_real_ else fit$loglik,
      usable=ok,chosen=FALSE,warning_count=nrow(ans$warnings),error=wrr_error_text(ans$error))
    if(nrow(ans$warnings))warnings[[method]]<-cbind(stage=paste0("ARMA_",method),ans$warnings)
    if(ok) {start<-list(mu=unname(pa["intercept"])*sc,ar1=unname(pa["ar1"]),ma1=unname(pa["ma1"]));break}
  }
  list(start=start,attempts=do.call(rbind,rows),warnings=wrr_bind(warnings,
       data.frame(stage=character(),message=character(),call=character())))
}
wrr_garch_margin <- function(x, train, model, arma_maxit=2000L, retry=TRUE) {
  if(any(!is.finite(x))||sum(train)<100L)stop("Invalid or insufficient marginal data")
  attempts<-warnings<-list();fits<-specs<-list()
  one<-function(label,start=list(),solver="hybrid") {
    spec<-rugarch::ugarchspec(variance.model=list(model=model,garchOrder=c(1,1)),
      mean.model=list(armaOrder=c(1,1),include.mean=TRUE),distribution.model="sstd",start.pars=start)
    ctrl<-if(solver=="nlminb")list(iter.max=2000L,eval.max=4000L) else list()
    ans<-wrr_capture(function()rugarch::ugarchfit(spec,x[train],solver=solver,
      solver.control=ctrl,fit.control=list(scale=1)))
    fit<-ans$value;code<-NA_integer_;ll<-NA_real_;ok<-FALSE
    if(!is.null(fit)) {
      chk<-wrr_capture(function() {
        code<-fit@fit$convergence;pa<-rugarch::coef(fit);sig<-as.numeric(rugarch::sigma(fit))
        ll<-as.numeric(rugarch::likelihood(fit))
        list(code=code,ll=ll,ok=isTRUE(code==0L)&&all(is.finite(pa))&&
          all(is.finite(sig)&sig>0)&&length(ll)==1L&&is.finite(ll))
      })
      if(is.null(chk$error)) {code<-chk$value$code;ll<-chk$value$ll;ok<-chk$value$ok}
      else ans$error<-chk$error
      ans$warnings<-rbind(ans$warnings,chk$warnings)
    }
    attempts[[label]]<<-data.frame(stage=label,solver=solver,
      initialization=if(length(start))"converged_training_ARMA" else "package_default",
      convergence=code,training_loglik=ll,usable=ok,chosen=FALSE,
      warning_count=nrow(ans$warnings),error=wrr_error_text(ans$error))
    if(nrow(ans$warnings))warnings[[label]]<<-cbind(stage=label,ans$warnings)
    if(ok) {fits[[label]]<<-fit;specs[[label]]<<-spec}
    list(ok=ok,arma_warning=any(grepl("arima",ans$warnings$call,ignore.case=TRUE)&
       grepl("convergence|optim",ans$warnings$message,ignore.case=TRUE)))
  }
  first<-one("GARCH_default")
  if(isTRUE(retry) && (!first$ok || first$arma_warning)) {
    pre<-wrr_arma_start(x[train],arma_maxit)
    for(i in seq_len(nrow(pre$attempts)))attempts[[pre$attempts$stage[i]]]<-pre$attempts[i,,drop=FALSE]
    if(nrow(pre$warnings))warnings[["ARMA_refit"]]<-pre$warnings
    if(!is.null(pre$start))one("GARCH_reinitialized",pre$start)
    if(!length(fits))one("GARCH_nlminb",if(is.null(pre$start))list() else pre$start,"nlminb")
  }
  tab<-do.call(rbind,attempts);warn<-wrr_bind(warnings,
    data.frame(stage=character(),message=character(),call=character()))
  fail<-function(msg)stop(structure(list(message=msg,call=NULL,attempts=tab,warning_log=warn),
    class=c("wrr_marginal_error","error","condition")))
  if(!length(fits))fail(paste("All numerical attempts failed for",model))
  # Select numerical attempts by TRAINING likelihood; validation chooses the family.
  ll<-vapply(fits,function(f)as.numeric(rugarch::likelihood(f)),numeric(1))
  chosen<-names(fits)[which.max(ll)];tab$chosen<-tab$stage==chosen
  fit<-fits[[chosen]];fixed<-specs[[chosen]];rugarch::setfixed(fixed)<-as.list(rugarch::coef(fit))
  filtered<-wrr_capture(function() {
    filt<-rugarch::ugarchfilter(fixed,data=x,n.old=sum(train),rec.init="all")
    mu<-as.numeric(rugarch::fitted(filt));sigma<-as.numeric(rugarch::sigma(filt))
    z<-as.numeric(rugarch::residuals(filt,standardize=TRUE));pa<-rugarch::coef(fit)
    u<-rugarch::pdist("sstd",z,mu=0,sigma=1,skew=pa["skew"],shape=pa["shape"])
    if(any(!is.finite(c(mu,sigma,z,u)))||any(sigma<=0)||any(u<0|u>1))
      stop("Non-finite or invalid conditional marginal output")
    list(u=wrr_clip(u),z=z,mu=mu,sigma=sigma,skew=unname(pa["skew"]),
      shape=unname(pa["shape"]),method=paste0(model,"_sstd"),fit=fit,fixed=fixed)
  })
  if(nrow(filtered$warnings))warn<-rbind(warn,cbind(stage="filter",filtered$warnings))
  if(!is.null(filtered$error))fail(conditionMessage(filtered$error))
  out<-filtered$value;out$attempts<-tab;out$warning_log<-warn;out$chosen_attempt<-chosen
  out
}
wrr_write_marginal_reports <- function(pits,cfg) {
  wrr_csv(pits$selection,file.path(cfg$out_dir,"tables/marginal_selection.csv"))
  if(!is.null(pits$attempts))wrr_csv(pits$attempts,file.path(cfg$out_dir,"tables/marginal_fit_attempts.csv"))
  if(!is.null(pits$warning_log))wrr_csv(pits$warning_log,file.path(cfg$out_dir,"tables/marginal_warnings.csv"))
}

wrr_margin_nll <- function(m, x, valid) {
  # Common proper predictive density criterion for the three GARCH candidates.
  density <- rugarch::ddist("sstd", x[valid], mu=m$mu[valid], sigma=m$sigma[valid],
                    skew=m$skew, shape=m$shape)
  if (any(!is.finite(density)) || any(density <= 0)) return(Inf)
  mean(-log(density))
}
wrr_margins <- function(dat,cfg) {
  selected<-choices<-attempts<-warnings<-list()
  empty_attempts<-data.frame(asset=character(),model=character(),stage=character(),solver=character(),
    initialization=character(),convergence=integer(),training_loglik=numeric(),usable=logical(),
    chosen=logical(),warning_count=integer(),error=character())
  empty_warnings<-data.frame(asset=character(),model=character(),stage=character(),message=character(),call=character())
  save_partial<-function()wrr_write_marginal_reports(list(selection=do.call(rbind,choices),
    attempts=wrr_bind(attempts,empty_attempts),warning_log=wrr_bind(warnings,empty_warnings)),cfg)
  for(a in colnames(dat$R)) {
    wrr_note("Marginal: ",a);x<-dat$R[,a]
    if(cfg$margin=="empirical") {
      selected[[a]]<-wrr_empirical_margin(x,cfg$margin_window)
      choices[[a]]<-data.frame(asset=a,model="rolling_empirical",validation_nll=NA_real_,selected=TRUE,
        chosen_attempt="none",final_convergence=NA_integer_,warning_count=0L,error="")
    } else {
      fit<-rows<-list()
      for(m in cfg$margin_candidates) {
        wrr_note("  Fit ",a," / ",m)
        ans<-wrr_capture(function()wrr_with_seed(wrr_seed(cfg$seed,a,m),
          wrr_garch_margin(x,dat$split=="train",m,cfg$arma_maxit,cfg$garch_retry)))
        obj<-if(is.null(ans$error))ans$value else ans$error
        at<-obj$attempts;wa<-obj$warning_log
        key<-paste(a,m,sep="/")
        if(!is.null(at)&&nrow(at))attempts[[key]]<-cbind(asset=a,model=m,at)
        if(!is.null(wa)&&nrow(wa))warnings[[key]]<-cbind(asset=a,model=m,wa)
        if(nrow(ans$warnings))warnings[[paste0(key,"/outer")]]<-cbind(asset=a,model=m,stage="outer",ans$warnings)
        good<-is.null(ans$error);score<-Inf;err<-wrr_error_text(ans$error)
        if(good) {
          sc<-wrr_capture(function()wrr_margin_nll(ans$value,x,dat$split=="validation"))
          if(nrow(sc$warnings))warnings[[paste0(key,"/score")]]<-cbind(asset=a,model=m,stage="density",sc$warnings)
          if(is.null(sc$error))score<-sc$value else err<-conditionMessage(sc$error)
          good<-is.finite(score)
          if(!good&&err=="")err<-"Non-finite validation density loss"
        }
        fit[m]<-list(if(good)ans$value else NULL)
        nw<-if(length(warnings))sum(vapply(warnings[startsWith(names(warnings),key)],nrow,integer(1))) else 0L
        rows[[m]]<-data.frame(asset=a,model=m,validation_nll=score,
          chosen_attempt=if(good)ans$value$chosen_attempt else "none",
          final_convergence=if(good)0L else NA_integer_,warning_count=nw,error=err)
        wrr_note("  ",a," / ",m,": ",if(good)"accepted" else "rejected",
          "; warnings=",nw,if(good)paste0("; attempt=",ans$value$chosen_attempt) else paste0("; ",err))
      }
      tab<-do.call(rbind,rows);best<-which.min(tab$validation_nll)
      tab$selected<-FALSE
      if(is.finite(tab$validation_nll[best]))tab$selected[best]<-TRUE
      choices[[a]]<-tab;save_partial()
      if(!is.finite(tab$validation_nll[best]))stop("All GARCH candidates failed for ",a,"; inspect marginal_fit_attempts.csv")
      selected[[a]]<-fit[[tab$model[best]]]
      wrr_note("Selected ",a,": ",tab$model[best])
    }
    save_partial()
  }
  U<-sapply(selected,`[[`,"u");Z<-sapply(selected,`[[`,"z")
  keep<-seq_len(nrow(U))>cfg$margin_window
  if(any(!is.finite(U[keep,,drop=FALSE]))||any(!is.finite(Z[keep,,drop=FALSE])))
    stop("Invalid marginal values after burn-in; refusing to silently delete evaluation dates")
  if(!any(keep&dat$split=="train"))stop("Insufficient training history after burn-in")
  list(U=U[keep,,drop=FALSE],Z=Z[keep,,drop=FALSE],R=dat$R[keep,,drop=FALSE],
    raw=dat$raw[keep,,drop=FALSE],v=dat$v[keep],dates=dat$dates[keep],split=dat$split[keep],
    original_row=which(keep),models=selected,selection=do.call(rbind,choices),bounds=dat$bounds,
    attempts=wrr_bind(attempts,empty_attempts),warning_log=wrr_bind(warnings,empty_warnings))
}

wrr_margin_inverse <- function(m,p,original_t) {
  if (m$method == "rolling_empirical")
    return(as.numeric(quantile(m$x[(original_t-m$window):(original_t-1)],p,type=8,names=FALSE)))
  rugarch::qdist("sstd",p,mu=m$mu[original_t],sigma=m$sigma[original_t],skew=m$skew,shape=m$shape)
}
wrr_regimes <- function(v,train,probs) {
  b <- as.numeric(quantile(v[train],probs,type=8))
  list(state=ifelse(v<=b[1],"calm",ifelse(v<=b[2],"medium","stress")),cuts=b)
}
wrr_event <- function(U,cols,q=.1,k=3L) {
  U <- as.matrix(U)
  if (any(!is.finite(U[,cols,drop=FALSE]))) stop("Missing scenario/PIT values")
  as.integer(rowSums(U[,cols,drop=FALSE]<=q)>=k)
}
wrr_weighted_ranks <- function(U,p) {
  # Weighted mid-CDF with ties grouped. This changes the estimation target:
  # it is the copula of the tilted empirical law, not a new marginal PIT forecast.
  p <- p/sum(p); out <- matrix(NA_real_,nrow(U),ncol(U),dimnames=dimnames(U))
  for (j in seq_len(ncol(U))) {
    vals <- sort(unique(U[,j])); g <- match(U[,j],vals)
    masses <- as.numeric(rowsum(p,g,reorder=TRUE)); mids <- cumsum(masses)-masses/2
    out[,j] <- mids[g]
  }
  wrr_clip(out)
}
wrr_weighted_ks <- function(u,p) {
  vals <- sort(unique(u)); g <- match(u,vals); mass <- as.numeric(rowsum(p/sum(p),g,reorder=TRUE))
  max(abs(cumsum(mass)-vals),abs(cumsum(mass)-mass-vals))
}
wrr_weights <- function(U,v,states,past,t,cfg,variant="WRR_direct") {
  n <- length(past); age <- (max(past)-past)/max(1,n)
  rec <- if(variant=="WRR_no_recency") 0 else cfg$recency_decay
  sb <- if(variant=="WRR_no_regime") 1 else cfg$same_regime_bonus
  vd <- if(variant=="WRR_no_vix_distance") 0 else cfg$vix_distance_decay
  tb <- if(variant=="WRR_no_tail") 0 else cfg$tail_bonus
  sc <- sd(v[past]); if(!is.finite(sc)||sc<1e-6) sc <- 1
  hit <- wrr_event(U[past,,drop=FALSE],intersect(cfg$common_assets,colnames(U)),cfg$q,cfg$k)
  lw <- -rec*age + log(ifelse(states[past]==states[t],sb,1)) - vd*abs(v[past]-v[t])/sc + log1p(tb*hit)
  p <- exp(lw-max(lw)); p <- p/sum(p); raw_ess <- 1/sum(p^2)
  # Fixed ESS floor blends toward uniform weights; it never uses current returns.
  mix <- 0
  if(raw_ess < cfg$ess_min_fraction*n) {
    mix <- uniroot(function(a)1/sum(((1-a)*p+a/n)^2)-cfg$ess_min_fraction*n,c(0,1))$root
    p <- (1-mix)*p+mix/n
  }
  list(p=p,raw_ess=raw_ess,ess=1/sum(p^2),uniform_mix=mix)
}
wrr_vine_fit <- function(U,weights=NULL,family="tail",cfg) {
  families <- if(family=="tail")cfg$family_tail else cfg$family_general
  if (is.null(weights)) weights <- rep(1,nrow(U)) else weights <- weights/sum(weights)*nrow(U)
  VineCopula::RVineStructureSelect(U,familyset=families,weights=weights,
       selectioncrit="AIC",indeptest=FALSE,rotations=TRUE,trunclevel=cfg$trunclevel,
       treecrit="tau",progress=FALSE,cores=1)
}
wrr_gaussian_ll <- function(U,R) {
  Z <- qnorm(wrr_clip(as.matrix(U))); R <- wrr_pd(R); d <- ncol(Z)
  -.5*as.numeric(determinant(R,logarithm=TRUE)$modulus)-.5*rowSums((Z%*%(solve(R)-diag(d)))*Z)
}
wrr_model_ll <- function(model,U) {
  U <- as.matrix(U)
  if(model$type=="gaussian") return(wrr_gaussian_ll(U,model$R))
  if(model$type=="t") return(as.numeric(copula::dCopula(U,model$cop,log=TRUE)))
  if(model$type=="vine") {
    d <- as.numeric(VineCopula::RVinePDF(U,model$fit))
    if(length(d)!=nrow(U)||any(!is.finite(d))||any(d<=0)) stop("Invalid vine density")
    return(log(d))
  }
  if(model$type=="mixture") {
    ll <- sapply(model$members,wrr_model_ll,U=U); ll <- matrix(ll,nrow=nrow(U))
    m <- apply(ll,1,max); return(m+log(rowMeans(exp(ll-m))))
  }
  stop("Unknown model type")
}
wrr_sim <- function(model,n,assets,seed) wrr_with_seed(seed, {
  d <- length(assets)
  if(model$type=="gaussian") U <- pnorm(matrix(rnorm(n*d),n,d)%*%chol(wrr_pd(model$R)))
  else if(model$type=="t") U <- copula::rCopula(n,model$cop)
  else if(model$type=="vine") U <- VineCopula::RVineSim(n,model$fit)
  else if(model$type=="mixture") {
    group <- sample.int(length(model$members),n,replace=TRUE); U <- matrix(NA_real_,n,d)
    for(g in unique(group)) U[group==g,] <- wrr_sim(model$members[[g]],sum(group==g),assets,wrr_seed(seed,g))
  } else stop("Unknown simulation model")
  U <- as.matrix(U);colnames(U) <- assets
  if(any(!is.finite(U)))stop("Invalid simulated scenarios")
  wrr_clip(U)
})
wrr_fit_window <- function(U,v,states,past,t,cfg,variant,params,seed) {
  uniform <- variant %in% c("rolling_vine","rolling_matched")
  w <- if(uniform)list(p=rep(1/length(past),length(past)),raw_ess=length(past),ess=length(past),uniform_mix=0) else
    wrr_weights(U,v,states,past,t,cfg,variant)
  X <- U[past,,drop=FALSE]
  raw <- variant %in% c("WRR_raw","WRR_legacy_resampled_B1","rolling_vine")
  V <- if(raw)X else wrr_weighted_ranks(X,w$p)
  if(grepl("resampled_B",variant)) {
    B <- as.integer(sub(".*_B","",variant)); members <- vector("list",B); uniques <- integer(B)
    for(b in seq_len(B)) {
      # Nested seeds: B=1 and B=5 are prefixes of B=20 at the same refit.
      id <- wrr_with_seed(wrr_seed(seed,"resample",b),sample.int(length(past),cfg$resample_n,TRUE,w$p))
      members[[b]] <- list(type="vine",fit=wrr_vine_fit(V[id,,drop=FALSE],NULL,params$family,cfg))
      uniques[b] <- length(unique(id))
    }
    model <- list(type="mixture",members=members)
  } else { B <- 0L; uniques <- length(past)
    model <- list(type="vine",fit=wrr_vine_fit(V,w$p,params$family,cfg))
  }
  diag <- data.frame(t=t,n_window=length(past),sample_size=if(B>0)cfg$resample_n else length(past),
       resampling_replicates=B,unique_mean=mean(uniques),raw_ess=w$raw_ess,ess=w$ess,
       uniform_mix=w$uniform_mix,max_weight=max(w$p),
       weighted_ks_before=mean(vapply(seq_len(ncol(X)),function(j)wrr_weighted_ks(X[,j],w$p),numeric(1))),
       weighted_ks_after=mean(vapply(seq_len(ncol(V)),function(j)wrr_weighted_ks(V[,j],w$p),numeric(1))))
  list(model=model,diagnostics=diag)
}
wrr_roll <- function(U,v,states,positions,cfg,variant,params,seed,cache_dir=NULL) {
  fits <- list(); map <- integer(length(positions)); diagnostics <- list(); last <- -Inf; state <- NA_character_
  for(i in seq_along(positions)) {
    t <- positions[i]
    changed <- cfg$refit_on_regime_change && !identical(state,states[t])
    if(i==1 || t-last >= params$refit || changed) {
      if(length(fits)%%cfg$progress_refits==0L)wrr_note(variant," refit ",length(fits)+1L,"; forecast index ",t)
      past <- seq.int(max(1L,t-params$window),t-1L)
      if(length(past)<params$window)stop("Insufficient rolling history")
      calc <- function() wrr_fit_window(U,v,states,past,t,cfg,variant,params,wrr_seed(seed,t))
      key <- wrr_hash(list(U[past,,drop=FALSE],v[c(past,t)],states[c(past,t)],cfg,variant,params,seed,t))
      ans <- if(is.null(cache_dir))calc() else wrr_cache(file.path(cache_dir,paste0(t,".rds")),key,calc,cfg$resume)
      fits[[length(fits)+1L]] <- ans$model; diagnostics[[length(fits)]] <- ans$diagnostics
      last <- t; state <- states[t]
    }
    map[i] <- length(fits)
  }
  list(fits=fits,map=map,positions=positions,diagnostics=do.call(rbind,diagnostics))
}
wrr_roll_ll <- function(roll,U) {
  ll <- numeric(length(roll$positions))
  for(j in seq_along(roll$fits)) {
    id <- which(roll$map==j);ll[id] <- wrr_model_ll(roll$fits[[j]],U[roll$positions[id],,drop=FALSE])
  }
  if(any(!is.finite(ll)))stop("Non-finite daily likelihood")
  ll
}
wrr_tune <- function(U,v,states,pos,cfg,variant,seed,root) {
  rows <- list()
  for(g in seq_len(nrow(cfg$grid))) {
    p <- as.list(cfg$grid[g,]); wrr_note("Tune ",variant," ",g,"/",nrow(cfg$grid))
    roll <- wrr_roll(U,v,states,pos,cfg,variant,p,seed,file.path(root,paste0("grid",g)))
    ll <- wrr_roll_ll(roll,U)
    rows[[g]] <- cbind(cfg$grid[g,,drop=FALSE],validation_nll=mean(-ll),n_validation=length(ll))
  }
  tab <- do.call(rbind,rows); best <- order(tab$validation_nll,tab$window,tab$refit,tab$family)[1]
  list(params=as.list(tab[best,c("window","refit","family")]),table=tab)
}
wrr_dcc <- function(U,train,valid,cfg) {
  # Gaussian scores are qnorm(PIT), never raw skew-t GARCH residuals.
  Z <- qnorm(wrr_clip(U)); S <- wrr_pd(cor(Z[train,,drop=FALSE]))
  filter <- function(a,b,upto) {
    Q<-S; rr<-vector("list",upto)
    for(t in seq_len(upto)) {
      if(t>1)Q<-(1-a-b)*S+a*tcrossprod(Z[t-1,])+b*Q
      rr[[t]]<-wrr_pd(cov2cor(Q))
    };rr
  }
  loss<-numeric(nrow(cfg$dcc_grid));lastv<-max(valid)
  for(g in seq_len(nrow(cfg$dcc_grid))) {
    rr<-filter(cfg$dcc_grid$a[g],cfg$dcc_grid$b[g],lastv)
    loss[g]<-mean(vapply(valid,function(t)-wrr_gaussian_ll(U[t,,drop=FALSE],rr[[t]]),numeric(1)))
  }
  b<-which.min(loss);p<-cfg$dcc_grid[b,]
  list(R=filter(p$a,p$b,nrow(U)),tuning=cbind(cfg$dcc_grid,validation_nll=loss))
}

wrr_energy <- function(X,y) {
  # Independent disjoint pairs give an unbiased MC estimate of the distance term.
  M<-nrow(X); m<-floor(M/2); if(m<1)stop("At least two scenarios required")
  mean(sqrt(rowSums(sweep(X,2,y,"-")^2))) -
    .5*mean(sqrt(rowSums((X[seq_len(m),,drop=FALSE]-X[m+seq_len(m),,drop=FALSE])^2)))
}
wrr_variogram <- function(X,y) {
  pairs<-combn(seq_along(y),2)
  sum(apply(pairs,2,function(ij)(abs(y[ij[1]]-y[ij[2]])^.5-
                   mean(abs(X[,ij[1]]-X[,ij[2]])^.5))^2))
}
wrr_scenario_metrics <- function(S,u,cfg) {
  Z<-qnorm(wrr_clip(S));z<-qnorm(wrr_clip(u)); common<-intersect(cfg$common_assets,colnames(S))
  main_hit<-wrr_event(matrix(u,1,dimnames=list(NULL,names(u))),common,cfg$q,cfg$k)
  h<-wrr_event(S,common,cfg$q,cfg$k);prob<-mean(h)
  base<-data.frame(energy=wrr_energy(Z,z),variogram=wrr_variogram(Z,z),
    energy_common=wrr_energy(Z[,common,drop=FALSE],z[common]),
    variogram_common=wrr_variogram(Z[,common,drop=FALSE],z[common]),
    probability=prob,observed=main_hit,brier=(prob-main_hit)^2,
    brier_mc_adjusted=(prob-main_hit)^2-prob*(1-prob)/(nrow(S)-1),
    probability_mc_se=sqrt(prob*(1-prob)/nrow(S)))
  grid<-expand.grid(q=c(.05,.10),k=c(2L,3L,4L))
  grid$probability<-grid$observed<-NA_real_
  for(i in seq_len(nrow(grid))) {
    grid$probability[i]<-mean(wrr_event(S,common,grid$q[i],grid$k[i]))
    grid$observed[i]<-wrr_event(matrix(u,1,dimnames=list(NULL,names(u))),common,grid$q[i],grid$k[i])
  }
  grid$brier<-(grid$probability-grid$observed)^2
  if("VIX"%in%colnames(S)) {
    # VIX rises in typical equity stress; never count its lower tail as ETF distress.
    base$joint6_vix_upper_probability<-mean(h==1L & S[,"VIX"]>=1-cfg$q)
    base$joint6_vix_upper_observed<-as.integer(main_hit==1L & u["VIX"]>=1-cfg$q)
  }
  list(base=base,grid=grid)
}
wrr_returns <- function(S,pits,t) {
  out<-S
  for(a in colnames(S))out[,a]<-wrr_margin_inverse(pits$models[[a]],S[,a],pits$original_row[t])
  out
}
wrr_portfolio <- function(S,pits,t,cfg) {
  # Equal-weight, daily rebalanced six-ETF diagnostic. VIX is not a tradable holding.
  x<-wrr_returns(S[,cfg$common_assets,drop=FALSE],pits,t)
  sim<-rowMeans(expm1(x));obs<-mean(expm1(pits$raw[t,cfg$common_assets]))
  q<-unname(quantile(sim,cfg$risk_alpha,type=8));e<-mean(sim[sim<=q]);hit<-as.integer(obs<=q)
  data.frame(portfolio_return=obs,VaR=q,ES=e,violation=hit,
    quantile_loss=(cfg$risk_alpha-as.integer(obs<q))*(obs-q),
    es_identification=e-q+(q-obs)*hit/cfg$risk_alpha,
    n_es_draws=sum(sim<=q))
}
wrr_covar <- function(S,pits,t,cfg,condition="XLE",target="XLF") {
  # Condition on the marginal PIT threshold, not a newly estimated sample quantile.
  y<-wrr_margin_inverse(pits$models[[target]],S[,target],pits$original_row[t])
  idx<-S[,condition]<=cfg$risk_alpha; yc<-y[idx]
  q<-if(length(yc)>10)unname(quantile(yc,cfg$risk_alpha,type=8)) else NA_real_
  n_es<-if(is.finite(q))sum(yc<=q) else 0L
  reliable<-n_es>=cfg$covar_min_joint_tail
  data.frame(condition=condition,target=target,CoVaR=if(reliable)q else NA_real_,
       CoES=if(reliable)mean(yc[yc<=q]) else NA_real_,n_condition=sum(idx),
       n_joint_tail=n_es,reliable=reliable)
}
wrr_hac <- function(d,lag=10L) {
  n<-length(d);if(n<5||any(!is.finite(d)))return(c(se=NA_real_,p=NA_real_))
  z<-d-mean(d);v<-sum(z*z)/n; lag<-min(lag,n-2L)
  for(k in seq_len(lag))v<-v+2*(1-k/(lag+1))*sum(z[(k+1):n]*z[1:(n-k)])/n
  se<-sqrt(max(0,v)/n)
  c(se=se,p=if(se>0)2*pnorm(-abs(mean(d)/se)) else if(mean(d)==0)1 else 0)
}
wrr_block_indices <- function(n,block) {
  starts<-sample.int(n,ceiling(n/block),TRUE)
  as.integer(unlist(lapply(starts,function(s)((s-1+seq_len(block)-1)%%n)+1)))[seq_len(n)]
}
wrr_inference <- function(daily,cfg) {
  # Align dates explicitly; no model-dependent na.rm deletion.
  ref<-daily[daily$model=="WRR_direct",];rows<-list()
  metrics<-c("nll","energy_common","variogram_common","brier","quantile_loss")
  for(m in setdiff(unique(daily$model),"WRR_direct")) {
    b<-daily[daily$model==m,];i<-match(ref$date,b$date)
    if(anyNA(i)||length(i)!=nrow(b))stop("Models have different evaluation dates")
    for(metric in metrics) {
      d<-ref[[metric]]-b[[metric]][i]; if(any(!is.finite(d)))stop("Non-finite paired losses")
      hac<-wrr_hac(d,cfg$hac_lag)
      for(block in cfg$block_lengths) {
        block<-min(block,length(d))
        bs<-wrr_with_seed(wrr_seed(cfg$seed,m,metric,block),
              replicate(cfg$inference_B,mean(d[wrr_block_indices(length(d),block)])))
        ci<-quantile(bs,c(.025,.975),names=FALSE)
        rows[[length(rows)+1L]]<-data.frame(benchmark=m,metric=metric,n=length(d),
          difference=mean(d),hac_se=hac["se"],p_two_sided=hac["p"],
          block=block,bootstrap_B=cfg$inference_B,ci_low=ci[1],ci_high=ci[2])
      }
    }
  }
  tab<-do.call(rbind,rows);tab$p_holm<-NA_real_
  # One declared family: all reference-vs-benchmark tests and five metrics per universe.
  for(block in unique(tab$block)) {
    i<-which(tab$block==block);tab$p_holm[i]<-p.adjust(tab$p_two_sided[i],"holm")
  };tab
}
wrr_pit_diagnostics <- function(pits,cfg) {
  rows<-list();qs<-c(.01,.05,.10,.50,.90,.95,.99)
  for(a in colnames(pits$U)) for(s in c("validation","test")) {
    id<-pits$split==s;u<-pits$U[id,a];z<-pits$Z[id,a];n<-length(u)
    if(n<10)next
    ks<-suppressWarnings(ks.test(u,"punif",exact=FALSE))
    for(q in qs) {
      hit<-as.numeric(u<=q);d<-hit-q;hac<-wrr_hac(d,cfg$hac_lag)
      bs<-wrr_with_seed(wrr_seed(cfg$seed,a,s,q),replicate(cfg$inference_B,
           mean(hit[wrr_block_indices(n,min(20L,n))])))
      ci<-quantile(bs,c(.025,.975),names=FALSE)
      rows[[length(rows)+1L]]<-data.frame(asset=a,split=s,n=n,q=q,hit_rate=mean(hit),
        hit_ci_low=ci[1],hit_ci_high=ci[2],calibration_hac_p=hac["p"],
        ks_stat=unname(ks$statistic),ks_p_iid_reference=ks$p.value,
        lb_z_p=Box.test(z,lag=min(10,n-2),type="Ljung-Box")$p.value,
        lb_z2_p=Box.test(z^2,lag=min(10,n-2),type="Ljung-Box")$p.value)
    }
  }
  tab<-do.call(rbind,rows)
  tab$calibration_holm_p<-ave(tab$calibration_hac_p,tab$split,FUN=function(x)p.adjust(x,"holm"))
  tab
}
wrr_network <- function(U,states,train,q=.1) {
  rows<-list();pairs<-combn(colnames(U),2,simplify=FALSE)
  for(s in c("calm","medium","stress")) for(pa in pairs) {
    X<-U[train & states==s,pa,drop=FALSE];n<-nrow(X);if(n<10)next
    ni<-sum(X[,1]<=q);nj<-sum(X[,2]<=q);nij<-sum(X[,1]<=q & X[,2]<=q)
    rows[[length(rows)+1L]]<-data.frame(regime=s,asset_i=pa[1],asset_j=pa[2],n=n,
      kendall_tau=cor(X[,1],X[,2],method="kendall"),n_i=ni,n_j=nj,n_both=nij,
      finite_q_cond_j_given_i=if(ni>0)nij/ni else NA_real_,
      finite_q_cond_i_given_j=if(nj>0)nij/nj else NA_real_,
      cohit_rate=nij/n)
  };do.call(rbind,rows)
}
wrr_summary <- function(daily) {
  met<-c("nll","energy","variogram","energy_common","variogram_common","brier",
         "brier_mc_adjusted","probability","observed","quantile_loss","violation","es_identification")
  out<-aggregate(daily[met],list(model=daily$model),mean)
  out$n<-as.integer(table(daily$model)[out$model]);out$mean_log_score<--out$nll
  out[order(out$nll),]
}
wrr_fit_base <- function(U,train,states,cfg) {
  X<-U[train,,drop=FALSE]
  ga<-list(type="gaussian",R=wrr_pd(cor(qnorm(X))))
  tc<-copula::tCopula(dim=ncol(X),dispstr="un",df=6,df.fixed=FALSE)
  fit<-copula::fitCopula(tc,X,method="ml",estimate.variance=FALSE)
  if(!is.null(fit@fitting.stats$convergence) && fit@fitting.stats$convergence != 0)
    stop("Student-t copula optimizer did not converge")
  tv<-list(type="t",cop=fit@copula)
  sv<-list(type="vine",fit=wrr_vine_fit(X,NULL,"general",cfg))
  regimes<-list();info<-list()
  for(s in c("calm","medium","stress")) {
    idx<-train & states==s;n<-sum(idx)
    regimes[[s]]<-if(n>=cfg$min_regime_n)list(type="vine",fit=wrr_vine_fit(U[idx,,drop=FALSE],NULL,"general",cfg)) else sv
    info[[s]]<-data.frame(state=s,n=n,static_fallback=n<cfg$min_regime_n)
  }
  list(gaussian=ga,t_copula=tv,static_vine=sv,regimes=regimes,regime_info=do.call(rbind,info))
}
wrr_get_model <- function(models,name,i,t,state) {
  if(name=="dcc_gaussian")return(list(type="gaussian",R=models$dcc$R[[t]]))
  if(name=="regime_vine")return(models$base$regimes[[state]])
  if(name%in%c("gaussian","t_copula","static_vine"))return(models$base[[name]])
  r<-models$rolling[[name]]; r$fits[[r$map[i]]]
}
wrr_origin_audits <- function(U,pits,states,test,params,cfg,root,seed) {
  # Origins fixed by chronological spacing, not ex-post event magnitude.
  ids<-unique(round(seq(1,length(test),length.out=min(cfg$audit_origins,length(test)))))
  mc<-uncertainty<-covars<-list()
  for(i in ids) {
    t<-test[i];wrr_note("Audit origin ",as.character(pits$dates[t]))
    past<-seq.int(t-params$window,t-1L)
    primary<-wrr_fit_window(U,pits$v,states,past,t,cfg,"WRR_direct",params,wrr_seed(seed,t))$model
    u<-U[t,];names(u)<-colnames(U)
    for(M in cfg$scenario_sensitivity)for(rep in seq_len(cfg$mc_repeats)) {
      S<-wrr_sim(primary,M,colnames(U),wrr_seed(seed,t,"MC",M,rep))
      z<-wrr_scenario_metrics(S,u,cfg)$base
      z$censored_nll<-if(z$observed==1)-wrr_model_ll(primary,U[t,,drop=FALSE]) else -log1p(-pmin(1-1e-10,z$probability))
      mc[[length(mc)+1L]]<-cbind(data.frame(date=pits$dates[t],M=M,rep=rep),z)
    }
    for(rep in seq_len(cfg$covar_repeats)) {
      S<-wrr_sim(primary,cfg$covar_scenarios,colnames(U),wrr_seed(seed,t,"covar",rep))
      for(pa in list(c("XLE","XLF"),c("ICLN","XLF"),c("TAN","XLE"),c("GLD","XLF"))) {
        z<-wrr_covar(S,pits,t,cfg,pa[1],pa[2])
        covars[[length(covars)+1L]]<-cbind(data.frame(date=pits$dates[t],regime=states[t],rep=rep),z)
      }
    }
    # Circular block bootstrap of local (PIT,VIX,state) tuples. Refit structure and
    # parameters; margins, threshold cutoffs and chosen hyperparameters stay fixed.
    # These intervals are conditional local estimation uncertainty, not full-pipeline CIs.
    for(b in seq_len(cfg$estimation_B)) {
      if(b==1L || b%%25L==0L || b==cfg$estimation_B)wrr_note("Local estimation bootstrap ",b,"/",cfg$estimation_B)
      z<-wrr_with_seed(wrr_seed(seed,t,"estimation",b), {
        idx<-past[wrr_block_indices(length(past),min(cfg$estimation_block,length(past)))]
        Us<-rbind(U[idx,,drop=FALSE],U[t,,drop=FALSE]);vs<-c(pits$v[idx],pits$v[t]);ss<-c(states[idx],states[t])
        ans<-tryCatch({
          model<-wrr_fit_window(Us,vs,ss,seq_along(past),length(past)+1,cfg,"WRR_direct",params,wrr_seed(seed,t,b))$model
          S<-wrr_sim(model,cfg$estimation_scenarios,colnames(U),wrr_seed(seed,t,"est_scen",b))
          data.frame(probability=mean(wrr_event(S,cfg$common_assets,cfg$q,cfg$k)),
                     log_score=wrr_model_ll(model,U[t,,drop=FALSE]),error="")
        },error=function(e)data.frame(probability=NA_real_,log_score=NA_real_,error=conditionMessage(e)))
        ans
      })
      uncertainty[[length(uncertainty)+1L]]<-cbind(data.frame(date=pits$dates[t],b=b),z)
    }
    # Save after every audit origin to preserve completed diagnostics.
    wrr_csv(do.call(rbind,mc),file.path(root,"tables/scenario_mc_replicates.csv"))
    wrr_csv(do.call(rbind,covars),file.path(root,"tables/covar_mc_replicates.csv"))
    wrr_csv(do.call(rbind,uncertainty),file.path(root,"tables/local_estimation_bootstrap_draws.csv"))
  }
  un<-do.call(rbind,uncertainty);summ<-do.call(rbind,lapply(split(un,un$date),function(x) {
    good<-is.finite(x$probability);ok<-mean(good)>=.95
    ci<-if(ok)quantile(x$probability[good],c(.025,.975),names=FALSE) else c(NA,NA)
    data.frame(date=x$date[1],B=nrow(x),success=sum(good),failure_rate=mean(!good),
               probability_low=ci[1],probability_high=ci[2],interval_reportable=ok)
  }))
  wrr_csv(summ,file.path(root,"tables/local_estimation_intervals.csv"))
  co<-do.call(rbind,covars)
  group<-interaction(co$date,co$condition,co$target,drop=TRUE)
  cos<-do.call(rbind,lapply(split(co,group),function(x)data.frame(date=x$date[1],regime=x$regime[1],
    condition=x$condition[1],target=x$target[1],CoVaR=if(all(x$reliable))mean(x$CoVaR) else NA_real_,
    CoES=if(all(x$reliable))mean(x$CoES) else NA_real_,
    CoVaR_mc_sd=if(all(x$reliable))sd(x$CoVaR) else NA_real_,
    CoES_mc_sd=if(all(x$reliable))sd(x$CoES) else NA_real_,
    min_joint_tail=min(x$n_joint_tail),repeats=nrow(x))))
  wrr_csv(cos,file.path(root,"tables/covar_mc_summary.csv"))
  invisible(summ)
}

wrr_universe <- function(pits,assets,label,cfg,light=FALSE) {
  root<-file.path(cfg$out_dir,label);dir.create(root,recursive=TRUE,showWarnings=FALSE)
  U<-pits$U[,assets,drop=FALSE];train<-pits$split=="train";valid<-which(pits$split=="validation")
  test<-which(pits$split=="test");if(is.finite(cfg$max_test))test<-head(test,cfg$max_test)
  if(max(cfg$grid$window)>=min(valid))stop("Tuning windows exceed available training PIT history")
  rg<-wrr_regimes(pits$v,train,cfg$regime_probs);states<-rg$state;seed<-wrr_seed(cfg$seed,label)
  wrr_csv(data.frame(date=pits$dates,split=pits$split,lagged_vix=pits$v,regime=states),file.path(root,"tables/regimes.csv"))
  wrr_csv(data.frame(prob=cfg$regime_probs,cut=rg$cuts),file.path(root,"tables/regime_cutoffs.csv"))
  key<-wrr_hash(list(U,pits$v,states,cfg))
  base<-wrr_cache(file.path(root,"cache/base.rds"),key,
    function()wrr_with_seed(seed,wrr_fit_base(U,train,states,cfg)),cfg$resume)
  dcc<-wrr_cache(file.path(root,"cache/dcc.rds"),key,
    function()wrr_dcc(U,train,valid,cfg),cfg$resume)
  wrr_csv(base$regime_info,file.path(root,"tables/regime_model_fallbacks.csv"))
  wrr_csv(dcc$tuning,file.path(root,"tables/dcc_validation.csv"))
  tun<-wrr_cache(file.path(root,"cache/tune_wrr.rds"),key,
    function()wrr_tune(U,pits$v,states,valid,cfg,"WRR_direct",seed,file.path(root,"cache/tune_wrr")),cfg$resume)
  tunroll<-wrr_cache(file.path(root,"cache/tune_rolling.rds"),key,
    function()wrr_tune(U,pits$v,states,valid,cfg,"rolling_vine",seed,file.path(root,"cache/tune_rolling")),cfg$resume)
  wrr_csv(tun$table,file.path(root,"tables/tuning_wrr.csv"));wrr_csv(tunroll$table,file.path(root,"tables/tuning_rolling.csv"))
  wrr_csv(as.data.frame(tun$params),file.path(root,"tables/selected_wrr.csv"))
  p<-tun$params;rolling<-list();refits<-list()
  variants<-c("WRR_direct","rolling_vine","rolling_matched")
  if(!light)variants<-c(variants,"WRR_raw","WRR_no_tail","WRR_no_recency","WRR_no_regime",
                          "WRR_no_vix_distance","WRR_other_family","WRR_legacy_resampled_B1")
  for(name in variants) {
    wrr_note(label," fit ",name)
    pp<-if(name=="rolling_vine")tunroll$params else p
    if(name=="WRR_other_family")pp$family<-if(p$family=="tail")"general" else "tail"
    rolling[[name]]<-wrr_roll(U,pits$v,states,test,cfg,name,pp,seed,file.path(root,"cache",name))
  }
  if(!light) {
    Bmax<-max(cfg$ensemble_B);name<-paste0("WRR_resampled_B",Bmax)
    big<-wrr_roll(U,pits$v,states,test,cfg,name,p,seed,file.path(root,"cache",name))
    for(B in cfg$ensemble_B) {
      small<-big;small$fits<-lapply(big$fits,function(x)list(type="mixture",members=x$members[seq_len(B)]))
      small$diagnostics$resampling_replicates<-B
      # Mean uniqueness inherited from Bmax is not a statistic for the prefix: omit it.
      if(B!=Bmax)small$diagnostics$unique_mean<-NA_real_
      rolling[[paste0("WRR_resampled_B",B)]]<-small
    }
    for(probs in cfg$regime_sensitivity) {
      rs<-wrr_regimes(pits$v,train,probs)$state;name<-paste0("WRR_threshold_",paste(round(100*probs),collapse="_"))
      rolling[[name]]<-wrr_roll(U,pits$v,rs,test,cfg,"WRR_direct",p,seed,file.path(root,"cache",name))
    }
  }
  for(name in names(rolling))refits[[name]]<-cbind(model=name,rolling[[name]]$diagnostics)
  refit<-do.call(rbind,refits);refit$date<-pits$dates[refit$t]
  wrr_csv(refit,file.path(root,"tables/refit_weight_diagnostics.csv"))
  models<-list(base=base,dcc=dcc,rolling=rolling)
  namesm<-c("gaussian","t_copula","static_vine","regime_vine","dcc_gaussian",names(rolling))
  daily<-grid<-list()
  for(i in seq_along(test)) {
    t<-test[i];if(i%%25==1)wrr_note(label," scenario day ",i,"/",length(test))
    # Per-day cache: interrupted runs continue without dropping old/new daily rows.
    day<-wrr_cache(file.path(root,"cache/daily",paste0(t,".rds")),
      wrr_hash(list(key,namesm,t,p,tunroll$params)),function() {
        d<-g<-list()
        for(name in namesm) {
          model<-wrr_get_model(models,name,i,t,states[t]);u<-U[t,];names(u)<-assets
          ll<-wrr_model_ll(model,U[t,,drop=FALSE])
          S<-wrr_sim(model,cfg$scenarios,assets,wrr_seed(seed,t,name,"scenario"))
          met<-wrr_scenario_metrics(S,u,cfg)
          risk<-wrr_portfolio(S,pits,t,cfg)
          # Censored log score on the six-asset event in the FULL copula space.
          # Its region probability is MC approximated; assess M sensitivity separately.
          pc<-pmin(1-1e-10,pmax(1e-10,met$base$probability))
          cens<-if(met$base$observed==1)-ll else -log1p(-pc)
          d[[name]]<-cbind(data.frame(date=pits$dates[t],regime=states[t],model=name,nll=-ll,
                              censored_nll=cens),met$base,risk)
          g[[name]]<-cbind(data.frame(date=pits$dates[t],regime=states[t],model=name),met$grid)
        }
        list(daily=do.call(rbind,d),grid=do.call(rbind,g))
      },cfg$resume)
    daily[[i]]<-day$daily;grid[[i]]<-day$grid
  }
  daily<-do.call(rbind,daily);grid<-do.call(rbind,grid)
  sumtab<-wrr_summary(daily);inf<-wrr_inference(daily,cfg)
  wrr_csv(daily,file.path(root,"tables/daily_forecasts.csv"));wrr_csv(grid,file.path(root,"tables/tail_grid_daily.csv"))
  wrr_csv(sumtab,file.path(root,"tables/model_comparison.csv"));wrr_csv(inf,file.path(root,"tables/paired_inference.csv"))
  wrr_csv(aggregate(daily[c("nll","energy_common","brier","observed")],
      list(model=daily$model,regime=daily$regime),mean),file.path(root,"tables/by_regime.csv"))
  wrr_csv(as.data.frame(table(daily[daily$model=="WRR_direct",c("regime","observed")])),
      file.path(root,"tables/regime_event_counts.csv"))
  wrr_csv(wrr_network(U,states,train,cfg$q),file.path(root,"tables/finite_threshold_network.csv"))
  wrr_figures(pits,states,daily,root,cfg)
  if(cfg$run_audits && !light)wrr_origin_audits(U,pits,states,test,p,cfg,root,seed)
  list(daily=daily,summary=sumtab,inference=inf,params=p,root=root)
}
wrr_figures <- function(pits,states,daily,root,cfg) {
  dir.create(file.path(root,"figures"),showWarnings=FALSE)
  pdf(file.path(root,"figures/fig01_vix_regimes.pdf"),width=9,height=3.8)
  cols<-c(calm="#197278",medium="#E6AA38",stress="#B33939")
  par(mar=c(4.2,4.3,.6,.7));plot(pits$dates,pits$v,type="l",col="#AAB3B8",xlab="Date",ylab="Lagged VIX level")
  points(pits$dates,pits$v,pch=16,cex=.25,col=cols[states])
  abline(v=as.numeric(as.Date(c(cfg$train_end,cfg$valid_end))),lty=2,col="#313D47")
  legend("topleft",names(cols),col=cols,pch=16,bty="n",horiz=TRUE);dev.off()
  p<-daily[daily$model=="WRR_direct",]
  pdf(file.path(root,"figures/fig02_joint_event_probability.pdf"),width=9,height=3.8)
  par(mar=c(4.2,4.3,.6,.7));plot(p$date,p$probability,type="l",xlab="Date",
    ylab="Six-ETF joint distress probability",col="#197278")
  points(p$date[p$observed==1],p$probability[p$observed==1],pch=16,col="#B33939");dev.off()
  invisible(NULL)
}
wrr_cross_universe <- function(a,b,cfg) {
  rows<-list()
  for(m in intersect(unique(a$model),unique(b$model)))for(metric in c("energy_common","variogram_common","brier","quantile_loss")) {
    x<-a[a$model==m,];y<-b[b$model==m,];j<-match(x$date,y$date)
    if(anyNA(j))stop("Cross-universe dates differ")
    d<-y[[metric]][j]-x[[metric]]
    bs<-wrr_with_seed(wrr_seed(cfg$seed,m,metric,"universes"),replicate(cfg$inference_B,
      mean(d[wrr_block_indices(length(d),min(20,length(d)))])))
    ci<-quantile(bs,c(.025,.975),names=FALSE);h<-wrr_hac(d,cfg$hac_lag)
    rows[[length(rows)+1L]]<-data.frame(model=m,metric=metric,n=length(d),
      extended_minus_main=mean(d),p=h["p"],ci_low=ci[1],ci_high=ci[2])
  }
  out<-do.call(rbind,rows);out$p_holm<-p.adjust(out$p,"holm");out
}
wrr_tex_escape <- function(x) {
  x<-gsub("_","\\_",as.character(x),fixed=TRUE)
  x<-gsub("%","\\%",x,fixed=TRUE);x<-gsub("&","\\&",x,fixed=TRUE);x
}
wrr_tex_table <- function(tab,file,digits=4) {
  vals<-lapply(tab,function(x)if(is.numeric(x))ifelse(is.finite(x),formatC(x,digits=digits,format="f"),"--") else wrr_tex_escape(x))
  t<-as.data.frame(vals,stringsAsFactors=FALSE)
  lines<-c(paste0("\\begin{tabular}{",paste0(rep("l",ncol(t)),collapse=""),"}"),"\\toprule",
    paste0(paste(wrr_tex_escape(names(t)),collapse=" & ")," \\\\"),"\\midrule",
    apply(t,1,function(r)paste0(paste(r,collapse=" & ")," \\\\")),"\\bottomrule","\\end{tabular}")
  dir.create(dirname(file),recursive=TRUE,showWarnings=FALSE);writeLines(lines,file)
}
wrr_exports <- function(results,pits,cfg) {
  root<-file.path(cfg$out_dir,"paper");dir.create(root,recursive=TRUE,showWarnings=FALSE)
  sm<-table(factor(pits$split,levels=c("train","validation","test")))
  macros<-c(paste0("\\newcommand{\\NAligned}{",nrow(pits$U),"}"),
    paste0("\\newcommand{\\NTrain}{",sm[1],"}"),paste0("\\newcommand{\\NValidation}{",sm[2],"}"),
    paste0("\\newcommand{\\NTest}{",nrow(results$main_no_vix$daily[results$main_no_vix$daily$model=="WRR_direct",]),"}"),
    paste0("\\newcommand{\\SampleFirst}{",min(pits$dates),"}"),paste0("\\newcommand{\\SampleLast}{",max(pits$dates),"}"))
  writeLines(macros,file.path(root,"sample_macros.tex"))
  for(u in names(results)) {
    r<-results[[u]];t<-r$summary
    wrr_tex_table(t[c("model","mean_log_score","energy_common","variogram_common","brier")],
      file.path(root,paste0(u,"_comparison.tex")))
    inf<-r$inference
    keep<-inf$block==min(20,nrow(r$daily[r$daily$model=="WRR_direct",])) & inf$metric%in%c("nll","energy_common","brier") &
      inf$benchmark%in%c("gaussian","t_copula","dcc_gaussian","static_vine","regime_vine","rolling_vine","rolling_matched")
    wrr_tex_table(inf[keep,c("benchmark","metric","difference","p_holm","ci_low","ci_high")],
      file.path(root,paste0(u,"_inference.tex")))
  }
  x<-results$main_no_vix$summary;y<-results$extended_with_vix$summary
  r<-x[x$model=="WRR_direct",];e<-y[y$model=="WRR_direct",]
  narrative<-sprintf(paste0("On the test sample, WRR-TVine has a mean copula log score of %.4f in the six-variable universe ",
    "and %.4f in the seven-variable universe. Its six-ETF Brier scores are %.5f and %.5f, respectively. ",
    "The density scores refer to different dimensions and must not be compared across universes. ",
    "Paired inference and matched ablations, rather than these point estimates alone, determine the scope of the contribution."),
    r$mean_log_score,e$mean_log_score,r$brier,e$brier)
  writeLines(narrative,file.path(root,"results_summary.tex"))
  invisible(root)
}
wrr_run <- function(cfg=wrr_config("full")) {
  wrr_check_packages(cfg);dir.create(cfg$out_dir,recursive=TRUE,showWarnings=FALSE)
  old_hash_dir <- options(wrr.hash_dir = file.path(cfg$out_dir, "cache", "hash_tmp"))
  on.exit(options(old_hash_dir), add = TRUE)
  started<-proc.time()[3];finished<-FALSE
  cfg$package_versions<-c(R=R.version.string,sapply(c("VineCopula","copula",if(cfg$margin=="garch_select")"rugarch"),function(p)as.character(packageVersion(p))))
  env<-environment(wrr_run);fn<-ls(envir=env,pattern="^wrr_")
  fn<-fn[vapply(fn,function(n)is.function(get(n,envir=env)),logical(1))]
  # Hash canonical code text, not mutable source-reference environments.
  cfg$engine_signature<-wrr_hash(lapply(fn,function(n)deparse(body(get(n,envir=env)),width.cutoff=500L)))
  on.exit({
    writeLines(paste("Elapsed seconds:",proc.time()[3]-started),file.path(cfg$out_dir,"elapsed_time.txt"))
    capture.output(sessionInfo(),file=file.path(cfg$out_dir,"sessionInfo.txt"))
    if(!finished)writeLines(c("INCOMPLETE: a stage failed. Completed checkpoints are retained.",
       "Read the R console/run log; fix the failure and rerun the same configuration."),file.path(cfg$out_dir,"RUN_STATUS.txt"))
  },add=TRUE)
  if(!isTRUE(cfg$mode%in%c("full","smoke")))stop("Invalid mode")
  writeLines(c("Historical assessment: the 2024-2025 outcomes have previously been examined.",
    "Configuration selection uses 2022-2023 only. Never choose seeds or specifications from test outcomes.",
    paste("Mode:",cfg$mode),paste("Started:",Sys.time())),file.path(cfg$out_dir,"RUN_STATUS.txt"))
  saveRDS(cfg,file.path(cfg$out_dir,"configuration.rds"));dput(cfg,file=file.path(cfg$out_dir,"configuration.R"))
  wrr_note("Stage: input data and marginal calibration")
  p<-wrr_prices(cfg);dat<-wrr_data(p,cfg)
  key<-wrr_hash(list(dat,cfg$margin,cfg$margin_candidates,cfg$margin_window,cfg$seed,cfg$arma_maxit,cfg$garch_retry,cfg$engine_signature,cfg$package_versions))
  pits<-wrr_cache(file.path(cfg$out_dir,"cache/margins.rds"),key,function()wrr_margins(dat,cfg),cfg$resume)
  wrr_write_marginal_reports(pits,cfg)
  wrr_csv(pits$bounds,file.path(cfg$out_dir,"tables/preprocessing_bounds.csv"))
  wrr_csv(data.frame(date=pits$dates,split=pits$split,pits$U,check.names=FALSE),file.path(cfg$out_dir,"tables/marginal_PIT.csv"))
  wrr_csv(data.frame(date=pits$dates,split=pits$split,pits$raw,check.names=FALSE),file.path(cfg$out_dir,"tables/untrimmed_returns.csv"))
  wrr_note("Stage: PIT diagnostics")
  diag<-wrr_pit_diagnostics(pits,cfg);wrr_csv(diag,file.path(cfg$out_dir,"tables/PIT_calibration.csv"))
  wrr_csv(data.frame(split=c("raw_returns","train","validation","test","aligned"),
       n=c(nrow(dat$R),sum(pits$split=="train"),sum(pits$split=="validation"),sum(pits$split=="test"),nrow(pits$U))),
       file.path(cfg$out_dir,"tables/sample_sizes.csv"))
  desc<-do.call(rbind,lapply(colnames(pits$raw),function(a) {
    x<-pits$raw[,a];data.frame(asset=a,n=length(x),mean=mean(x),sd=sd(x),min=min(x),q05=quantile(x,.05),
                            median=median(x),q95=quantile(x,.95),max=max(x))
  }));wrr_csv(desc,file.path(cfg$out_dir,"tables/descriptive_returns.csv"))
  wrr_note("Stage: copula validation, test forecasts and uncertainty")
  results<-list(main_no_vix=wrr_universe(pits,cfg$common_assets,"main_no_vix",cfg),
      extended_with_vix=wrr_universe(pits,cfg$symbols,"extended_with_vix",cfg))
  wrr_note("Stage: common-target comparisons and table exports")
  cross<-wrr_cross_universe(results$main_no_vix$daily,results$extended_with_vix$daily,cfg)
  wrr_csv(cross,file.path(cfg$out_dir,"tables/common_six_cross_universe.csv"))
  wrr_exports(results,pits,cfg)
  wrr_tex_table(desc,file.path(cfg$out_dir,"paper/descriptive.tex"))
  wrr_tex_table(diag[diag$split=="test" & diag$q==.05,c("asset","n","hit_rate","hit_ci_low","hit_ci_high","calibration_holm_p")],
                file.path(cfg$out_dir,"paper/PIT_calibration.tex"))
  # Predeclared sensitivities: rebuild margins; no reuse of old empirical numbers.
  if(cfg$run_robustness) for(kind in c("empirical_margin","train_winsor")) {
    wrr_note("Sensitivity: ",kind);cc<-cfg;cc$out_dir<-file.path(cfg$out_dir,"sensitivity",kind)
    cc$run_audits<-FALSE;cc$run_robustness<-FALSE
    if(kind=="empirical_margin")cc$margin<-"empirical" else cc$preprocessing<-"train_winsor"
    dd<-wrr_data(p,cc);kk<-wrr_hash(list(dd,cc$margin,cc$margin_candidates,cc$margin_window,cc$seed,cc$arma_maxit,cc$garch_retry,cc$engine_signature,cc$package_versions))
    pp<-wrr_cache(file.path(cc$out_dir,"cache/margins.rds"),kk,function()wrr_margins(dd,cc),cc$resume)
    wrr_write_marginal_reports(pp,cc)
    wrr_csv(wrr_pit_diagnostics(pp,cc),file.path(cc$out_dir,"tables/PIT_calibration.csv"))
    for(uni in c("main_no_vix","extended_with_vix")) {
      aa<-if(uni=="main_no_vix")cc$common_assets else cc$symbols
      wrr_universe(pp,aa,uni,cc,light=TRUE)
    }
  }
  capture.output(sessionInfo(),file=file.path(cfg$out_dir,"sessionInfo.txt"))
  writeLines(c(paste("COMPLETED",cfg$mode,Sys.time()),
    if(cfg$mode=="smoke")"SMOKE OUTPUTS ARE SOFTWARE CHECKS, NOT PAPER RESULTS." else
       "Empirical tables generated. Inspect calibration and uncertainty diagnostics before interpretation."),
    file.path(cfg$out_dir,"RUN_STATUS.txt"))
  finished<-TRUE
  wrr_note("Finished: ",normalizePath(cfg$out_dir))
  invisible(results)
}

wrr_run_logged <- function(cfg,log_file=paste0(toupper(cfg$mode),"_RUN_CONSOLE.log"),runner=wrr_run) {
  dir.create(cfg$out_dir,recursive=TRUE,showWarnings=FALSE)
  lock<-file.path(cfg$out_dir,"RUN_LOCK")
  if(!dir.create(lock,showWarnings=FALSE))stop("Output folder is locked: ",lock,
    ". Do not run two writers. If an earlier session crashed, first verify it is stopped, then remove ONLY RUN_LOCK.")
  on.exit(unlink(lock,recursive=TRUE),add=TRUE)
  writeLines(c(paste("PID",Sys.getpid()),paste("Started",Sys.time())),file.path(lock,"owner.txt"))
  old_options<-options(wrr.progress_file=file.path(cfg$out_dir,"PROGRESS.txt"))
  on.exit(options(old_options),add=TRUE)
  con<-file(log_file,open="at",encoding="UTF-8");depth<-sink.number(type="output")
  sink(con,split=TRUE)
  old_log<-options(wrr.log_connection=con)
  on.exit({options(old_log);while(sink.number(type="output")>depth)sink();close(con)},add=TRUE)
  warning_file<-file.path(cfg$out_dir,"general_warnings.csv")
  record_warning<-function(w) {
    row<-data.frame(time=format(Sys.time(),"%Y-%m-%d %H:%M:%S"),
      message=conditionMessage(w),call=paste(deparse(conditionCall(w)),collapse=" "))
    write.table(row,warning_file,sep=",",row.names=FALSE,col.names=!file.exists(warning_file),
      append=file.exists(warning_file),qmethod="double",fileEncoding="UTF-8")
    wrr_note("WARNING [recorded]: ",conditionMessage(w));invokeRestart("muffleWarning")
  }
  wrr_note("START ",toupper(cfg$mode),"; output=",normalizePath(cfg$out_dir,winslash="/"))
  value<-tryCatch(withCallingHandlers(runner(cfg),warning=record_warning,message=function(m) {
      wrr_note(trimws(conditionMessage(m)));invokeRestart("muffleMessage")
    }),error=function(e) {wrr_note("ERROR: ",conditionMessage(e));stop(e)},
       interrupt=function(e) {wrr_note("INTERRUPTED: checkpoints retained");stop(e)})
  wrr_note("END ",toupper(cfg$mode),"; log=",normalizePath(log_file,winslash="/"))
  invisible(value)
}
wrr_package_results <- function(root=file.path("results","full"),log_file="FULL_RUN_CONSOLE.log",
                                archive="WRR_TVINE_FULL_RESULTS.zip") {
  if(dir.exists(file.path(root,"RUN_LOCK")))stop("A run is active or a lock remains; do not package partial outputs")
  status<-file.path(root,"RUN_STATUS.txt")
  if(!file.exists(status)||!any(grepl("^COMPLETED full",readLines(status))))stop("Full run has not completed")
  if(!requireNamespace("zip",quietly=TRUE))stop("Install zip before packaging: install.packages('zip')")
  root_abs<-normalizePath(root,winslash="/",mustWork=TRUE)
  archive_abs<-file.path(normalizePath(dirname(archive),winslash="/",mustWork=TRUE),basename(archive))
  if(startsWith(archive_abs,paste0(root_abs,"/")))stop("Keep the results ZIP outside the output folder")
  if(file.exists(log_file)) {
    src<-normalizePath(log_file,winslash="/",mustWork=TRUE)
    dst<-file.path(root_abs,"run_console.log")
    if(src!=dst&&!file.copy(src,dst,overwrite=TRUE))stop("Could not copy the console log into results")
  }
  relative<-list.files(root_abs,recursive=TRUE,full.names=FALSE)
  relative<-relative[!grepl("(^|/)cache/",relative)]
  files<-paste0(basename(root_abs),"/",relative)
  zip::zip(archive_abs,files=files,root=dirname(root_abs),include_directories=FALSE,mode="mirror")
  entries<-zip::zip_list(archive_abs)$filename
  stopifnot(!anyDuplicated(entries),setequal(entries,gsub("\\\\","/",files)))
  wrr_note("Results archive: ",archive_abs)
  invisible(archive_abs)
}
wrr_launch <- function(mode=c("full","smoke"),prices_csv="data/prices.csv") {
  mode<-match.arg(mode);cfg<-wrr_config(mode);cfg$prices_csv<-prices_csv
  log_file<-paste0(toupper(mode),"_RUN_CONSOLE.log")
  result<-wrr_run_logged(cfg,log_file)
  if(mode=="full")wrr_package_results(cfg$out_dir,log_file)
  invisible(result)
}

