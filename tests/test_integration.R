# Synthetic integration test only; does not establish empirical forecasting accuracy.
source("R/wrr_tvine.R", local = TRUE, encoding = "UTF-8")
cfg<-wrr_config("smoke");cfg$margin<-"empirical";cfg$margin_window<-40L
cfg$family_tail<-cfg$family_general<-c(0L,1L)
cfg$grid<-data.frame(window=120L,refit=2L,family="tail")
cfg$dcc_grid<-data.frame(a=.03,b=.94);cfg$regime_sensitivity<-list(c(.55,.80))
cfg$inference_B<-29L;cfg$estimation_B<-2L;cfg$audit_origins<-2L;cfg$max_test<-4L
cfg$run_robustness<-FALSE;cfg$download<-FALSE
cfg$out_dir<-"test_outputs_synthetic";cfg$prices_csv<-file.path(cfg$out_dir,"synthetic_prices.csv")
cfg$start<-"2019-01-01";cfg$end<-"2022-01-01"
dates<-seq(as.Date("2019-01-01"),by="day",length.out=430L)
cfg$train_end<-as.character(dates[301]);cfg$valid_end<-as.character(dates[381]);cfg$test_start<-as.character(dates[382])
set.seed(513);f<-rnorm(430);ret<-sapply(1:6,function(j).01*(.45*f+rnorm(430)))
prices<-100*exp(apply(ret,2,cumsum));colnames(prices)<-cfg$common_assets
p<-data.frame(date=dates,prices,VIX=exp(3+.2*sin(1:430/19)+.04*rnorm(430)))
wrr_csv(p,cfg$prices_csv)
r<-wrr_run(cfg)
stopifnot(nrow(r$main_no_vix$daily)>0,all(is.finite(r$main_no_vix$daily$nll)))
stopifnot(all(r$main_no_vix$daily$observed==rep(r$main_no_vix$daily$observed[seq(1,nrow(r$main_no_vix$daily),length(unique(r$main_no_vix$daily$model)))],each=length(unique(r$main_no_vix$daily$model)))))
stopifnot(all(c("WRR_direct","WRR_no_tail","rolling_vine","WRR_resampled_B2")%in%r$main_no_vix$summary$model))
# Exact repeat must reuse the cached objects, as well as preserve scores.
cache_files<-list.files(cfg$out_dir,pattern="[.]rds$",recursive=TRUE,full.names=TRUE)
cache_files<-cache_files[grepl("/cache/",cache_files)]
cache_times<-file.info(cache_files)$mtime
r2<-wrr_run(cfg)
stopifnot(identical(r$main_no_vix$daily,r2$main_no_vix$daily))
stopifnot(identical(r$extended_with_vix$daily,r2$extended_with_vix$daily))
stopifnot(identical(cache_times,file.info(cache_files)$mtime))
cat("SYNTHETIC INTEGRATION PASSED; these are not empirical results.\n")
