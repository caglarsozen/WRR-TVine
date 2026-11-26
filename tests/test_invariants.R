# Run from the package folder: source("tests/test_invariants.R")
source("R/wrr_tvine.R", local = TRUE, encoding = "UTF-8")
checks <- character()
check <- function(x,label) {stopifnot(isTRUE(x)); checks<<-c(checks,label);cat("PASS",label,"\n")}
cfg<-wrr_config("smoke");cfg$common_assets<-c("A","B","C");cfg$k<-2L
set.seed(2);U<-matrix(runif(900),300,3,dimnames=list(NULL,c("A","B","C")))
v<-20+sin(seq_len(300)/13);states<-ifelse(v>20,"stress","calm");past<-100:199;t<-200L
w<-wrr_weights(U,v,states,past,t,cfg)
check(abs(sum(w$p)-1)<1e-12 && all(w$p>0),"normalized positive weights")
V<-wrr_weighted_ranks(U[past,],w$p)
check(all(V>0 & V<1),"weighted PIT support")
check(all(vapply(1:3,function(j)wrr_weighted_ks(V[,j],w$p)<=max(w$p)/2+1e-10,logical(1))),"weighted midpoint CDF bound")
U2<-U;U2[t:nrow(U),]<-1-U2[t:nrow(U),]
check(identical(w,wrr_weights(U2,v,states,past,t,cfg)),"weights ignore current and future returns")
v2<-v;v2[(t+1):length(v2)]<-999
check(identical(w,wrr_weights(U,v2,states,past,t,cfg)),"weights ignore future VIX")
state2<-states;state2[(t+1):length(v2)]<-"other"
check(identical(w,wrr_weights(U,v,state2,past,t,cfg)),"weights ignore future states")
set.seed(1);rng<-.Random.seed;tmp<-wrr_with_seed(123,runif(8))
check(identical(rng,.Random.seed),"local RNG restores caller state")
check(identical(tmp,wrr_with_seed(123,runif(8))),"deterministic independent RNG streams")
check(abs(wrr_gaussian_ll(matrix(c(.1,.3,.9),1),diag(3)))<1e-12,"independence Gaussian density equals one")
X<-matrix(c(0,1,2,3),ncol=1)
check(isTRUE(all.equal(wrr_energy(X,1),mean(abs(X-1))-.5*mean(c(2,2)))),"energy score disjoint-pair estimator")
ids<-wrr_block_indices(31,10)
check(length(ids)==31 && all(ids%in%1:31),"circular block bootstrap length and support")
E<-rbind(c(.01,.02,.99),c(.01,.9,.9));colnames(E)<-cfg$common_assets
check(identical(wrr_event(E,cfg$common_assets,.1,2L),c(1L,0L)),"joint crash event count")
m<-wrr_empirical_margin(seq_len(60),20)
check(all(is.na(m$u[1:20])) && all(m$u[21:60]<1),"rolling empirical PIT burn-in and endpoints")
m2<-wrr_empirical_margin(c(seq_len(50),rep(1e6,10)),20)
check(identical(m$u[1:50],m2$u[1:50]),"empirical margins have prefix invariance")
p<-tempfile(fileext="rds");val<-wrr_cache(p,"first",function()1)
check(identical(wrr_cache(p,"second",function()2),2),"checkpoint key invalidation");unlink(p)
tied<-matrix(c(.1,.1,.2,.4),ncol=1);q<-wrr_weighted_ranks(tied,c(.1,.2,.3,.4))[,1]
check(max(abs(q-c(.15,.15,.45,.8)))<1e-12,"weighted ranks handle ties exactly")
# CDF-tilting counterexample: changing historical weights can distort margins.
check(wrr_weighted_ks(c(.1,.3,.6,.9),c(.7,.1,.1,.1))>.3,"weighted marginal distortion is detectable")

# All thresholds belong to the declared family, including undefined tests.
pit_seed <- wrr_with_seed(321L, matrix(runif(80L * 7L, .02, .98), ncol = 7L))
pit_seed[, 1L] <- 0.22 + 0.76 * pit_seed[, 1L]
colnames(pit_seed) <- LETTERS[1:7]
pit_check <- wrr_pit_diagnostics(
  list(U = pit_seed, Z = qnorm(pit_seed), split = rep(c("validation", "test"), each = 40L)),
  cfg
)
constant <- pit_check$hit_rate %in% c(0, 1)
check(
  all(is.na(pit_check$calibration_hac_p[constant])) &&
    all(is.na(pit_check$calibration_holm_p[constant])) &&
    all(pit_check$hit_ci_low[constant] == pit_check$hit_rate[constant]) &&
    all(pit_check$hit_ci_high[constant] == pit_check$hit_rate[constant]) &&
    all(is.finite(pit_check$calibration_hac_p[!constant])),
  "constant PIT hits retain frequencies and intervals without HAC significance"
)
for (period in c("validation", "test")) {
  period_pit <- pit_check[pit_check$split == period, ]
  p <- period_pit$calibration_hac_p
  adjusted <- p.adjust(p, "holm", n = 49L)
  reduced <- p.adjust(p, "holm", n = sum(!is.na(p)))
  check(
    nrow(period_pit) == 49L && any(adjusted > reduced, na.rm = TRUE) &&
      identical(period_pit$calibration_holm_p, adjusted),
    paste("PIT Holm adjustment retains the 49-test family in", period)
  )
}

# Nondegenerate event tilting and the corresponding proper-score penalties.
b <- c(.10,.20,.30,.40); K <- c(0,1,2,3); A <- as.integer(K>=2)
a <- sum(b*A); lambda <- 2; tilted <- b*(1+lambda*A)/(1+lambda*a)
a_plus <- sum(tilted*A)
check(abs(a_plus-(1+lambda)*a/(1+lambda*a))<1e-12 &&
      sum(tilted*K)>sum(b*K),"event tilting increases event and aggregate hit mass")
kl <- sum(b*log(b/tilted))
check(kl>0 && abs(kl-(log1p(lambda*a)-a*log1p(lambda)))<1e-12,
      "positive logarithmic-loss penalty under the base law")
check(abs(sum(b*((a_plus-A)^2-(a-A)^2))-(a_plus-a)^2)<1e-12,
      "squared event-probability loss identity")
# A cache may contain a damaged file after an interrupted write.
cache <- tempfile(fileext=".rds")
wrr_cache(cache,"key",function()42)
check(identical(wrr_cache(cache,"key",function()stop("unexpected computation")),42),
      "compatible checkpoint reuse")
writeLines("incomplete checkpoint",cache)
check(identical(wrr_cache(cache,"key",function()43),43),"unreadable checkpoint recovery")
unlink(cache)
# Hashes remain stable when the dedicated hashing directory has to be recreated.
hashdir <- tempfile(pattern="wrr_hash_test_")
old_hash <- options(wrr.hash_dir=hashdir)
hash1 <- wrr_hash(list(a=seq_len(5),b=pi))
unlink(hashdir,recursive=TRUE)
hash2 <- wrr_hash(list(a=seq_len(5),b=pi))
check(identical(hash1,hash2),"stable checkpoint hash after directory recreation")
options(old_hash);unlink(hashdir,recursive=TRUE)
full <- wrr_config("full"); smoke <- wrr_config("smoke")
check(nrow(full$grid)==12L && full$margin_window==750L &&
      full$scenarios==10000L && full$estimation_B==199L && full$inference_B==1999L &&
      identical(full$ensemble_B,c(1L,5L,20L)) && full$run_robustness,
      "full experiment settings")
check(full$out_dir!=smoke$out_dir && smoke$max_test==8L && !smoke$run_robustness,
      "smoke outputs remain separate from full experiment")
cat("TOTAL",length(checks),"checks passed\n")
