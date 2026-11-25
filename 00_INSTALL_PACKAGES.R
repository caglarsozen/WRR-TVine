# Run once in R / RStudio. Only declared required dependencies are installed.
options(repos=c(CRAN="https://cloud.r-project.org"),timeout=600)
pkgs<-c("quantmod","rugarch","VineCopula","copula","zip")
missing<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))install.packages(missing,dependencies=c("Depends","Imports","LinkingTo"))
stopifnot(all(vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)))
print(sapply(pkgs,function(p)as.character(packageVersion(p))))
