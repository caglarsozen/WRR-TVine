# Computational environment

The recorded reference experiment used R 4.5.0 on Windows with `VineCopula` 2.6.1, `copula` 1.1-6, and `rugarch` 1.5-4. The complete recorded R session is in [reference_session_info.txt](reference_session_info.txt). The input snapshot fingerprint is recorded separately.

`00_INSTALL_PACKAGES.R` installs missing packages from CRAN and prints the installed versions. It does not replace already installed packages merely because their releases differ. Select an appropriate R environment and archived package releases if exact numerical reproduction is required. Optimizer, platform, BLAS/LAPACK and market-data differences can affect numerical results.

The analysis utilities require Python 3 plus the packages listed in [analysis/requirements.txt](../analysis/requirements.txt). Every full run records its own configuration, package information, session details, elapsed time and completion status in `results/full`.
