.this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.root <- if (!is.null(.this_file) && length(.this_file) == 1L && nzchar(.this_file)) {
  dirname(normalizePath(.this_file))
} else {
  getwd()
}
options(sslmm.root = .root)
source(file.path(.root, "config.R"))

source(file.path(.root, "scripts", "01_parameter_recovery.R"))
source(file.path(.root, "scripts", "02_efficiency.R"))
source(file.path(.root, "scripts", "03_fixed_missingness.R"))
source(file.path(.root, "scripts", "04_latent_robustness.R"))

# Full Monte Carlo studies.
run_parameter_recovery()
run_efficiency()
run_fixed_missingness()
run_latent_robustness()

# The numerical leading-order ARE calculation is intentionally optional because
# numerical expected-information Hessians are expensive. Run separately with:
# source("scripts/05_theory_are.R")
# run_theory_are(n_mc = 20000)
