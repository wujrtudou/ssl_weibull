.this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
.root <- if (!is.null(.this_file) && length(.this_file) == 1L && nzchar(.this_file)) {
  dirname(dirname(normalizePath(.this_file)))
} else {
  getwd()
}
options(sslmm.root = .root)
source(file.path(.root, "config.R"))
source(file.path(.root, "R", "00_utils.R"))
source_sslmm(.root)

run_theory_are <- function(n_mc = 20000L, overwrite = FALSE, control = SIM_CONTROL) {
  sum_dir <- ensure_dir(file.path(.root, "results", "summary"))
  out <- file.path(sum_dir, "theory_are_numerical.csv")
  if (file.exists(out) && !overwrite) return(read.csv(out))

  grid <- expand.grid(p = c(1L,5L), Delta = c(1,2,3),
                      alpha = c(0,0.10,0.20,0.30,0.40), KEEP.OUT.ATTRS = FALSE)
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    p <- grid$p[i]; Delta <- grid$Delta[i]; alpha <- grid$alpha[i]
    theta <- make_theta(p, Delta, Sigma = diag(p), pi1 = 0.5)
    xi1 <- 3; gamma <- 0.30
    xi0 <- calibrate_xi0(theta, xi1, gamma)
    message(sprintf("Numerical ARE: p=%d Delta=%g alpha=%g", p, Delta, alpha))
    are <- try(estimate_leading_are(theta, xi0, xi1, alpha, n_mc = n_mc,
                                    seed = control$seed + i * 1000L,
                                    entropy_floor = control$entropy_floor), silent = TRUE)
    if (inherits(are, "try-error")) are <- NA_real_
    rows[[i]] <- data.frame(p=p, Delta=Delta, alpha=alpha, gamma=gamma, xi1=xi1,
                            xi0=xi0, ARE=are, n_mc=n_mc)
  }
  ans <- do.call(rbind, rows)
  write.csv(ans, out, row.names = FALSE)
  ans
}
