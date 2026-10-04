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

run_smoke_test <- function() {
  ctl <- SIM_CONTROL
  ctl$n_starts_pc <- 1L
  ctl$n_starts_recorded <- 1L
  ctl$latent_alpha_starts <- c(0, 0.1)
  ctl$test_n <- 500L
  ctl$latent_curve_eval_n <- 500L
  ctl$nlminb_control <- list(eval.max=500L, iter.max=300L, rel.tol=1e-7, x.tol=1e-7)

  theta <- make_theta(2, 2, Sigma = diag(2), pi1 = .5)
  xi1 <- 3
  xi0 <- calibrate_xi0(theta, xi1, .30)
  d <- make_partially_labelled_data(200, theta, .10, xi0, xi1)
  pc <- fit_pc(d, n_starts=1, control=ctl$nlminb_control)
  rec <- fit_recorded_mixed(d, n_starts=1, control=ctl$nlminb_control, pc_fit=pc)
  lat <- fit_latent_mixed(d, alpha_starts=c(0,.1), control=ctl$nlminb_control, pc_fit=pc)
  mar <- fit_mar(d, n_starts=1, control=ctl$nlminb_control, pc_fit=pc)
  ora <- fit_oracle_latent(d, .10, c(xi0=xi0, xi1=xi1), n_starts=1,
                           control=ctl$nlminb_control, pc_fit=pc)
  cc <- fit_cc(d)
  ans <- data.frame(method=c("CC","PC/MCAR","recorded","latent","MAR","oracle"),
                    convergence=c(cc$convergence,pc$convergence,rec$convergence,lat$convergence,
                                  mar$convergence,ora$convergence),
                    error=c(classification_metrics(cc,theta)["error"],
                            classification_metrics(pc,theta)["error"],
                            classification_metrics(rec,theta)["error"],
                            classification_metrics(lat,theta)["error"],
                            classification_metrics(mar,theta)["error"],
                            classification_metrics(ora,theta)["error"]))
  print(ans)
  invisible(ans)
}
