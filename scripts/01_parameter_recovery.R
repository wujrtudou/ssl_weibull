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

worker_parameter_recovery <- function(cfg, rep_id, control = SIM_CONTROL) {
  d <- make_partially_labelled_data(cfg$n, cfg$theta, cfg$alpha, cfg$xi0, cfg$xi1,
                                    control$entropy_floor)
  pc <- fit_pc(d, n_starts = control$n_starts_pc, control = control$nlminb_control,
               entropy_floor = control$entropy_floor, seed = 1000L + rep_id)
  mixed <- fit_recorded_mixed(d, n_starts = control$n_starts_recorded,
                              control = control$nlminb_control,
                              entropy_floor = control$entropy_floor,
                              pc_fit = pc, compute_hessian = control$compute_hessian_recovery,
                              seed = 2000L + rep_id)
  cc <- fit_cc(d)
  test <- simulate_gaussian2(control$test_n, cfg$theta)

  tr <- truth_mixed_vector(cfg$theta, cfg$xi0, cfg$xi1, cfg$alpha)
  est <- reported_mixed_vector(mixed)
  sev <- recorded_se_vectors_with_n(mixed, cfg$p, cfg$n)
  se <- sev$parameter
  if (is.null(se)) se <- setNames(rep(NA_real_, length(est)), names(est))
  se <- se[names(est)]

  beta_true <- theta_to_beta(cfg$theta)
  beta_est <- mixed$beta
  beta_se <- sev$beta
  if (is.null(beta_se)) beta_se <- setNames(rep(NA_real_, length(beta_true)), names(beta_true))

  cmix <- classification_metrics(mixed, cfg$theta, test)
  cpc <- classification_metrics(pc, cfg$theta, test)
  ccc <- classification_metrics(cc, cfg$theta, test)

  list(
    rep_id = rep_id,
    truth_param = tr,
    est_param = est,
    se_param = se,
    beta_truth = beta_true,
    beta_est = beta_est,
    beta_se = beta_se,
    mixed_error = cmix["error"], mixed_test_error = cmix["test_error"],
    pc_error = cpc["error"], pc_test_error = cpc["test_error"],
    cc_error = ccc["error"], cc_test_error = ccc["test_error"],
    mixed_convergence = mixed$convergence,
    pc_convergence = pc$convergence,
    valid_information = !is.null(mixed$vcov_internal),
    loglik = mixed$loglik
  )
}

run_parameter_recovery <- function(B = SIM_CONTROL$B, n_cores = SIM_CONTROL$n_cores,
                                   overwrite = FALSE, control = SIM_CONTROL) {
  cfgs <- lapply(parameter_recovery_configs(), prepare_config_truth, recovery = TRUE)
  raw_dir <- ensure_dir(file.path(.root, "results", "raw", "parameter_recovery"))
  sum_dir <- ensure_dir(file.path(.root, "results", "summary"))
  summaries <- list()
  class_summaries <- list()

  for (i in seq_along(cfgs)) {
    cfg <- cfgs[[i]]
    message("Running ", cfg$id)
    out <- file.path(raw_dir, paste0(cfg$id, ".rds"))
    res <- checkpoint_mc(cfg, B, worker_parameter_recovery, out,
                         n_cores = n_cores, batch_size = control$batch_size,
                         seed_base = control$seed + i * 100000L,
                         overwrite = overwrite, control = control)
    summaries[[i]] <- summarise_recovery_config(res, cfg)
    class_summaries[[i]] <- summarise_recovery_classification_config(res, cfg)
  }
  sm <- do.call(rbind, summaries)
  # Do not interpret regular Wald coverage for alpha at the boundary.
  sm$coverage95[sm$parameter == "alpha" & sm$truth == 0] <- NA_real_
  write.csv(sm, file.path(sum_dir, "parameter_recovery_summary.csv"), row.names = FALSE)
  csm <- do.call(rbind, class_summaries)
  write.csv(csm, file.path(sum_dir, "parameter_recovery_classification.csv"), row.names = FALSE)
  invisible(sm)
}
