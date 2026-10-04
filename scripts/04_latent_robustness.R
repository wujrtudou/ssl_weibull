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

worker_latent_robustness <- function(cfg, rep_id, control = SIM_CONTROL) {
  d <- make_partially_labelled_data(cfg$n, cfg$theta, cfg$alpha, cfg$xi0, cfg$xi1,
                                    control$entropy_floor)
  pc <- fit_pc(d, n_starts = control$n_starts_pc, control = control$nlminb_control,
               entropy_floor = control$entropy_floor, seed = 1000L + rep_id)
  cc <- fit_cc(d)
  recorded <- fit_recorded_mixed(d, n_starts = control$n_starts_recorded,
                                 control = control$nlminb_control,
                                 entropy_floor = control$entropy_floor,
                                 pc_fit = pc, compute_hessian = FALSE,
                                 seed = 2000L + rep_id)
  latent <- fit_latent_mixed(d, alpha_starts = control$latent_alpha_starts,
                             control = control$nlminb_control,
                             entropy_floor = control$entropy_floor,
                             alpha_upper = control$alpha_upper,
                             pc_fit = pc, compute_hessian = FALSE,
                             seed = 3000L + rep_id)
  oracle <- fit_oracle_latent(d, alpha = cfg$alpha,
                             xi = c(xi0 = cfg$xi0, xi1 = cfg$xi1),
                             n_starts = control$n_starts_pc,
                             control = control$nlminb_control,
                             entropy_floor = control$entropy_floor,
                             pc_fit = pc, seed = 4000L + rep_id)
  mar <- fit_mar(d, n_starts = control$n_starts_recorded,
                 control = control$nlminb_control,
                 entropy_floor = control$entropy_floor,
                 pc_fit = pc, seed = 5000L + rep_id)
  mcar <- pc

  c_cc <- classification_metrics(cc, cfg$theta)
  c_rec <- classification_metrics(recorded, cfg$theta)
  c_lat <- classification_metrics(latent, cfg$theta)
  c_orc <- classification_metrics(oracle, cfg$theta)
  c_mar <- classification_metrics(mar, cfg$theta)
  c_mca <- classification_metrics(mcar, cfg$theta)

  Dr <- composite_curve_rmse(latent, cfg$theta, cfg$alpha,
                             c(xi0 = cfg$xi0, xi1 = cfg$xi1),
                             n_eval = control$latent_curve_eval_n,
                             entropy_floor = control$entropy_floor)
  bi <- sqrt(mean((latent$beta - theta_to_beta(cfg$theta))^2))
  si <- start_instability(latent)

  list(rep_id = rep_id,
       latent_alpha = latent$alpha,
       latent_xi0 = latent$xi["xi0"], latent_xi1 = latent$xi["xi1"],
       recorded_alpha = recorded$alpha,
       recorded_xi0 = recorded$xi["xi0"], recorded_xi1 = recorded$xi["xi1"],
       Dr = Dr,
       beta_true = theta_to_beta(cfg$theta),
       latent_beta = latent$beta,
       latent_beta_rmse = bi,
       cc_error = c_cc["error"], recorded_error = c_rec["error"],
       oracle_error = c_orc["error"], latent_error = c_lat["error"],
       mar_error = c_mar["error"], mcar_error = c_mca["error"],
       cc_excess = c_cc["excess"], recorded_excess = c_rec["excess"],
       oracle_excess = c_orc["excess"], latent_excess = c_lat["excess"],
       mar_excess = c_mar["excess"], mcar_excess = c_mca["excess"],
       latent_convergence = latent$convergence,
       recorded_convergence = recorded$convergence,
       loglik_start_spread = si["loglik_spread"],
       alpha_start_spread = si["alpha_spread"])
}

run_latent_robustness <- function(B = SIM_CONTROL$B, n_cores = SIM_CONTROL$n_cores,
                                  overwrite = FALSE, control = SIM_CONTROL) {
  cfgs <- lapply(latent_configs(), prepare_config_truth, recovery = FALSE)
  raw_dir <- ensure_dir(file.path(.root, "results", "raw", "latent_robustness"))
  sum_dir <- ensure_dir(file.path(.root, "results", "summary"))
  summaries <- list()
  beta_summaries <- list()
  for (i in seq_along(cfgs)) {
    cfg <- cfgs[[i]]
    message("Running ", cfg$id)
    out <- file.path(raw_dir, paste0(cfg$id, ".rds"))
    res <- checkpoint_mc(cfg, B, worker_latent_robustness, out,
                         n_cores = n_cores, batch_size = control$batch_size,
                         seed_base = control$seed + 30000000L + i * 100000L,
                         overwrite = overwrite, control = control)
    summaries[[i]] <- summarise_latent_config(res, cfg)
    beta_summaries[[i]] <- summarise_latent_beta_config(res, cfg)
  }
  sm <- do.call(rbind, summaries)
  write.csv(sm, file.path(sum_dir, "latent_robustness_summary.csv"), row.names = FALSE)
  bsm <- do.call(rbind, beta_summaries)
  write.csv(bsm, file.path(sum_dir, "latent_beta_recovery.csv"), row.names = FALSE)
  invisible(sm)
}
