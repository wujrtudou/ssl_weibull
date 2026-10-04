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

worker_fixed_missingness <- function(cfg, rep_id, control = SIM_CONTROL) {
  d <- make_partially_labelled_data(cfg$n, cfg$theta, cfg$alpha, cfg$xi0, cfg$xi1,
                                    control$entropy_floor)
  pc0 <- fit_pc(d, n_starts = control$n_starts_pc, control = control$nlminb_control,
                entropy_floor = control$entropy_floor, seed = 1000L + rep_id)
  mixed <- fit_recorded_mixed(d, n_starts = control$n_starts_recorded,
                              control = control$nlminb_control,
                              entropy_floor = control$entropy_floor,
                              pc_fit = pc0, compute_hessian = FALSE,
                              seed = 2000L + rep_id)
  cc <- fit_cc(d)
  cm <- classification_metrics(mixed, cfg$theta)
  ccmet <- classification_metrics(cc, cfg$theta)
  list(rep_id = rep_id,
       mixed_excess = cm["excess"], cc_excess = ccmet["excess"],
       mixed_error = cm["error"], cc_error = ccmet["error"],
       mixed_convergence = mixed$convergence)
}

run_fixed_missingness <- function(B = SIM_CONTROL$B, n_cores = SIM_CONTROL$n_cores,
                                  overwrite = FALSE, control = SIM_CONTROL,
                                  bootstrap_R = 2000L) {
  cfgs <- lapply(fixed_missing_configs(), prepare_config_truth, recovery = FALSE)
  raw_dir <- ensure_dir(file.path(.root, "results", "raw", "fixed_missingness"))
  sum_dir <- ensure_dir(file.path(.root, "results", "summary"))
  summaries <- list()
  for (i in seq_along(cfgs)) {
    cfg <- cfgs[[i]]
    message("Running ", cfg$id, "; MAR gamma=", signif(cfg$gamma, 4))
    out <- file.path(raw_dir, paste0(cfg$id, ".rds"))
    res <- checkpoint_mc(cfg, B, worker_fixed_missingness, out,
                         n_cores = n_cores, batch_size = control$batch_size,
                         seed_base = control$seed + 20000000L + i * 100000L,
                         overwrite = overwrite, control = control)
    x <- summarise_efficiency_config(res, cfg, bootstrap_R = bootstrap_R)
    x$Gamma <- cfg$Gamma
    x$gamma_mar <- cfg$gamma
    summaries[[i]] <- x
  }
  sm <- do.call(rbind, summaries)
  write.csv(sm, file.path(sum_dir, "fixed_missingness_summary.csv"), row.names = FALSE)
  invisible(sm)
}
