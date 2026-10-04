DEFAULT_CORES <- parallel::detectCores(logical = FALSE)
if (is.na(DEFAULT_CORES)) DEFAULT_CORES <- 1L

SIM_CONTROL <- list(
  B = 500L,
  n_cores = min(8L, max(1L, DEFAULT_CORES - 1L)),
  batch_size = 25L,
  seed = 20260807L,
  n_starts_recorded = 3L,
  n_starts_pc = 3L,
  latent_alpha_starts = c(0, 0.05, 0.15, 0.30),
  test_n = 10000L,
  latent_curve_eval_n = 5000L,
  compute_hessian_recovery = TRUE,
  nlminb_control = list(eval.max = 2000L, iter.max = 1000L, rel.tol = 1e-9, x.tol = 1e-8),
  entropy_floor = 1e-300,
  alpha_upper = 1 - 1e-8
)
