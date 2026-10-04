# ============================================================
# make_simulation_figures_complete.R
# Complete plotting script for the Monte Carlo simulation section
#
# Project layout assumed:
#
# C:/Users/uqjwu15/Desktop/ssl_mm_simulation_windows_8core/
#   R/
#     00_utils.R
#     01_model.R
#     02_missingness.R
#     03_likelihoods.R
#     05_metrics.R
#     06_theory_are.R
#   raw/
#     ... .rds files (subfolders are OK)
#   results/
#     efficiency_summary.csv
#     fixed_missingness_summary.csv
#     observed_ARE_for_plot.csv   [created automatically]
#     figures/                    [created automatically]
#
# Output:
#   results/figures/sim_finite_sample.pdf
#   results/figures/sim_efficiency.pdf
#   results/figures/sim_fixed_missingness.pdf
#   results/figures/sim_latent.pdf
#
# PNG versions (600 dpi) are also written.
# ============================================================


rm(list = ls())


# ============================================================
# 0. PROJECT PATHS
# ============================================================

ROOT <- "C:/Users/uqjwu15/Desktop/ssl_mm_simulation_windows_8core"

R_DIR       <- file.path(ROOT, "R")
RESULTS_DIR <- file.path(ROOT, "results")
RAW_DIR     <- file.path(RESULTS_DIR, "raw")
FIG_DIR     <- file.path(RESULTS_DIR, "figures")

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(R_DIR)) {
  stop("R directory not found: ", R_DIR)
}

if (!dir.exists(RAW_DIR)) {
  stop("raw directory not found: ", RAW_DIR)
}


# ============================================================
# 1. PACKAGES
# ============================================================

pkgs <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "purrr",
  "stringr",
  "readr",
  "patchwork",
  "tibble"
)

missing_pkgs <- pkgs[
  !vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs)) {
  stop(
    "Please install these packages first: ",
    paste(missing_pkgs, collapse = ", ")
  )
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(readr)
library(patchwork)
library(tibble)


# ============================================================
# 2. SOURCE ALL R FILES
# ============================================================
#
# Your project stores the model, likelihood, fitter, metric,
# simulation-helper, and ARE functions under ROOT/R.
#
# We source every .R file in alphabetical order.  This avoids
# hard-coding names such as "01_model.R" when your actual file
# may be named, for example, "01_model(2).R".
# ============================================================

r_files <- sort(
  list.files(
    R_DIR,
    pattern = "\\.R$",
    full.names = TRUE
  )
)

if (!length(r_files)) {
  stop("No .R files found in: ", R_DIR)
}

message("Sourcing R files from: ", R_DIR)

for (f in r_files) {
  message("  source: ", basename(f))
  source(f, local = FALSE)
}

if (!exists("estimate_leading_are", mode = "function")) {
  stop(
    "After sourcing ROOT/R, estimate_leading_are() was not found. ",
    "Please check that 06_theory_are.R is present in: ", R_DIR
  )
}


# ============================================================
# 3. GENERAL HELPERS
# ============================================================

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

near_num <- function(x, y, tol = 1e-8) {
  abs(as.numeric(x) - as.numeric(y)) < tol
}

find_one_file <- function(filename) {
  
  candidates <- c(
    file.path(RESULTS_DIR, filename),
    file.path(ROOT, filename)
  )
  
  hit <- candidates[file.exists(candidates)]
  
  if (length(hit)) {
    return(hit[1])
  }
  
  all_files <- list.files(
    ROOT,
    recursive = TRUE,
    full.names = TRUE
  )
  
  hit <- all_files[basename(all_files) == filename]
  
  if (!length(hit)) {
    stop(
      "Cannot find ", filename,
      "\nSearched under: ", ROOT
    )
  }
  
  hit[1]
}


save_plot <- function(p, filename, width, height) {
  
  pdf_file <- file.path(
    FIG_DIR,
    paste0(filename, ".pdf")
  )
  
  png_file <- file.path(
    FIG_DIR,
    paste0(filename, ".png")
  )
  
  ggsave(
    filename = pdf_file,
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = cairo_pdf
  )
  
  ggsave(
    filename = png_file,
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 600
  )
  
  message("Saved: ", pdf_file)
  message("Saved: ", png_file)
}

theme_paper <- function(base_size = 10.5) {
  
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(
        fill = "grey95",
        colour = "black",
        linewidth = 0.4
      ),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "plain"),
      legend.position = "bottom",
      legend.title = element_text(face = "plain"),
      plot.tag = element_text(face = "bold"),
      plot.margin = margin(6, 8, 6, 6)
    )
}

# ------------------------------------------------------------
# Publication palette (color-blind friendly)
# ------------------------------------------------------------

COL_BLUE   <- "#0072B2"
COL_ORANGE <- "#D55E00"
COL_GREEN  <- "#009E73"
COL_PURPLE <- "#CC79A7"
COL_SKY    <- "#56B4E9"
COL_NAVY   <- "#003B5C"
COL_GREY   <- "#666666"
COL_LIGHT  <- "#E6E6E6"


# Rename a column to a canonical name if necessary.
rename_first_existing <- function(dat, target, candidates) {
  
  if (target %in% names(dat)) {
    return(dat)
  }
  
  hit <- candidates[candidates %in% names(dat)]
  
  if (!length(hit)) {
    stop(
      "Cannot find a column for '", target, "'.\n",
      "Available columns:\n",
      paste(names(dat), collapse = ", ")
    )
  }
  
  names(dat)[names(dat) == hit[1]] <- target
  dat
}


# ============================================================
# 4. FIGURE 1
# FINITE-SAMPLE CLASSIFICATION GAIN
#
# y = err(ignore) - err(observed)
#
# Positive values favor the full likelihood with observed
# missing-label indicators.
# ============================================================

# Parameter-recovery RDS files are stored under results/raw/parameter_recovery.
RECOVERY_RAW_DIR <- file.path(
  RAW_DIR,
  "parameter_recovery"
)

if (!dir.exists(RECOVERY_RAW_DIR)) {
  stop(
    "Parameter-recovery raw directory not found:\n",
    RECOVERY_RAW_DIR
  )
}

recovery_files <- list.files(
  RECOVERY_RAW_DIR,
  pattern = "^recovery_.*\\.rds$",
  full.names = TRUE
)

if (!length(recovery_files)) {
  stop(
    "No recovery_*.rds files found in:\n",
    RECOVERY_RAW_DIR
  )
}

message(
  "Found ", length(recovery_files),
  " parameter-recovery RDS files in ",
  RECOVERY_RAW_DIR
)


parse_recovery_name <- function(f) {
  
  nm <- basename(f)
  
  m <- str_match(
    nm,
    paste0(
      "^recovery_",
      "n([0-9]+)_",
      "D([0-9.]+)_",
      "g([0-9.]+)_",
      "x([0-9.]+)_",
      "a([0-9.]+)",
      "\\.rds$"
    )
  )
  
  if (any(is.na(m))) {
    stop("Unexpected recovery filename: ", nm)
  }
  
  tibble(
    file = f,
    n = as.integer(m[2]),
    Delta = as.numeric(m[3]),
    gamma = as.numeric(m[4]),
    xi1 = as.numeric(m[5]),
    alpha = as.numeric(m[6])
  )
}

extract_recovery_replications <- function(f) {
  
  meta <- parse_recovery_name(f)
  res <- readRDS(f)
  
  if (!is.list(res)) {
    stop("Expected a list of replications in: ", f)
  }
  
  vals <- map_dfr(
    seq_along(res),
    function(i) {
      
      x <- res[[i]]
      
      if (!is.null(x$error)) {
        return(NULL)
      }
      
      obs_err <- x$mixed_error %||% NA_real_
      ign_err <- x$pc_error %||% NA_real_
      
      if (!is.finite(obs_err) || !is.finite(ign_err)) {
        return(NULL)
      }
      
      tibble(
        replication = i,
        error_reduction = ign_err - obs_err
      )
    }
  )
  
  if (!nrow(vals)) {
    return(tibble())
  }
  
  bind_cols(
    meta[
      rep(1, nrow(vals)),
      c("n", "Delta", "gamma", "xi1", "alpha")
    ],
    vals
  )
}

recovery_raw <- map_dfr(
  recovery_files,
  extract_recovery_replications
)

if (!nrow(recovery_raw)) {
  stop("No usable recovery replications were found.")
}


# ---------- Panel A: sample size ----------

fig1_A <- recovery_raw %>%
  filter(
    near_num(Delta, 2),
    near_num(gamma, 0.30),
    near_num(xi1, 3),
    near_num(alpha, 0.10),
    n %in% c(500L, 1000L)
  ) %>%
  mutate(
    panel = "A. Sample size",
    setting = factor(
      paste0("n=", n),
      levels = c("n=500", "n=1000")
    )
  )


# ---------- Panel B: class separation ----------

fig1_B <- recovery_raw %>%
  filter(
    n == 1000L,
    near_num(gamma, 0.30),
    near_num(xi1, 3),
    near_num(alpha, 0.10),
    Delta %in% c(1, 2, 3)
  ) %>%
  mutate(
    panel = "B. Class separation",
    setting = factor(
      paste0("\u0394=", Delta),
      levels = c("\u0394=1", "\u0394=2", "\u0394=3")
    )
  )


# ---------- Panel C: entropy-MAR proportion ----------

fig1_C <- recovery_raw %>%
  filter(
    n == 1000L,
    near_num(Delta, 2),
    near_num(xi1, 3),
    near_num(alpha, 0.10),
    gamma %in% c(0.15, 0.30, 0.45)
  ) %>%
  mutate(
    panel = "C. Entropy-MAR proportion",
    setting = factor(
      sprintf("\u03B3=%.2f", gamma),
      levels = c(
        "\u03B3=0.15",
        "\u03B3=0.30",
        "\u03B3=0.45"
      )
    )
  )


# ---------- Panel D: entropy dependence ----------

fig1_D <- recovery_raw %>%
  filter(
    n == 1000L,
    near_num(Delta, 2),
    near_num(gamma, 0.30),
    near_num(alpha, 0.10),
    xi1 %in% c(1, 3, 5)
  ) %>%
  mutate(
    panel = "D. Entropy dependence",
    setting = factor(
      paste0("\u03BE\u2081=", xi1),
      levels = c(
        "\u03BE\u2081=1",
        "\u03BE\u2081=3",
        "\u03BE\u2081=5"
      )
    )
  )


# ---------- Panel E: MCAR component ----------

fig1_E <- recovery_raw %>%
  filter(
    n == 1000L,
    near_num(Delta, 2),
    near_num(gamma, 0.30),
    near_num(xi1, 3),
    alpha %in% c(0, 0.10, 0.30)
  ) %>%
  mutate(
    panel = "E. MCAR component",
    setting = factor(
      sprintf("\u03B1=%.2f", alpha),
      levels = c(
        "\u03B1=0.00",
        "\u03B1=0.10",
        "\u03B1=0.30"
      )
    )
  )


fig1_dat <- bind_rows(
  fig1_A,
  fig1_B,
  fig1_C,
  fig1_D,
  fig1_E
)

fig1 <- ggplot(
  fig1_dat,
  aes(
    x = setting,
    y = error_reduction,
    fill = panel
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 3,
    linewidth = 0.45
  ) +
  geom_boxplot(
    width = 0.62,
    linewidth = 0.42,
    outlier.size = 0.45
  ) +
  facet_wrap(
    ~ panel,
    scales = "free_x",
    ncol = 2
  ) +
  scale_fill_manual(
    values = c(
      "A. Sample size" = COL_BLUE,
      "B. Class separation" = COL_ORANGE,
      "C. Entropy-MAR proportion" = COL_GREEN,
      "D. Entropy dependence" = COL_PURPLE,
      "E. MCAR component" = COL_SKY
    ),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = expression(
      err(hat(beta)[ig]) -
        err(hat(beta)[obs])
    )
  ) +
  theme_paper() +
  theme(
    axis.text.x = element_text(
      angle = 25,
      hjust = 1
    )
  )

save_plot(
  fig1,
  "sim_finite_sample",
  width = 8.2,
  height = 7.8
)


# ============================================================
# 5. FIGURE 2
# OBSERVED MISSING-LABEL RELATIVE EFFICIENCY
#
# Main curves:
#   simulated RE, n = 1000
#   ARE
#
# Additional open points:
#   n = 500, Delta = 2 only
#
# No bootstrap intervals.
# ============================================================

eff_file <- find_one_file(
  "efficiency_summary.csv"
)

eff <- read_csv(
  eff_file,
  show_col_types = FALSE
)

eff <- rename_first_existing(
  eff,
  "Delta",
  c("Delta", "delta", "D")
)

eff <- rename_first_existing(
  eff,
  "simulated_RE",
  c(
    "simulated_RE",
    "RE",
    "relative_efficiency"
  )
)

required_eff_cols <- c(
  "n",
  "p",
  "Delta",
  "alpha",
  "simulated_RE"
)

if (!all(required_eff_cols %in% names(eff))) {
  stop(
    "efficiency_summary.csv is missing required columns.\n",
    "Required: ",
    paste(required_eff_cols, collapse = ", ")
  )
}

eff <- eff %>%
  mutate(
    n = as.integer(n),
    p = as.integer(p),
    Delta = as.numeric(Delta),
    alpha = as.numeric(alpha),
    simulated_RE = as.numeric(simulated_RE)
  )


# ------------------------------------------------------------
# Observed-experiment ARE
#
# This uses estimate_leading_are() from R/06_theory_are.R.
# Results are cached in results/observed_ARE_for_plot.csv.
# ------------------------------------------------------------

are_cache <- file.path(
  RESULTS_DIR,
  "observed_ARE_for_plot.csv"
)

if (file.exists(are_cache)) {
  
  are_obs <- read_csv(
    are_cache,
    show_col_types = FALSE
  )
  
} else {
  
  message(
    "Observed ARE cache not found. ",
    "Computing ARE values now..."
  )
  
  are_grid <- expand_grid(
    p = c(1L, 5L),
    Delta = c(1, 2, 3),
    alpha = c(
      0,
      0.10,
      0.20,
      0.30,
      0.40
    )
  )
  
  # Increase if desired for a more stable numerical approximation.
  ARE_N_MC <- 20000L
  ARE_SEED <- 20260813L
  
  are_obs <- pmap_dfr(
    are_grid,
    function(p, Delta, alpha) {
      
      theta_true <- make_theta(
        p = p,
        Delta = Delta,
        Sigma = diag(p),
        pi1 = 0.5
      )
      
      xi1 <- 3
      gamma <- 0.30
      
      xi0 <- calibrate_xi0(
        theta_true,
        xi1,
        gamma
      )
      
      are_value <- estimate_leading_are(
        theta_true = theta_true,
        xi0 = xi0,
        xi1 = xi1,
        alpha = alpha,
        n_mc = ARE_N_MC,
        seed =
          ARE_SEED +
          1000L * p +
          100L * as.integer(Delta) +
          10L * as.integer(round(alpha * 10))
      )
      
      tibble(
        p = p,
        Delta = Delta,
        alpha = alpha,
        ARE = are_value
      )
    }
  )
  
  write_csv(
    are_obs,
    are_cache
  )
  
  message(
    "Saved ARE cache: ",
    are_cache
  )
}

are_obs <- are_obs %>%
  mutate(
    p = as.integer(p),
    Delta = as.numeric(Delta),
    alpha = as.numeric(alpha),
    ARE = as.numeric(ARE)
  )

sim1000 <- eff %>%
  filter(n == 1000L)

sim500 <- eff %>%
  filter(
    n == 500L,
    near_num(Delta, 2)
  )


# Separate legend labels for theoretical and simulated curves.

fig2 <- ggplot() +
  
  geom_hline(
    yintercept = 1,
    linetype = 3,
    linewidth = 0.45
  ) +
  
  geom_line(
    data = are_obs,
    aes(
      x = alpha,
      y = ARE,
      linetype = "ARE",
      colour = "ARE"
    ),
    linewidth = 0.80
  ) +
  
  geom_line(
    data = sim1000,
    aes(
      x = alpha,
      y = simulated_RE,
      group = 1,
      linetype = "Simulated RE",
      colour = "n = 1000"
    ),
    linewidth = 0.65
  ) +
  
  geom_point(
    data = sim1000,
    aes(
      x = alpha,
      y = simulated_RE,
      shape = "n = 1000",
      colour = "n = 1000"
    ),
    size = 2.2
  ) +
  
  geom_line(
    data = sim500,
    aes(
      x = alpha,
      y = simulated_RE,
      group = 1,
      linetype = "Simulated RE",
      colour = "n = 500"
    ),
    linewidth = 0.55
  ) +
  
  geom_point(
    data = sim500,
    aes(
      x = alpha,
      y = simulated_RE,
      shape = "n = 500",
      colour = "n = 500"
    ),
    size = 2.3,
    stroke = 0.85
  ) +
  
  facet_grid(
    rows = vars(p),
    cols = vars(Delta),
    labeller = labeller(
      p = function(x) {
        paste0("p = ", x)
      },
      Delta = function(x) {
        paste0("\u0394 = ", x)
      }
    )
  ) +
  
  scale_x_continuous(
    breaks = c(
      0,
      0.10,
      0.20,
      0.30,
      0.40
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "ARE" = 1,
      "Simulated RE" = 2
    ),
    name = NULL
  ) +
  
  scale_shape_manual(
    values = c(
      "n = 1000" = 16,
      "n = 500" = 1
    ),
    name = NULL
  ) +
  scale_colour_manual(
    values = c(
      "ARE" = COL_NAVY,
      "n = 1000" = COL_BLUE,
      "n = 500" = COL_ORANGE
    ),
    name = NULL
  ) +
  
  labs(
    x = expression(alpha),
    y = "Relative efficiency"
  ) +
  
  theme_paper()

save_plot(
  fig2,
  "sim_efficiency",
  width = 8.2,
  height = 8.0
)


# ============================================================
# 6. FIGURE 3
# FIXED EXPECTED MISSING-LABEL PROPORTION
#
# Gamma = 0.30 and 0.50 panels
# p = 1 and 5 curves
#
# No bootstrap intervals and no ARE.
# ============================================================

fixed_file <- find_one_file(
  "fixed_missingness_summary.csv"
)

fixed <- read_csv(
  fixed_file,
  show_col_types = FALSE
)

fixed <- rename_first_existing(
  fixed,
  "Gamma",
  c("Gamma", "gamma_total", "Gamma_target")
)

fixed <- rename_first_existing(
  fixed,
  "Delta",
  c("Delta", "delta", "D")
)

fixed <- rename_first_existing(
  fixed,
  "simulated_RE",
  c(
    "simulated_RE",
    "RE",
    "relative_efficiency"
  )
)

# gamma(alpha) may be stored under gamma or gamma_mar.
if (!"gamma_mar" %in% names(fixed)) {
  if ("gamma" %in% names(fixed)) {
    fixed$gamma_mar <- fixed$gamma
  } else {
    fixed$gamma_mar <- (
      as.numeric(fixed$Gamma) -
        as.numeric(fixed$alpha)
    ) /
      (
        1 -
          as.numeric(fixed$alpha)
      )
  }
}

fixed <- fixed %>%
  mutate(
    Gamma = as.numeric(Gamma),
    p = as.integer(p),
    alpha = as.numeric(alpha),
    gamma_mar = as.numeric(gamma_mar),
    simulated_RE = as.numeric(simulated_RE),
    p_lab = factor(
      p,
      levels = c(1, 5),
      labels = c("p = 1", "p = 5")
    )
  )

fig3 <- ggplot(
  fixed,
  aes(
    x = alpha,
    y = simulated_RE,
    group = p_lab,
    linetype = p_lab,
    shape = p_lab,
    colour = p_lab
  )
) +
  
  geom_hline(
    yintercept = 1,
    linetype = 3,
    linewidth = 0.45
  ) +
  
  geom_line(
    linewidth = 0.70
  ) +
  
  geom_point(
    size = 2.2
  ) +
  
  facet_wrap(
    ~ Gamma,
    nrow = 1,
    labeller = labeller(
      Gamma = function(x) {
        paste0("\u0393 = ", x)
      }
    )
  ) +
  
  scale_x_continuous(
    breaks = c(
      0,
      0.10,
      0.20
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "p = 1" = 1,
      "p = 5" = 2
    ),
    name = NULL
  ) +
  
  scale_shape_manual(
    values = c(
      "p = 1" = 16,
      "p = 5" = 1
    ),
    name = NULL
  ) +
  scale_colour_manual(
    values = c(
      "p = 1" = COL_BLUE,
      "p = 5" = COL_ORANGE
    ),
    name = NULL
  ) +
  
  labs(
    x = expression(alpha),
    y = "Simulated relative efficiency"
  ) +
  
  theme_paper()

save_plot(
  fig3,
  "sim_fixed_missingness",
  width = 7.4,
  height = 3.9
)


# ============================================================
# 7. FIGURE 4
# UNOBSERVED MISSING-LABEL INDICATORS
#
# Panel (a):
#   distribution of alpha-hat
#
# Panel (b):
#   simulated RE and ARE
#
# Boundary frequency is NOT computed or displayed.
# ============================================================

# Unobserved/latent RDS files are stored under results/raw/latent_robustness.
LATENT_RAW_DIR <- file.path(
  RAW_DIR,
  "latent_robustness"
)

if (!dir.exists(LATENT_RAW_DIR)) {
  stop(
    "Unobserved/latent raw directory not found:\n",
    LATENT_RAW_DIR
  )
}

latent_files <- list.files(
  LATENT_RAW_DIR,
  pattern = "^latent_.*\\.rds$",
  full.names = TRUE
)

if (!length(latent_files)) {
  stop(
    "No latent_*.rds files found in:\n",
    LATENT_RAW_DIR
  )
}

message(
  "Found ", length(latent_files),
  " unobserved/latent RDS files in ",
  LATENT_RAW_DIR
)


parse_latent_name <- function(f) {
  
  nm <- basename(f)
  
  m <- str_match(
    nm,
    paste0(
      "^latent_",
      "n([0-9]+)_",
      "D([0-9.]+)_",
      "a([0-9.]+)",
      "\\.rds$"
    )
  )
  
  if (any(is.na(m))) {
    stop("Unexpected latent filename: ", nm)
  }
  
  tibble(
    file = f,
    n = as.integer(m[2]),
    Delta = as.numeric(m[3]),
    alpha = as.numeric(m[4])
  )
}

extract_latent_replications <- function(f) {
  
  meta <- parse_latent_name(f)
  res <- readRDS(f)
  
  if (!is.list(res)) {
    stop("Expected a list of replications in: ", f)
  }
  
  vals <- map_dfr(
    seq_along(res),
    function(i) {
      
      x <- res[[i]]
      
      if (!is.null(x$error)) {
        return(NULL)
      }
      
      conv <- x$latent_convergence %||% NA_integer_
      ahat <- x$latent_alpha %||% NA_real_
      
      cc_excess <- x$cc_excess %||% NA_real_
      latent_excess <- x$latent_excess %||% NA_real_
      
      tibble(
        replication = i,
        convergence = as.integer(conv),
        alpha_hat = as.numeric(ahat),
        cc_excess = as.numeric(cc_excess),
        latent_excess = as.numeric(latent_excess)
      )
    }
  )
  
  if (!nrow(vals)) {
    return(tibble())
  }
  
  bind_cols(
    meta[
      rep(1, nrow(vals)),
      c("n", "Delta", "alpha")
    ],
    vals
  )
}

latent_all <- map_dfr(
  latent_files,
  extract_latent_replications
)

if (!nrow(latent_all)) {
  stop(
    "No usable latent/unobserved replications were found."
  )
}


# ------------------------------------------------------------
# Figure 4(a): alpha-hat distribution
#
# Use converged fits only.
# No boundary count/frequency is calculated.
# ------------------------------------------------------------

latent_alpha <- latent_all %>%
  filter(
    convergence == 0L,
    is.finite(alpha_hat)
  ) %>%
  mutate(
    n_lab = factor(
      n,
      levels = c(500, 1000),
      labels = c("n = 500", "n = 1000")
    ),
    alpha_lab = factor(
      sprintf("%.1f", alpha),
      levels = c(
        "0.0",
        "0.1",
        "0.2",
        "0.3"
      )
    )
  )

truth_alpha <- latent_alpha %>%
  distinct(
    Delta,
    alpha,
    alpha_lab
  )

fig4a <- ggplot(
  latent_alpha,
  aes(
    x = alpha_lab,
    y = alpha_hat,
    fill = n_lab
  )
) +
  
  geom_boxplot(
    position = position_dodge(
      width = 0.78
    ),
    width = 0.66,
    linewidth = 0.40,
    outlier.size = 0.40
  ) +
  
  geom_point(
    data = truth_alpha,
    aes(
      x = alpha_lab,
      y = alpha
    ),
    inherit.aes = FALSE,
    shape = 4,
    size = 2.1,
    stroke = 0.8
  ) +
  
  facet_wrap(
    ~ Delta,
    nrow = 1,
    labeller = labeller(
      Delta = function(x) {
        paste0("\u0394 = ", x)
      }
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "n = 500" = COL_SKY,
      "n = 1000" = COL_ORANGE
    ),
    name = NULL
  ) +
  
  labs(
    x = expression("True " * alpha),
    y = expression(hat(alpha))
  ) +
  
  theme_paper()


# ------------------------------------------------------------
# Figure 4(b): simulated RE
#
# RE is computed directly from the raw replications:
#
#   mean(CC excess error) / mean(unobserved excess error)
#
# using converged latent fits.
# ------------------------------------------------------------

latent_eff <- latent_all %>%
  filter(
    convergence == 0L,
    is.finite(cc_excess),
    is.finite(latent_excess)
  ) %>%
  group_by(
    n,
    Delta,
    alpha
  ) %>%
  summarise(
    simulated_RE =
      mean(cc_excess, na.rm = TRUE) /
      mean(latent_excess, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    n_lab = factor(
      n,
      levels = c(500, 1000),
      labels = c("n = 500", "n = 1000")
    )
  )


# ------------------------------------------------------------
# ARE values for the marginal / unobserved experiment.
#
# These are the values used in the current simulation table.
# If you later recompute the marginal-experiment ARE, replace
# only this small table.
# ------------------------------------------------------------

are_unobs <- tribble(
  ~Delta, ~alpha, ~ARE,
  1, 0.00, 1.367907,
  1, 0.10, 1.063920,
  1, 0.20, 0.866121,
  1, 0.30, 0.717745,
  2, 0.00, 2.427513,
  2, 0.10, 1.747027,
  2, 0.20, 1.404630,
  2, 0.30, 1.123807
)


fig4b <- ggplot() +
  
  geom_hline(
    yintercept = 1,
    linetype = 3,
    linewidth = 0.45
  ) +
  
  geom_line(
    data = are_unobs,
    aes(
      x = alpha,
      y = ARE,
      linetype = "ARE",
      colour = "ARE"
    ),
    linewidth = 0.80
  ) +
  
  geom_line(
    data = latent_eff,
    aes(
      x = alpha,
      y = simulated_RE,
      group = n_lab,
      linetype = n_lab,
      colour = n_lab
    ),
    linewidth = 0.60
  ) +
  
  geom_point(
    data = latent_eff,
    aes(
      x = alpha,
      y = simulated_RE,
      shape = n_lab,
      colour = n_lab
    ),
    size = 2.2
  ) +
  
  facet_wrap(
    ~ Delta,
    nrow = 1,
    labeller = labeller(
      Delta = function(x) {
        paste0("\u0394 = ", x)
      }
    )
  ) +
  
  scale_x_continuous(
    breaks = c(
      0,
      0.10,
      0.20,
      0.30
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "ARE" = 1,
      "n = 500" = 2,
      "n = 1000" = 3
    ),
    name = NULL
  ) +
  
  scale_shape_manual(
    values = c(
      "n = 500" = 1,
      "n = 1000" = 16
    ),
    name = NULL
  ) +
  scale_colour_manual(
    values = c(
      "ARE" = COL_NAVY,
      "n = 500" = COL_SKY,
      "n = 1000" = COL_ORANGE
    ),
    name = NULL
  ) +
  
  labs(
    x = expression(alpha),
    y = "Relative efficiency"
  ) +
  
  theme_paper()


# Combine 4(a) and 4(b).

fig4 <- fig4a + fig4b +
  plot_layout(
    ncol = 2,
    widths = c(1, 1),
    guides = "collect"
  ) +
  plot_annotation(
    tag_levels = "a"
  )

save_plot(
  fig4,
  "sim_latent",
  width = 10.2,
  height = 5.0
)


# ============================================================
# 8. FINAL CHECKS
# ============================================================

message("")
message("==============================================")
message("All figures completed.")
message("Output directory:")
message(FIG_DIR)
message("")
message("Files:")
message("  sim_finite_sample.pdf / .png")
message("  sim_efficiency.pdf / .png")
message("  sim_fixed_missingness.pdf / .png")
message("  sim_latent.pdf / .png")
message("==============================================")