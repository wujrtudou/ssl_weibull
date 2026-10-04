`%||%` <- function(x, y) if (is.null(x)) y else x

softplus <- function(x) {
  out <- numeric(length(x))
  hi <- x > 30
  lo <- x < -30
  mid <- !(hi | lo)
  out[hi] <- x[hi] + log1p(exp(-x[hi]))
  out[lo] <- exp(x[lo])
  out[mid] <- log1p(exp(x[mid]))
  out
}

log_sigmoid <- function(x) -softplus(-x)
log1m_sigmoid <- function(x) -softplus(x)

logsumexp2 <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

safe_log <- function(x, floor = 1e-300) log(pmax(x, floor))

rmse <- function(x, truth) sqrt(mean((x - truth)^2, na.rm = TRUE))

make_id <- function(...) {
  x <- paste(..., sep = "_")
  gsub("[^A-Za-z0-9_.-]", "", x)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

project_root <- function() {
  root <- getOption("sslmm.root", NULL)
  if (!is.null(root)) return(normalizePath(root, mustWork = FALSE))
  normalizePath(getwd(), mustWork = FALSE)
}

source_sslmm <- function(root = project_root()) {
  fs <- sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
  for (f in fs) source(f, local = .GlobalEnv)
  invisible(fs)
}

parallel_lapply <- function(X, FUN, n_cores = 1L, ...) {
  n_cores <- max(1L, as.integer(n_cores))
  if (n_cores <= 1L || .Platform$OS.type == "windows") {
    return(lapply(X, FUN, ...))
  }
  parallel::mclapply(X, FUN, ..., mc.cores = n_cores, mc.preschedule = FALSE)
}

numeric_jacobian <- function(fun, x, rel_step = 1e-5) {
  f0 <- fun(x)
  J <- matrix(NA_real_, nrow = length(f0), ncol = length(x))
  for (j in seq_along(x)) {
    h <- rel_step * max(1, abs(x[j]))
    xp <- xm <- x
    xp[j] <- xp[j] + h
    xm[j] <- xm[j] - h
    J[, j] <- (fun(xp) - fun(xm)) / (2 * h)
  }
  dimnames(J) <- list(names(f0), names(x))
  J
}

safe_inverse <- function(H, tol = 1e-10) {
  if (any(!is.finite(H))) return(NULL)
  H <- (H + t(H)) / 2
  ee <- try(eigen(H, symmetric = TRUE), silent = TRUE)
  if (inherits(ee, "try-error") || min(ee$values) <= tol) return(NULL)
  solve(H)
}

run_multistart_nlminb <- function(starts, objective, lower, upper, control = list(), ...) {
  fits <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    st <- pmin(pmax(starts[[i]], lower), upper)
    fits[[i]] <- try(
      nlminb(start = st, objective = objective, lower = lower, upper = upper,
             control = control, ...),
      silent = TRUE
    )
  }
  ok <- vapply(fits, function(z) !inherits(z, "try-error") && is.finite(z$objective), logical(1))
  if (!any(ok)) stop("All nlminb starts failed.")
  cand <- fits[ok]
  conv <- vapply(cand, function(z) isTRUE(z$convergence == 0L), logical(1))
  pool <- if (any(conv)) cand[conv] else cand
  best <- pool[[which.min(vapply(pool, `[[`, numeric(1), "objective"))]]
  all_summary <- do.call(rbind, lapply(cand, function(z) {
    data.frame(objective = z$objective, convergence = z$convergence,
               iterations = z$iterations %||% NA_integer_, stringsAsFactors = FALSE)
  }))
  list(best = best, all = cand, summary = all_summary)
}

checkpoint_mc <- function(cfg, B, worker, out_file, n_cores = 1L,
                          batch_size = 25L, seed_base = 1L, overwrite = FALSE, ...) {
  ensure_dir(dirname(out_file))
  if (overwrite && file.exists(out_file)) file.remove(out_file)
  ans <- if (file.exists(out_file)) readRDS(out_file) else vector("list", 0L)
  done <- length(ans)
  if (done >= B) return(ans[seq_len(B)])

  n_cores <- max(1L, as.integer(n_cores))
  dots <- list(...)

  # Windows does not support fork-based mclapply().  Use a persistent
  # PSOCK cluster for the whole checkpoint run instead, so workers are
  # created once per configuration rather than once per batch.
  cl <- NULL
  if (.Platform$OS.type == "windows" && n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    root <- project_root()
    parallel::clusterExport(cl, "root", envir = environment())
    parallel::clusterEvalQ(cl, {
      options(sslmm.root = root)
      source(file.path(root, "config.R"))
      source(file.path(root, "R", "00_utils.R"))
      source_sslmm(root)
      NULL
    })
  }

  while (done < B) {
    idx <- seq.int(done + 1L, min(B, done + batch_size))

    if (!is.null(cl)) {
      res <- parallel::parLapply(
        cl, idx,
        function(b, cfg, worker, seed_base, dots) {
          set.seed(as.integer(seed_base + b))
          tryCatch(
            do.call(worker, c(list(cfg = cfg, rep_id = b), dots)),
            error = function(e) list(rep_id = b, error = conditionMessage(e))
          )
        },
        cfg = cfg, worker = worker, seed_base = seed_base, dots = dots
      )
    } else {
      one <- function(b) {
        set.seed(as.integer(seed_base + b))
        tryCatch(
          do.call(worker, c(list(cfg = cfg, rep_id = b), dots)),
          error = function(e) list(rep_id = b, error = conditionMessage(e))
        )
      }
      res <- parallel_lapply(idx, one, n_cores = n_cores)
    }

    ans <- c(ans, res)
    saveRDS(ans, out_file)
    done <- length(ans)
    message(sprintf("%s: %d/%d replications saved", basename(out_file), done, B))
  }
  ans
}
