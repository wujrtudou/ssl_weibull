
root <- "C:/Users/uqjwu15/Desktop/ssl_mm_simulation_windows_8core"

setwd(root)

source(file.path(root, "scripts", "01_parameter_recovery.R"))
source(file.path(root, "scripts", "02_efficiency.R"))
source(file.path(root, "scripts", "03_fixed_missingness.R"))
source(file.path(root, "scripts", "04_latent_robustness.R"))


SIM_CONTROL$n_cores <- 8L

cat("Project:", root, "\n")
cat("Workers:", SIM_CONTROL$n_cores, "\n")
cat("B per configuration:", SIM_CONTROL$B, "\n")

start_time <- Sys.time()
cat("Start time:", format(start_time), "\n\n")

# 1. Parameter recovery
cat("\n========== 1/4 Parameter Recovery ==========\n")
run_parameter_recovery()

# 2. Efficiency
cat("\n========== 2/4 Efficiency ==========\n")
run_efficiency()

# 3. Fixed missingness
cat("\n========== 3/4 Fixed Missingness ==========\n")
run_fixed_missingness()

# 4. Latent robustness
cat("\n========== 4/4 Latent Robustness ==========\n")
run_latent_robustness()

end_time <- Sys.time()

cat("\n============================================\n")
cat("ALL SIMULATIONS FINISHED\n")
cat("Start :", format(start_time), "\n")
cat("Finish:", format(end_time), "\n")
cat(
  "Total hours:",
  round(as.numeric(difftime(end_time, start_time, units = "hours")), 2),
  "\n"
)
cat("============================================\n")