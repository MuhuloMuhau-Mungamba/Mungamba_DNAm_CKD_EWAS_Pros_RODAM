############################################
# Project: CKD EWAS (RODAM-Pros)
# Analysis: Power analysis
# Purpose: Reviewer 3 – Power calculation
# Author: Muhulo Muhau Mungamba
# Date: 2026-04-16

############################################
# -------- POWER CALCULATION --------
############################################

# Install if needed
install.packages("pwr")

library(pwr)

# Your sample size (approximate)
n_cases <- 92
n_controls <- 182
n_total <- n_cases + n_controls

# Assume small EWAS effect sizes (realistic)
effect_sizes <- c(0.02, 0.05, 0.1)  # small, moderate, larger

# Calculate power for each effect size
power_results <- lapply(effect_sizes, function(f2) {
  pwr.f2.test(
    u = 1,                 # predictor (methylation)
    v = n_total - 2,       # residual df
    f2 = f2,
    sig.level = 1e-7       # EWAS threshold (stringent)
  )
})

# Format output
power_table <- data.frame(
  Effect_Size_f2 = effect_sizes,
  Power = sapply(power_results, function(x) x$power)
)

print(power_table)
