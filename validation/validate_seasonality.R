# =============================================================================
# validate_seasonality.R
#
# Cross-check the "seasonality and the timing of control" explainer against
# malariasimulation (mrc-ide) 2.0.2.
#
# The explainer makes five claims that this script checks:
#   1. The SMC drug-prophylaxis curve in the toy is the malariasimulation
#      Weibull survival exp(-(t/scale)^shape) with the SP-AQ parameters.
#   2. In a highly seasonal setting, childhood clinical cases arrive in a single
#      peak that LAGS the rainfall / carrying-capacity peak.
#   3. SMC given as a block of monthly rounds suppresses childhood cases in a
#      SCALLOPED pattern (cases dip after each round, recover in the gaps).
#   4. The same rounds avert FEWER cases when they sit off the case peak than
#      when they are aligned with it (the timing / fixed-calendar point). Moving
#      the SMC calendar relative to a fixed season is equivalent to the season
#      shifting under a fixed calendar, so this is the climate-variability point.
#   5. At equal per-round coverage, RANDOM rounds (inter_round_rho = 0) reach more
#      children and avert more cases than the SAME children each round
#      (inter_round_rho = 1), because doses concentrated on the same children are
#      partly wasted and the never-covered share stays fully exposed.
#
# Run (from the repository root, so figures/ resolves):
#   & 'C:/Program Files/R-aarch64/R-4.5.2/bin/Rscript' validation/validate_seasonality.R
# =============================================================================

.libPaths('C:/Users/pwinskil/Documents/r_packages_arm64')
suppressMessages({
  library(malariasimulation)
  library(ggplot2)
})
set.seed(1)

year <- 365

# ----------------------------------------------------------------------------
# 1. SMC drug prophylaxis curve (exact, from the package)
#    malariasimulation utils.R: weibull_survival(t, shape, scale) = exp(-(t/scale)^shape)
#    drug_parameters.R:         SP_AQ_params = c(0.9, 0.32, 4.3, 38.1)
#                               = c(efficacy, rel_c, prophylaxis_shape, prophylaxis_scale)
# ----------------------------------------------------------------------------
proph_shape <- SP_AQ_params[3]   # 4.3
proph_scale <- SP_AQ_params[4]   # 38.1
weibull_survival <- function(t, shape, scale) exp(-((t / scale) ^ shape))
e_dose <- function(t) ifelse(t < 0, 0, weibull_survival(t, proph_shape, proph_scale))

t50 <- proph_scale * (log(2)) ^ (1 / proph_shape)              # protection = 0.5
t10 <- proph_scale * (log(10)) ^ (1 / proph_shape)             # protection = 0.1

cat("\n================ SMC (SP-AQ) PROPHYLAXIS CURVE ================\n")
cat(sprintf("Weibull survival exp(-(t/%.1f)^%.1f)\n", proph_scale, proph_shape))
prof_days <- c(0, 10, 20, 28, 35, 42, 50, 60)
print(data.frame(day = prof_days, protection = round(e_dose(prof_days), 3)),
      row.names = FALSE)
cat(sprintf("Half protection at %.1f d; 10%% protection at %.1f d\n", t50, t10))

# ----------------------------------------------------------------------------
# 2. A highly seasonal (Sahelian-type) profile
#    Fourier carrying capacity C(t) = max(floor, g0 + sum_k g_k cos + h_k sin).
# ----------------------------------------------------------------------------
seasonal <- list(
  model_seasonality = TRUE,
  g0 = 0.285, g = c(-0.325, -0.132, 0.104), h = c(-0.419, -0.158, -0.061),
  rainfall_floor = 0.001
)

human_population <- 20000
init_EIR         <- 50          # high baseline so the seasonal peak is pronounced
warmup_yr        <- 6           # reach seasonal equilibrium
analysis_yr      <- warmup_yr + 1
sim_length       <- analysis_yr * year

# SMC eligible age group: children 3 to 59 months (WHO). We treat and count cases in
# this same cohort, so the treated and monitored children match.
child_min <- round(0.25 * year)   # ~3 months (91 days)
child_max <- 5 * year             # 60 months
inc_col   <- paste0("n_inc_clinical_", child_min, "_", child_max)
n_col     <- paste0("n_age_",          child_min, "_", child_max)

base_params <- function() {
  p <- get_parameters(c(seasonal, list(
    human_population = human_population,
    clinical_incidence_rendering_min_ages = c(child_min, 0),
    clinical_incidence_rendering_max_ages = c(child_max, 200 * year)
  )))
  p <- set_species(p, list(gamb_params), 1)
  p <- set_equilibrium(p, init_EIR = init_EIR)
  p
}

peak_day <- peak_season_offset(base_params())   # day of the carrying-capacity peak
cat("\n================ SEASONAL PROFILE ================\n")
cat(sprintf("Rainfall / carrying-capacity peak (peak_season_offset): day %d\n", peak_day))

# ----------------------------------------------------------------------------
# 3. SMC schedules: a block of 4 monthly rounds, placed in the analysis year.
#    Aligned: block centred on the seasonal peak.
#    Shifted: same block moved 8 weeks later (off the peak) -> the timing test.
# ----------------------------------------------------------------------------
n_rounds     <- 4
spacing_days <- 28
smc_coverage <- 0.9
analysis_t0  <- (analysis_yr - 1) * year          # first timestep of the analysis year

block_start <- function(first_doy) analysis_t0 + first_doy + (0:(n_rounds - 1)) * spacing_days

# centre the 4-round block (~84 days wide) on the peak
aligned_first <- peak_day - round(((n_rounds - 1) * spacing_days) / 2)
aligned_rounds <- block_start(aligned_first)
shifted_rounds <- block_start(aligned_first + 8 * 7)   # 8 weeks late

add_smc <- function(p, rounds, cov = smc_coverage) {
  p <- set_drugs(p, list(SP_AQ_params))
  set_smc(
    p, drug = 1, timesteps = rounds,
    coverages = rep(cov,       n_rounds),
    min_ages  = rep(child_min, n_rounds),
    max_ages  = rep(child_max, n_rounds)
  )
}

# ----------------------------------------------------------------------------
# Run: baseline (no SMC), SMC aligned, SMC shifted
# ----------------------------------------------------------------------------
message("Running 3 scenarios (pop=", human_population, ", ", analysis_yr,
        " yr each). This takes a few minutes...")

message("  baseline (no SMC) ...")
out_base <- run_simulation(timesteps = sim_length, parameters = base_params())
message("  SMC aligned to peak ...")
out_algn <- run_simulation(timesteps = sim_length, parameters = add_smc(base_params(), aligned_rounds))
message("  SMC shifted 8 weeks late ...")
out_shft <- run_simulation(timesteps = sim_length, parameters = add_smc(base_params(), shifted_rounds))

# ----------------------------------------------------------------------------
# Extract the analysis-year daily 3-59 month case curve (per person per day), smoothed.
# ----------------------------------------------------------------------------
roll_mean <- function(x, k = 7) {
  n <- length(x); out <- numeric(n); half <- (k - 1) %/% 2
  for (i in seq_len(n)) { lo <- max(1, i - half); hi <- min(n, i + half); out[i] <- mean(x[lo:hi]) }
  out
}

year_curve <- function(out) {
  out$doy  <- ((out$timestep - 1) %% year) + 1
  out$yr   <- ((out$timestep - 1) %/% year) + 1
  d <- out[out$yr == analysis_yr, ]
  d <- d[order(d$doy), ]
  data.frame(doy = d$doy,
             inc = roll_mean(d[[inc_col]] / d[[n_col]]))
}

cb <- year_curve(out_base)
ca <- year_curve(out_algn)
cs <- year_curve(out_shft)

# case peak (baseline) and EIR peak in the analysis year
case_peak_doy <- cb$doy[which.max(cb$inc)]
ob <- out_base; ob$doy <- ((ob$timestep - 1) %% year) + 1
ob$yr <- ((ob$timestep - 1) %/% year) + 1
ob <- ob[ob$yr == analysis_yr, ]
eir_by_doy   <- tapply(ob$EIR_gamb, ob$doy, mean)
eir_peak_doy <- as.integer(names(eir_by_doy)[which.max(eir_by_doy)])

# annual 3-59 month cases and % averted
annual <- function(out) {
  out$yr <- ((out$timestep - 1) %/% year) + 1
  d <- out[out$yr == analysis_yr, ]
  sum(d[[inc_col]]) / mean(d[[n_col]])
}
a_base <- annual(out_base); a_algn <- annual(out_algn); a_shft <- annual(out_shft)

cat("\n================ SEASONAL CASE TIMING (analysis year) ================\n")
cat(sprintf("Rainfall peak day        : %d\n", peak_day))
cat(sprintf("EIR peak day             : %d  (lag vs rainfall = %+d d)\n",
            eir_peak_doy, eir_peak_doy - peak_day))
cat(sprintf("Case (3-59mo) peak day   : %d  (lag vs rainfall = %+d d)\n",
            case_peak_doy, case_peak_doy - peak_day))
cat(sprintf("SMC aligned rounds (doy) : %s\n",
            paste(aligned_rounds - analysis_t0, collapse = ", ")))
cat(sprintf("SMC shifted rounds (doy) : %s\n",
            paste(shifted_rounds - analysis_t0, collapse = ", ")))

cat("\n================ ANNUAL 3-59mo CASES AND % AVERTED ================\n")
print(data.frame(
  scenario = c("baseline", "SMC aligned", "SMC shifted +8wk"),
  cases_3_59mo = round(c(a_base, a_algn, a_shft), 3),
  pct_averted = round(100 * c(0, 1 - a_algn / a_base, 1 - a_shft / a_base), 1)
), row.names = FALSE)
cat("\nExpect: EIR and case peaks LAG the rainfall peak;\n",
    "        aligned SMC averts MORE than shifted SMC (timing / fixed-calendar point).\n", sep = "")

# ----------------------------------------------------------------------------
# 5. Coverage correlation between rounds (inter_round_rho on 'smc').
#    Same aligned rounds and coverage, but rho = 0 (random, independent rounds)
#    vs rho = 1 (same children each round). Expect random to reach more children
#    and avert more cases. Coverage held at a moderate level so the never-covered
#    share under rho = 1 is visible.
# ----------------------------------------------------------------------------
corr_cov <- 0.6
run_corr <- function(rho) {
  p  <- add_smc(base_params(), aligned_rounds, cov = corr_cov)
  cp <- get_correlation_parameters(p)
  cp$inter_round_rho("smc", rho)
  run_simulation(timesteps = sim_length, parameters = p, correlations = cp)
}
message("  SMC coverage ", corr_cov, ", rho = 0 (random rounds) ...")
out_rho0 <- run_corr(0)
message("  SMC coverage ", corr_cov, ", rho = 1 (same children) ...")
out_rho1 <- run_corr(1)

a_rho0 <- annual(out_rho0); a_rho1 <- annual(out_rho1)

cat("\n================ COVERAGE CORRELATION BETWEEN ROUNDS ================\n")
cat(sprintf("Aligned rounds, coverage %.0f%% per round.\n", corr_cov * 100))
print(data.frame(
  inter_round_rho = c("0 (random)", "1 (same children)"),
  cases_3_59mo = round(c(a_rho0, a_rho1), 3),
  pct_averted  = round(100 * c(1 - a_rho0 / a_base, 1 - a_rho1 / a_base), 1)
), row.names = FALSE)
cat("\nExpect: rho = 0 (random) averts MORE than rho = 1 (same children) at equal coverage.\n")

# ----------------------------------------------------------------------------
# 6. Overlay figure: malariasimulation vs the toy model
#    Panel A: SMC prophylaxis curve (toy == package, so a single curve).
#    Panel B: analysis-year 3-59mo cases, baseline vs SMC-aligned,
#             malariasimulation (solid) with the toy model (dashed) overlaid.
# ----------------------------------------------------------------------------
# ---- toy model, identical functions to seasonality.html ---------------------
# rise = fraction of the year the trough-to-peak rise spans (speed of take-off slider:
# 0.67 slow, 0.5 medium, 0.33 fast). rise = 0.5 is the symmetric ((1+cos)/2)^p curve.
toy_seasonal <- function(t, S, peak, rise = 0.5) {
  p <- 1 + 8 * S + 16 * S^4; floor <- (1 - S) ^ 1.5
  R <- rise * year; F <- year - R
  d <- (t - peak) %% year                                   # days since peak, 0..year
  phi <- ifelse(d <= F, pi * d / F, pi * (d - year) / R)    # falling 0..pi, then rising -pi..0
  s <- ((1 + cos(phi)) / 2) ^ p
  floor + (1 - floor) * s
}
# DRUG_EFF matches seasonality.html: a covered child's protection is 0.9 * e(t). The
# malariasimulation overlay run uses ~90% coverage with independent rounds, so nearly
# every child receives nearly every round and the population ceiling is ~the efficacy.
DRUG_EFF <- 0.9
toy_protection <- function(t, rounds_doy) {
  best <- 0
  # days since dose measured around the year (matches seasonality.html): a late round
  # carries protection into the start of the season, as the periodic case curve does.
  for (ti in rounds_doy) { dt <- (t - ti) %% year; e <- e_dose(dt); if (e > best) best <- e }
  DRUG_EFF * best
}

doy_grid    <- 1:year
toy_S       <- 0.85                       # highly seasonal, matching the profile used
toy_peak    <- case_peak_doy              # align the toy season to the modelled case peak
toy_rounds  <- aligned_rounds - analysis_t0
toy_c0 <- sapply(doy_grid, toy_seasonal,   S = toy_S, peak = toy_peak)
toy_P  <- sapply(doy_grid, toy_protection, rounds_doy = toy_rounds)
toy_c  <- toy_c0 * (1 - toy_P)

# normalise every curve to its own baseline peak so shapes are comparable
norm <- function(v, ref) v / max(ref)
plot_df <- rbind(
  data.frame(doy = cb$doy, inc = norm(cb$inc, cb$inc), src = "malariasimulation", line = "no SMC"),
  data.frame(doy = ca$doy, inc = norm(ca$inc, cb$inc), src = "malariasimulation", line = "SMC aligned"),
  data.frame(doy = doy_grid, inc = norm(toy_c0, toy_c0), src = "toy model", line = "no SMC"),
  data.frame(doy = doy_grid, inc = norm(toy_c,  toy_c0), src = "toy model", line = "SMC aligned")
)

round_df <- data.frame(doy = toy_rounds)

g <- ggplot(plot_df, aes(doy, inc, colour = line, linetype = src)) +
  geom_vline(data = round_df, aes(xintercept = doy),
             colour = "#e8804b", linetype = "dotted", linewidth = 0.4) +
  geom_line(linewidth = 0.85) +
  scale_colour_manual(values = c("no SMC" = "#3d6fb4", "SMC aligned" = "#7b5bd6"),
                      name = NULL) +
  scale_linetype_manual(values = c("malariasimulation" = "solid", "toy model" = "dashed"),
                        name = NULL) +
  labs(
    title = "Seasonality and SMC timing: toy model vs malariasimulation",
    subtitle = paste0("Highly seasonal profile, 4 monthly SMC rounds (dotted) centred on the ",
                      "case peak.\nChildhood (3-59 month) cases through the year, each normalised to ",
                      "its own no-SMC peak."),
    x = "Day of year", y = "3-59mo cases (x no-SMC peak)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

dir.create("figures", showWarnings = FALSE)
png_path <- file.path(getwd(), "figures", "seasonality_toy_vs_malariasim.png")
ggsave(png_path, g, width = 10, height = 5, dpi = 120)
cat("\nSaved overlay plot to: ", png_path, "\n", sep = "")
