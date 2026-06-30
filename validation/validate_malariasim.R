# =============================================================================
# validate_malariasim.R
#
# Verify a "toy" malaria rebound model against malariasimulation (mrc-ide) 2.0.2.
#
# Scenario: equilibrium warmup -> IRS sprayed annually for 4 years -> STOP ->
#           long follow-up. Repeated at low / moderate / high baseline EIR.
#
# Question: during IRS, does clinical incidence drop and stay BELOW the
#           pre-intervention level? After IRS stops, does it REBOUND, and does
#           it overshoot the pre-IRS baseline? How does rebound depend on EIR?
#
# Run:
#   & 'C:/Program Files/R-aarch64/R-4.5.2/bin/Rscript' validate_malariasim.R
# =============================================================================

.libPaths('C:/Users/pwinskil/Documents/r_packages_arm64')
suppressMessages({
  library(malariasimulation)
  library(ggplot2)
})
set.seed(1)

year <- 365

# ----------------------------------------------------------------------------
# Timeline (in years -> timesteps)
# ----------------------------------------------------------------------------
warmup_yr     <- 6     # reach equilibrium before any intervention
irs_years     <- 4     # annual IRS rounds
followup_yr   <- 16    # follow-up after IRS stops
total_yr      <- warmup_yr + irs_years + followup_yr   # 26 years
sim_length    <- total_yr * year

irs_start_yr  <- warmup_yr
irs_stop_yr   <- warmup_yr + irs_years
# annual spray rounds at the start of each IRS year
spray_timesteps <- (warmup_yr + 0:(irs_years - 1)) * year + 1

human_population <- 10000
baseline_EIRs    <- c(5, 30, 100)   # low / moderate / high

# ----------------------------------------------------------------------------
# Realistic IRS product efficacy (pirimiphos-methyl / Actellic-300CS-like),
# single Anopheles gambiae species. Parameters follow the IRS vector-control
# model of Sherrard-Smith et al. (2018) Nature Communications,
# 10.1038/s41467-018-07357-w (the S.I. referenced by ?set_spraying):
#   ls_* mortality, ks_* feeding-success, ms_* deterrence (logistic decay).
# Gives a strong initial kill that decays over ~6-12 months.
# ----------------------------------------------------------------------------
n_rounds <- length(spray_timesteps)
mat <- function(v) matrix(v, nrow = n_rounds, ncol = 1)   # 1 species column

irs_pars <- list(
  ls_theta = mat( 2.025),  ls_gamma = mat(-0.009),   # mortality
  ks_theta = mat(-2.222),  ks_gamma = mat( 0.008),   # feeding success
  ms_theta = mat(-1.232),  ms_gamma = mat(-0.009)    # deterrence
)
irs_coverage <- 0.8

# ----------------------------------------------------------------------------
# Helper: run one baseline-EIR scenario, return tidy time series of clinical
# incidence per person per year (all ages and 0-5y), aggregated to years.
# ----------------------------------------------------------------------------
run_scenario <- function(init_EIR) {
  p <- get_parameters(list(
    human_population = human_population,
    # render clinical incidence for 0-5y and all-ages (0-200y)
    clinical_incidence_rendering_min_ages = c(0,        0) * year,
    clinical_incidence_rendering_max_ages = c(5 * year, 200 * year)
  ))
  p <- set_equilibrium(p, init_EIR = init_EIR)

  p <- set_spraying(
    p,
    timesteps = spray_timesteps,
    coverages = rep(irs_coverage, n_rounds),
    ls_theta  = irs_pars$ls_theta, ls_gamma = irs_pars$ls_gamma,
    ks_theta  = irs_pars$ks_theta, ks_gamma = irs_pars$ks_gamma,
    ms_theta  = irs_pars$ms_theta, ms_gamma = irs_pars$ms_gamma
  )

  out <- run_simulation(timesteps = sim_length, parameters = p)
  out$year <- (out$timestep - 1) %/% year + 1   # 1-based year index

  # Aggregate to calendar years: incidence per person per year =
  #   (sum of new clinical cases in year) / (mean population in age band)
  agg <- function(inc_col, n_col) {
    cases <- tapply(out[[inc_col]], out$year, sum)
    popn  <- tapply(out[[n_col]],   out$year, mean)
    as.numeric(cases / popn)
  }

  yrs <- sort(unique(out$year))
  # Annual mean of the model's EIR output (EIR_gamb). Only used below as a
  # RATIO (during-IRS / pre-IRS) to estimate the toy's control strength S, so
  # the absolute scale is irrelevant; we keep the model's native units.
  eir_yr <- as.numeric(tapply(out$EIR_gamb, out$year, mean))
  data.frame(
    EIR      = init_EIR,
    year     = yrs,
    inc_all  = agg("n_inc_clinical_0_73000", "n_age_0_73000"),
    inc_0_5  = agg("n_inc_clinical_0_1825",  "n_age_0_1825"),
    eir_real = eir_yr
  )
}

# ----------------------------------------------------------------------------
# Run all scenarios
# ----------------------------------------------------------------------------
message("Running ", length(baseline_EIRs), " scenarios (pop=",
        human_population, ", ", total_yr, " yr each). This takes a few minutes...")

results <- do.call(rbind, lapply(baseline_EIRs, function(e) {
  message("  init_EIR = ", e, " ...")
  run_scenario(e)
}))
results$EIR_lab <- factor(paste0("EIR = ", results$EIR),
                          levels = paste0("EIR = ", baseline_EIRs))

# ----------------------------------------------------------------------------
# Numeric summary: pre-IRS, min-during-IRS, peak-after-stop, rebound multiple
# Use whole years. Pre-IRS = mean of last 2 warmup years (equilibrium).
# ----------------------------------------------------------------------------
summ <- do.call(rbind, lapply(baseline_EIRs, function(e) {
  d <- results[results$EIR == e, ]
  pre_idx   <- d$year %in% c(warmup_yr - 1, warmup_yr)          # equilibrium
  irs_idx   <- d$year %in% (irs_start_yr + 1):irs_stop_yr        # during IRS
  post_idx  <- d$year >  irs_stop_yr                            # after stop

  pre_all   <- mean(d$inc_all[pre_idx])
  min_all   <- min(d$inc_all[irs_idx])
  peak_all  <- max(d$inc_all[post_idx])
  peak_yr   <- d$year[post_idx][which.max(d$inc_all[post_idx])] - irs_stop_yr

  pre_05    <- mean(d$inc_0_5[pre_idx])
  min_05    <- min(d$inc_0_5[irs_idx])
  peak_05   <- max(d$inc_0_5[post_idx])

  data.frame(
    EIR = e,
    pre_IRS_all      = round(pre_all, 3),
    min_during_all   = round(min_all, 3),
    peak_after_all   = round(peak_all, 3),
    rebound_mult_all = round(peak_all / pre_all, 3),
    yrs_to_peak      = peak_yr,
    pre_IRS_0_5      = round(pre_05, 3),
    min_during_0_5   = round(min_05, 3),
    peak_after_0_5   = round(peak_05, 3),
    rebound_mult_0_5 = round(peak_05 / pre_05, 3)
  )
}))

cat("\n================ CLINICAL INCIDENCE SUMMARY (cases / person / yr) ================\n")
cat("(IRS sprayed yrs ", irs_start_yr, "-", irs_stop_yr,
    "; pre-IRS = equilibrium; min during IRS; peak after stop)\n\n", sep = "")
print(summ, row.names = FALSE)
cat("\nrebound_mult = peak_after_stop / pre_IRS  (>1 means overshoot of baseline)\n")
cat("yrs_to_peak  = years after IRS stop until the post-IRS peak\n")

# ----------------------------------------------------------------------------
# Plot: all-ages clinical incidence over time, IRS markers + pre-IRS baseline
# ----------------------------------------------------------------------------
pre_lines <- data.frame(EIR = summ$EIR, pre = summ$pre_IRS_all)
pre_lines$EIR_lab <- factor(paste0("EIR = ", pre_lines$EIR),
                            levels = levels(results$EIR_lab))

g <- ggplot(results, aes(year, inc_all)) +
  annotate("rect", xmin = irs_start_yr, xmax = irs_stop_yr,
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.3) +
  geom_hline(data = pre_lines, aes(yintercept = pre),
             linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = c(irs_start_yr, irs_stop_yr),
             linetype = "dotted", colour = "steelblue") +
  geom_line(aes(colour = "All ages"), linewidth = 0.9) +
  geom_line(aes(y = inc_0_5, colour = "Children 0-5y"),
            linewidth = 0.7, alpha = 0.8) +
  facet_wrap(~EIR_lab, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c("All ages" = "firebrick",
                                 "Children 0-5y" = "darkorange"),
                      name = NULL) +
  labs(
    title = "malariasimulation 2.0.2: clinical incidence under IRS then withdrawal",
    subtitle = paste0("IRS (Actellic-like) sprayed annually yr ", irs_start_yr,
                      "-", irs_stop_yr, " at ", irs_coverage*100,
                      "% coverage, then stopped.\nDotted = IRS start/stop; ",
                      "dashed grey = pre-IRS (equilibrium) incidence."),
    x = "Year", y = "Clinical incidence (cases per person per year)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

# (standalone malariasim-only figure is not saved; the overlay below is the single
#  diagnostic we keep. `g` remains available if a standalone plot is ever wanted.)

# =============================================================================
# OVERLAY DIAGNOSTIC: toy rebound model vs malariasimulation
#
# Plot CLINICAL INCIDENCE RELATIVE TO THE PRE-IRS LEVEL (pre-IRS = 1.0) against
# YEARS SINCE IRS START, one panel per baseline EIR. malariasim curve solid;
# toy-model curve dashed overlaid. Lets us see whether the toy reproduces the
# suppression-during-IRS and rebound-after-stop SHAPE.
# =============================================================================

# ---- toy model (protection model from rebounds.html; incidence derived here
#      only, as a SHAPE check against malariasimulation) ------------------------
K              <- 28
g_imm          <- function(x) x^2 / (x^2 + K^2)
rate_of        <- function(Tr) 0.08 + 0.20 * (Tr / 100)   # immunity relaxation rate
loss_half_life <- function(v) 3.0 * (0.17 / 3.0) ^ (v / 100)

# The web lesson plots total PROTECTION (immunity + intervention); here we derive
# a clinical-incidence proxy rr = (1-J)(1-c*I)/(1-c*B) from the same immunity stock
# purely to compare SHAPE (suppress -> overshoot, larger at higher transmission)
# against malariasimulation. Magnitudes are illustrative, not fitted.
c_prot <- 0.78
toy_rr <- function(Tr, S, D, lossV, t_grid) {
  hl      <- loss_half_life(lossV)
  k_decay <- log(2) / hl
  B       <- g_imm(Tr)
  r       <- rate_of(Tr)
  dt      <- 0.02
  ts <- seq(min(t_grid), max(t_grid), by = dt)
  I  <- numeric(length(ts)); J <- numeric(length(ts)); I[1] <- B
  for (i in seq_along(ts)) {
    t  <- ts[i]
    Jt <- if (t >= 0 && t < D) S
          else if (t >= D)     S * exp(-k_decay * (t - D))
          else                 0
    J[i] <- Jt
    Ieq  <- g_imm(Tr * (1 - Jt))
    if (i > 1) I[i] <- I[i - 1] + r * (Ieq - I[i - 1]) * dt
  }
  I  <- pmax(0, I)
  rr <- (1 - J) * (1 - c_prot * I) / max(1e-6, 1 - c_prot * B)
  approx(ts, rr, xout = t_grid)$y    # sample onto the requested grid
}

# IRS protection decays over months -> IRS-like loss slider.
# lossV = 85 gives half-life loss_half_life(85) ~ 0.26 yr (~3 months).
toy_lossV <- 85
cat(sprintf("\nToy overlay: IRS-like loss slider lossV = %d (half-life %.2f yr).\n",
            toy_lossV, loss_half_life(toy_lossV)))

# ---- build overlay data: one toy curve per baseline EIR --------------------
# Time axis = years since IRS start (warmup excluded -> negative pre-IRS years).
D_toy <- irs_years   # toy control duration = number of IRS years sprayed

overlay <- do.call(rbind, lapply(baseline_EIRs, function(e) {
  d <- results[results$EIR == e, ]
  # years since IRS start (whole-year index minus the warmup)
  d$t <- d$year - irs_start_yr

  # pre-IRS reference incidence = equilibrium (last 2 warmup years)
  pre_idx <- d$year %in% c(warmup_yr - 1, warmup_yr)
  pre_inc <- mean(d$inc_all[pre_idx])

  # realised transmission reduction during IRS, from THIS run's mean EIR:
  #   S = 1 - (mean EIR while IRS active) / (pre-IRS EIR)
  pre_eir_idx <- pre_idx
  irs_eir_idx <- d$year %in% (irs_start_yr + 1):irs_stop_yr
  pre_eir     <- mean(d$eir_real[pre_eir_idx])
  irs_eir     <- mean(d$eir_real[irs_eir_idx])
  S_run       <- max(0, min(0.95, 1 - irs_eir / pre_eir))

  # toy curve on the same years-since-start grid; Tr = baseline EIR
  t_grid  <- d$t
  rr_toy  <- toy_rr(Tr = e, S = S_run, D = D_toy, lossV = toy_lossV,
                    t_grid = t_grid)

  cat(sprintf("  EIR %3d: during-IRS / pre-IRS EIR ratio = %.3f -> S = %.2f\n",
              e, irs_eir / pre_eir, S_run))

  rbind(
    data.frame(EIR = e, t = d$t, rel = d$inc_all / pre_inc,
               src = "malariasimulation"),
    data.frame(EIR = e, t = t_grid, rel = rr_toy, src = "toy model")
  )
}))
overlay$EIR_lab <- factor(paste0("baseline EIR = ", overlay$EIR),
                          levels = paste0("baseline EIR = ", baseline_EIRs))

# trim to the plotted window: a couple of years before IRS through follow-up
overlay <- overlay[overlay$t >= -2, ]

g2 <- ggplot(overlay, aes(t, rel, colour = src, linetype = src)) +
  geom_hline(yintercept = 1, colour = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(0, irs_years), linetype = "dotted",
             colour = "steelblue") +
  geom_line(linewidth = 0.9) +
  facet_wrap(~EIR_lab, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = c("malariasimulation" = "firebrick",
                                 "toy model" = "#2f8f7e"), name = NULL) +
  scale_linetype_manual(values = c("malariasimulation" = "solid",
                                   "toy model" = "dashed"), name = NULL) +
  labs(
    title = "Toy rebound model vs malariasimulation: clinical incidence relative to pre-IRS",
    subtitle = paste0("Dotted lines = IRS start (yr 0) and stop (yr ", irs_years,
                      "); horizontal line = pre-IRS level (= 1.0).\n",
                      "Toy: Tr = baseline EIR, D = ", D_toy,
                      " yr, S from realised EIR drop, IRS-like loss (lossV = ",
                      toy_lossV, ")."),
    x = "Years since IRS start", y = "Clinical incidence (x pre-IRS level)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

dir.create("figures", showWarnings = FALSE)
png_path2 <- file.path(getwd(), "figures", "toy_vs_malariasim_overlay.png")
ggsave(png_path2, g2, width = 12, height = 5, dpi = 120)
cat("\nSaved overlay plot to: ", png_path2, "\n", sep = "")
