# =============================================================================
# validate_nmf.R
#
# Static R/ggplot reproduction of the interactive toy model in nmf.html, as an
# independent cross-check of the JavaScript. It re-implements the same quantities
# through a single year (day-of-year t = 0..364):
#
#   g(t)  seasonal shape (single peak fixed mid-year, peak = 1)
#   C(t)  true clinical malaria       = c(u) * g(t)          [level set by transmission]
#   P(t)  prevalence of infection     = v(u) * ghat(t)       [level set by transmission]
#            ghat = broadened/lagged/floored low-pass of g (relaxation, timescale D)
#   M(t)  true malaria cases recorded = treat * C(t)
#   X(t)  falsely-attributed cases    = treat * kappa * n * P(t)   [coincident positives]
#   R(t)  routine-recorded total      = M(t) + X(t)
#
# The point of the model: non-malaria fevers occur at a roughly CONSTANT rate
# through the year, so the coincidental component X follows the slow reservoir
# P(t) and comes to dominate the recorded count in the low season, where true
# clinical malaria has collapsed. Two secondary properties: the over-count grows
# with transmission (prevalence outpaces clinical incidence as immunity limits
# disease), and treatment seeking scales the recorded volume but cancels from
# every over-count.
#
# All parameters are [ASSUMED] illustrative choices (a toy model). The scale
# kappa and the level mappings are set so the whole-year share of recorded cases
# not caused by malaria lands in the broad range reported by Dalrymple et al.
# 2017 (eLife 6:e29198), who estimate ~28% of malaria-positive fevers in
# under-fives were causally attributable to malaria in 2014. The non-malaria
# fever structure mirrors mrc-ide/malariasimulation PR #372. Not data-fitted.
#
# Run from the repository root:
#   Rscript validation/validate_nmf.R
# =============================================================================

## ---- libraries -------------------------------------------------------------
.libPaths("C:/Users/pwinskil/Documents/r_packages_arm64")  # user's package lib
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

## ===========================================================================
## 1. PARAMETERS  (identical to the JavaScript in nmf.html)
## ===========================================================================
YEAR     <- 365
DT       <- 1
PEAK_DAY <- 182   # case peak fixed at mid-year
D_INF    <- 80    # [ASSUMED] mean duration of detectable infection (days)
NMF_REF  <- 6     # [ASSUMED] reference non-malaria fever rate (episodes/child/yr)
KAPPA    <- 0.9   # [ASSUMED] relative scale for (fever rate x prevalence)

## --- Colours (match the web page) -------------------------------------------
col_clin  <- "#7b5bd6"  # true malaria cases
col_warm  <- "#e8804b"  # falsely attributed / over-count
col_ink   <- "#26303b"  # recorded line
col_blue  <- "#3d6fb4"
col_muted <- "#6b7785"

## ===========================================================================
## 2. MODEL FUNCTIONS  (identical logic to the JavaScript)
## ===========================================================================

## Seasonal shape; symmetric single peak fixed at mid-year (same as seasonality.html).
seasonal_shape <- function(t, S) {
  p     <- 1 + 8 * S + 16 * S^4
  floor <- (1 - S)^1.5
  phi   <- 2 * pi * (t - PEAK_DAY) / YEAR
  s     <- ((1 + cos(phi)) / 2)^p
  floor + (1 - floor) * s
}

## Transmission intensity (u = T/100) -> peak levels of clinical incidence and
## prevalence. Clinical saturates sooner, so prevalence outpaces it as T rises.
clin_level <- function(u) u / (u + 0.15)
prev_level <- function(u) 0.9 * u / (u + 0.5)

## Broadened/lagged/floored seasonal shape for prevalence: relax towards g with
## timescale D_INF, iterate to a periodic steady state, normalise to a peak of 1.
prev_shape <- function(S) {
  idx <- 0:(round(YEAR / DT) - 1)
  g   <- seasonal_shape(idx * DT, S)
  A <- 0.3
  for (rep in 1:8) for (i in seq_along(g)) A <- A + (g[i] - A) / D_INF * DT
  out <- numeric(length(g))
  for (i in seq_along(g)) { A <- A + (g[i] - A) / D_INF * DT; out[i] <- A }
  out / max(out)
}

## ---- simulate one scenario -------------------------------------------------
## T        transmission intensity (0..100)
## S        seasonality (0..1)
## nmfPerYr non-malaria fever rate (episodes/child/yr)
## treatPct treatment-seeking share (0..100)
simulate_nmf <- function(T, S, nmfPerYr, treatPct) {
  u     <- T / 100
  nNorm <- nmfPerYr / NMF_REF
  tau   <- treatPct / 100
  cL <- clin_level(u); vL <- prev_level(u)
  ghat <- prev_shape(S)
  idx  <- 0:(length(ghat) - 1)
  g    <- seasonal_shape(idx * DT, S)
  C    <- cL * g               # true clinical incidence
  P    <- vL * ghat            # prevalence of detectable infection
  M    <- tau * C              # true malaria cases recorded
  X    <- tau * KAPPA * nNorm * P   # falsely-attributed (coincident) cases
  R    <- M + X
  cmax <- max(C)
  peak <- C >= 0.5 * cmax      # peak season = clinical incidence >= half its max
  list(
    t = idx * DT, M = M, R = R, X = X, cmax = cmax,
    peakRatio  = if (sum(M[peak])  > 0) sum(R[peak])  / sum(M[peak])  else NA_real_,
    lowRatio   = if (sum(M[!peak]) > 0) sum(R[!peak]) / sum(M[!peak]) else NA_real_,
    falseShare = sum(X) / sum(R)
  )
}

## ===========================================================================
## 3. NUMERIC CROSS-CHECK  (compare against the web-page readouts)
## ===========================================================================
cat("\n=== validate_nmf.R : toy non-malaria-fever over-counting model ===\n\n")

# default web settings: transmission 55, seasonality 0.75, 6 fevers/child/yr, 60% seeking care
base <- simulate_nmf(55, 0.75, 6, 60)

cat("Defaults (T=55, S=0.75, NMF=6/yr, treat=60%) -- match the web readouts:\n")
cat(sprintf("  peak-season recorded / true : %.1fx   (web: 1.5x)\n", base$peakRatio))
cat(sprintf("  low-season  recorded / true : %.1fx   (web: 2.6x)\n", base$lowRatio))
cat(sprintf("  recorded cases not malaria  : %2.0f%%    (web: 51%%)\n\n", 100 * base$falseShare))

# scan the non-malaria fever rate
cat("Non-malaria fever rate scan (default transmission & season, treat=60%):\n")
cat(sprintf("  %-9s %-9s %-9s %-7s\n", "NMF/yr", "peak x", "low x", "not-mal"))
for (nmf in c(0, 2, 4, 6, 8, 12)) {
  s <- simulate_nmf(55, 0.75, nmf, 60)
  cat(sprintf("  %-9d %-9.2f %-9.2f %2.0f%%\n", nmf, s$peakRatio, s$lowRatio, 100 * s$falseShare))
}
cat("\n")

# scan transmission intensity
cat("Transmission-intensity scan (default season, NMF=6/yr, treat=60%):\n")
cat(sprintf("  %-9s %-9s %-9s %-7s\n", "T", "peak x", "low x", "not-mal"))
for (T in c(20, 40, 55, 75, 90)) {
  s <- simulate_nmf(T, 0.75, 6, 60)
  cat(sprintf("  %-9d %-9.2f %-9.2f %2.0f%%\n", T, s$peakRatio, s$lowRatio, 100 * s$falseShare))
}
cat("\n")

## ---- assertions ------------------------------------------------------------
# default readouts match the web page (to the shown precision)
stopifnot(round(base$peakRatio, 1) == 1.5)
stopifnot(round(base$lowRatio,  1) == 2.6)
stopifnot(round(base$falseShare * 100) == 51)

# the over-count is DIFFERENTIAL: worse in the low season than the peak season
stopifnot(base$lowRatio > base$peakRatio)

# no non-malaria fevers -> no inflation at all (recorded == true everywhere)
zero <- simulate_nmf(55, 0.75, 0, 60)
stopifnot(abs(zero$peakRatio - 1) < 1e-9, abs(zero$lowRatio - 1) < 1e-9, zero$falseShare < 1e-9)

# over-count rises monotonically with the fever rate and with transmission intensity
share_by_nmf <- sapply(c(0, 2, 4, 6, 8, 12), function(n) simulate_nmf(55, 0.75, n, 60)$falseShare)
stopifnot(all(diff(share_by_nmf) > 0))
share_by_T <- sapply(c(20, 40, 55, 75, 90), function(T) simulate_nmf(T, 0.75, 6, 60)$falseShare)
stopifnot(all(diff(share_by_T) > 0))

# the differential widens as the season concentrates (low-season over-count grows with S)
low_by_S <- sapply(c(0.4, 0.6, 0.8, 0.95), function(S) simulate_nmf(55, S, 6, 60)$lowRatio)
stopifnot(all(diff(low_by_S) > 0))

# treatment seeking scales recorded volume but cancels from every over-count
for (ts in c(10, 40, 100)) {
  s <- simulate_nmf(55, 0.75, 6, ts)
  stopifnot(abs(s$peakRatio  - base$peakRatio)  < 1e-9,
            abs(s$lowRatio   - base$lowRatio)   < 1e-9,
            abs(s$falseShare - base$falseShare) < 1e-9)
}

# Dalrymple anchor: at high transmission with a moderate fever rate, the majority
# of recorded cases are not caused by malaria (they estimate ~72% coincident)
dal <- simulate_nmf(90, 0.75, 10, 60)
stopifnot(dal$falseShare > 0.5)
cat(sprintf("Dalrymple-anchor scenario (T=90, NMF=10/yr): %.0f%% of recorded cases not caused by malaria.\n",
            100 * dal$falseShare))
cat("Treatment-seeking invariance confirmed: over-counts identical at 10%%, 40%%, 100%%.\n")

cat("\nAll assertions passed.\n\n")

## ===========================================================================
## 4. PLOTS
## ===========================================================================
month_starts <- cumsum(c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))
month_labs   <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

# Panel A: the default year - true, recorded, and the falsely-attributed band
# (CASE_SCALE = 100 maps the model's relative rate onto cases /1000 children /month,
#  matching the web-page axis; a display choice only, it does not affect the readouts)
CASE_SCALE <- 100
dfA <- tibble(t = base$t, true = base$M * CASE_SCALE, recorded = base$R * CASE_SCALE)
pA <- ggplot(dfA, aes(t)) +
  geom_ribbon(aes(ymin = true, ymax = recorded), fill = col_warm, alpha = 0.55) +
  geom_area(aes(y = true), fill = col_clin, alpha = 0.13) +
  geom_line(aes(y = recorded), colour = col_ink,  linewidth = 0.9) +
  geom_line(aes(y = true),     colour = col_clin, linewidth = 0.9) +
  scale_x_continuous(breaks = month_starts, labels = month_labs, expand = c(0, 0)) +
  labs(title = "True vs routine-recorded malaria (default setting)",
       subtitle = "purple = true malaria cases; dark = recorded total; orange = falsely attributed",
       x = NULL, y = "cases per 1000 children / month") +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(colour = col_muted, size = 9),
        panel.grid.minor = element_blank())

# Panel B: peak- vs low-season over-count as the fever rate rises
nmf_grid <- 0:12
dfB <- bind_rows(lapply(nmf_grid, function(n) {
  s <- simulate_nmf(55, 0.75, n, 60)
  tibble(nmf = n, `low season` = s$lowRatio, `peak season` = s$peakRatio)
})) |>
  pivot_longer(-nmf, names_to = "season", values_to = "ratio")
pB <- ggplot(dfB, aes(nmf, ratio, colour = season)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = col_muted) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.4) +
  scale_colour_manual(values = c("low season" = col_warm, "peak season" = col_blue)) +
  labs(title = "Over-count grows with the fever rate, faster in the low season",
       x = "non-malaria fevers per child per year", y = "recorded / true", colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", panel.grid.minor = element_blank())

fig <- pA / pB + plot_annotation(
  caption = "Illustrative toy model (nmf.html cross-check). Numbers show the shape of the relationship, not real-world impact.")

dir.create("figures", showWarnings = FALSE)
ggsave("figures/nmf_validation.png", fig, width = 8, height = 8, dpi = 130)
cat("Saved figures/nmf_validation.png\n")
