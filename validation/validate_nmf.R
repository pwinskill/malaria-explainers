# =============================================================================
# validate_nmf.R
#
# Static R/ggplot reproduction of the interactive toy model in nmf.html, as an
# independent cross-check of the JavaScript. It re-implements the same three
# quantities through a single year (day-of-year t = 0..364):
#
#   C(t)  true seasonal clinical malaria in young children (peak fixed at 1)
#   A(t)  prevalence of detectable (asymptomatic) infection: the periodic steady
#         state of an SIS relation  dA/dt = beta*C(t)*(1-A) - A/D  driven by C,
#         so A is a broadened, lagged, floored version of C
#   R(t)  routine-recorded cases = C(t) + kappa*n*A(t), where the second term is
#         coincidental positives (non-malaria fevers that happen to test positive)
#
# The point of the model: non-malaria fevers occur at a roughly CONSTANT rate
# through the year, so the coincidental component follows the slow reservoir A(t)
# and comes to dominate the recorded count in the low season, where true clinical
# malaria has collapsed. Routine confirmed-case counts therefore over-state
# malaria, most heavily off-peak.
#
# All parameters are [ASSUMED] illustrative choices (a toy model). The scale
# kappa and reservoir mapping are set so the whole-year share of recorded cases
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
YEAR    <- 365
DT      <- 1
D_INF   <- 80    # [ASSUMED] mean duration of detectable infection (days)
NMF_REF <- 6     # [ASSUMED] reference non-malaria fever rate (episodes/child/yr)
KAPPA   <- 1.0   # [ASSUMED] relative scale for (fever rate x prevalence)
RISE_FRAC <- c(0.67, 0.5, 0.33)   # [slow, medium, fast], takeoff index 0..2

## --- Colours (match the web page) -------------------------------------------
col_clin  <- "#7b5bd6"  # true clinical malaria
col_warm  <- "#e8804b"  # falsely attributed / over-count
col_ink   <- "#26303b"  # recorded line
col_blue  <- "#3d6fb4"
col_muted <- "#6b7785"

## ===========================================================================
## 2. MODEL FUNCTIONS  (identical logic to the JavaScript)
## ===========================================================================

## Seasonal true-clinical curve; peak fixed at 1 (same as seasonality.html).
seasonal_cases <- function(t, S, peakDay, riseFrac) {
  p     <- 1 + 8 * S + 16 * S^4
  floor <- (1 - S)^1.5
  R <- riseFrac * YEAR; Ff <- YEAR - R
  d   <- ((t - peakDay) %% YEAR + YEAR) %% YEAR
  phi <- ifelse(d <= Ff, pi * d / Ff, pi * (d - YEAR) / R)
  s   <- ((1 + cos(phi)) / 2)^p
  floor + (1 - floor) * s
}

## Transmission-intensity slider (0..100) -> SIS transmission coefficient beta.
beta_of <- function(E) 0.003 * exp(E / 33)

## Prevalence of detectable infection: periodic steady state of
## dA/dt = beta*C(t)*(1-A) - A/D_INF, iterated around the year to convergence.
asymp_prev <- function(S, peakDay, riseFrac, beta) {
  idx <- 0:(round(YEAR / DT) - 1)
  C   <- seasonal_cases(idx * DT, S, peakDay, riseFrac)
  A <- 0.3
  for (rep in 1:6) for (i in seq_along(C)) A <- A + (beta * C[i] * (1 - A) - A / D_INF) * DT
  out <- numeric(length(C))
  for (i in seq_along(C)) { A <- A + (beta * C[i] * (1 - A) - A / D_INF) * DT; out[i] <- A }
  out
}

## ---- simulate one scenario -------------------------------------------------
## S        seasonality (0..1)
## peakDay  day-of-year of the case peak (1..365)
## riseFrac fraction of the year the rise into the peak spans
## E        transmission intensity (0..100)
## nmfPerYr non-malaria fever rate (episodes/child/yr)
simulate_nmf <- function(S, peakDay, riseFrac, E, nmfPerYr) {
  nNorm <- nmfPerYr / NMF_REF
  A     <- asymp_prev(S, peakDay, riseFrac, beta_of(E))
  idx   <- 0:(length(A) - 1)
  c0    <- seasonal_cases(idx * DT, S, peakDay, riseFrac)
  cmax  <- max(c0)
  fa    <- KAPPA * nNorm * A          # falsely-attributed (coincident) cases
  rec   <- c0 + fa
  peak  <- c0 >= 0.5 * cmax           # peak season = true incidence >= half its max
  list(
    t = idx * DT, c0 = c0, rec = rec, A = A, fa = fa, cmax = cmax,
    annualRatio = sum(rec) / sum(c0),
    peakRatio   = if (sum(c0[peak])  > 0) sum(rec[peak])  / sum(c0[peak])  else NA_real_,
    lowRatio    = if (sum(c0[!peak]) > 0) sum(rec[!peak]) / sum(c0[!peak]) else NA_real_,
    falseShare  = sum(fa) / sum(rec)
  )
}

## ===========================================================================
## 3. NUMERIC CROSS-CHECK  (compare against the web-page readouts)
## ===========================================================================
cat("\n=== validate_nmf.R : toy non-malaria-fever over-counting model ===\n\n")

# default web settings: seasonality 0.75, peak day 258 (Sep), medium take-off,
# transmission 55, 6 non-malaria fevers/child/yr
base <- simulate_nmf(0.75, 258, RISE_FRAC[2], 55, 6)

cat(sprintf("Defaults (S=0.75, peak=258, medium, E=55, NMF=6/yr) -- match the web readouts:\n"))
cat(sprintf("  peak-season recorded / true : %.1fx   (web: 1.4x)\n", base$peakRatio))
cat(sprintf("  low-season  recorded / true : %.1fx   (web: 2.4x)\n", base$lowRatio))
cat(sprintf("  recorded cases not malaria  : %2.0f%%    (web: 47%%)\n", 100 * base$falseShare))
cat(sprintf("  (annual recorded / true     : %.2fx)\n\n", base$annualRatio))

# scan the non-malaria fever rate at the default setting
cat("Non-malaria fever rate scan (default season & transmission):\n")
cat(sprintf("  %-9s %-9s %-9s %-7s\n", "NMF/yr", "peak x", "low x", "not-mal"))
for (nmf in c(0, 2, 4, 6, 8, 12)) {
  s <- simulate_nmf(0.75, 258, RISE_FRAC[2], 55, nmf)
  cat(sprintf("  %-9d %-9.2f %-9.2f %2.0f%%\n", nmf, s$peakRatio, s$lowRatio, 100 * s$falseShare))
}
cat("\n")

# scan transmission intensity (reservoir size) at the default fever rate
cat("Transmission-intensity scan (default season, NMF=6/yr):\n")
cat(sprintf("  %-9s %-11s %-9s %-9s %-7s\n", "E", "mean A(t)", "peak x", "low x", "not-mal"))
for (E in c(20, 40, 55, 75, 90)) {
  s <- simulate_nmf(0.75, 258, RISE_FRAC[2], E, 6)
  cat(sprintf("  %-9d %-11.3f %-9.2f %-9.2f %2.0f%%\n", E, mean(s$A), s$peakRatio, s$lowRatio, 100 * s$falseShare))
}
cat("\n")

## ---- assertions ------------------------------------------------------------
# default readouts match the web page (to the shown precision)
stopifnot(round(base$peakRatio, 1) == 1.4)
stopifnot(round(base$lowRatio,  1) == 2.4)
stopifnot(round(base$falseShare * 100) == 47)

# the over-count is DIFFERENTIAL: worse in the low season than the peak season
stopifnot(base$lowRatio > base$peakRatio)

# no non-malaria fevers -> no inflation at all (recorded == true everywhere)
zero <- simulate_nmf(0.75, 258, RISE_FRAC[2], 55, 0)
stopifnot(abs(zero$peakRatio - 1) < 1e-9, abs(zero$lowRatio - 1) < 1e-9, zero$falseShare < 1e-9)

# over-count rises monotonically with the fever rate and with transmission intensity
share_by_nmf <- sapply(c(0, 2, 4, 6, 8, 12), function(n) simulate_nmf(0.75, 258, RISE_FRAC[2], 55, n)$falseShare)
stopifnot(all(diff(share_by_nmf) > 0))
share_by_E <- sapply(c(20, 40, 55, 75, 90), function(E) simulate_nmf(0.75, 258, RISE_FRAC[2], E, 6)$falseShare)
stopifnot(all(diff(share_by_E) > 0))

# the differential widens as the season concentrates (low-season over-count grows with S)
low_by_S <- sapply(c(0.4, 0.6, 0.8, 0.95), function(S) simulate_nmf(S, 258, RISE_FRAC[2], 55, 6)$lowRatio)
stopifnot(all(diff(low_by_S) > 0))

# Dalrymple anchor: at high transmission with a moderate fever rate, the majority
# of recorded cases are not caused by malaria (they estimate ~72% coincident)
dal <- simulate_nmf(0.75, 258, RISE_FRAC[2], 90, 8)
stopifnot(dal$falseShare > 0.5)
cat(sprintf("Dalrymple-anchor scenario (E=90, NMF=8/yr): %.0f%% of recorded cases not caused by malaria.\n", 100 * dal$falseShare))

cat("\nAll assertions passed.\n\n")

## ===========================================================================
## 4. PLOTS
## ===========================================================================
month_starts <- cumsum(c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30))
month_labs   <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

# Panel A: the default year - true, recorded, and the falsely-attributed band
dfA <- tibble(t = base$t, true = base$c0, recorded = base$rec)
pA <- ggplot(dfA, aes(t)) +
  geom_ribbon(aes(ymin = true, ymax = recorded), fill = col_warm, alpha = 0.55) +
  geom_area(aes(y = true), fill = col_clin, alpha = 0.13) +
  geom_line(aes(y = recorded), colour = col_ink,  linewidth = 0.9) +
  geom_line(aes(y = true),     colour = col_clin, linewidth = 0.9) +
  scale_x_continuous(breaks = month_starts, labels = month_labs, expand = c(0, 0)) +
  labs(title = "True vs routine-recorded malaria (default setting)",
       subtitle = "purple = true clinical malaria; dark = recorded; orange = falsely attributed",
       x = NULL, y = "relative cases") +
  theme_minimal(base_size = 11) +
  theme(plot.subtitle = element_text(colour = col_muted, size = 9),
        panel.grid.minor = element_blank())

# Panel B: peak- vs low-season over-count as the fever rate rises
nmf_grid <- 0:12
dfB <- bind_rows(lapply(nmf_grid, function(n) {
  s <- simulate_nmf(0.75, 258, RISE_FRAC[2], 55, n)
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
