# =============================================================================
# validate_cascade.R
#
# Static R/ggplot reproduction of the five plots in nonlinearities.html, as an
# independent cross-check of the interactive model. It re-implements the same
# cascade
#
#        bed-net coverage  ->  EIR  ->  { PfPR , clinical incidence }
#
# in R. Since the explainer was simplified, both the EIR -> PfPR and the
# EIR -> clinical-incidence relationships come from the canonical Griffin-model
# equilibrium (mrc-ide/malariaEquilibrium), the same solver the age-distribution
# explainer uses; only the coverage -> EIR step is a separate mechanistic model.
#
# Parameter provenance is labelled throughout as one of:
#   [ASSUMED]   a fixed model choice (entomology, gonotrophic cycle, scenario)
#   [MODEL]     the Griffin equilibrium (malariaEquilibrium fitted parameters)
#   [DATA]      field observations (Smith et al. 2005 EIR-PfPR), plotted on the PfPR~EIR panel
#
# Run from the repository root:
#   & 'C:/Program Files/R-aarch64/R-4.5.2/bin/Rscript' validation/validate_cascade.R
# =============================================================================

## ---- libraries -------------------------------------------------------------
.libPaths("C:/Users/pwinskil/Documents/r_packages_arm64")  # user's package lib
suppressPackageStartupMessages({
  library(malariaEquilibrium)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(scales)
})

## ===========================================================================
## 1. PARAMETERS
## ===========================================================================

## --- Net effect on the vector, from experimental-hut trials -----------------
## Generic new pyrethroid net, no insecticide resistance (a best case).
beta <- 0.70 # [ASSUMED/EHT] blood-feeding inhibition: fraction of bites a net
#               prevents on its user
mu <- 0.45 # [ASSUMED/EHT] per-encounter mortality: fraction of mosquitoes
#               feeding on a net user that are killed
gono <- 3 # [ASSUMED]     gonotrophic cycle length (days); bridges the
#               per-encounter mortality to a daily survival

## --- Baseline vector / parasite biology (Ross-Macdonald) --------------------
p0 <- 0.90 # [ASSUMED] baseline daily mosquito survival (no nets)
eip <- 10 # [ASSUMED] extrinsic incubation period (days)

## --- Cascade controls -------------------------------------------------------
use_max <- 0.80 # [ASSUMED] coverage capped at 80% (as in the web page)

## --- Scenario plotted (matches the web page defaults) -----------------------
eir_base <- 20 # [ASSUMED] baseline annual EIR with no nets
cov_mark <- 0.50 # [ASSUMED] coverage at which the operating-point marker sits

## --- Griffin equilibrium settings (must match the JS port) ------------------
pset <- load_parameter_set() # [MODEL] standard Griffin fitted parameters
ft <- 0                      # [ASSUMED] treatment coverage off, matching the toy
age <- seq(0, 80, 0.1)       # uniform 0.1y grid to 80; the equilibrium is
#                              grid-dependent, so this MUST match the JS AGE grid

## --- Colours (match the web page) -------------------------------------------
col_accent <- "#2f8f7e"
col_blue <- "#3d6fb4"
col_clin <- "#7b5bd6"
col_warm <- "#e8804b"
col_ink <- "#26303b"
col_grey <- "#6b7785"  # matches CSS --muted
age_cols <- c("under 5" = col_clin, "5-15" = col_accent, "15+" = col_blue)

## ===========================================================================
## 2. MODEL FUNCTIONS  (identical logic to the JavaScript)
## ===========================================================================

## Survival-leverage term of vectorial capacity: p^n / -ln(p)
## (expected infectious bites a mosquito delivers given daily survival p).
f_surv <- function(p) {
  p^eip / (-log(p))
}
f0 <- f_surv(p0)

## Combined bed-net effect on transmission, relative to no nets, at coverage cov:
##   biting   a -> a0 (1 - beta*cov)            ... enters squared
##   survival p -> p0 (1 - mu*cov)^(1/gono)     ... per-cycle killing over the cycle
g_net <- function(cov) {
  a <- 1 - beta * cov
  p <- p0 * (1 - mu * cov)^(1 / gono)
  a * a * (f_surv(p) / f0)
}

## EIR with nets at coverage cov.
eir_at <- function(cov, eir0 = eir_base) {
  eir0 * g_net(cov)
}

## Aggregate one Griffin equilibrium into the quantities the explainer shows.
## Mirrors equilibrium() in nonlinearities.html / age-distribution.html:
##   pfpr        microscopy prevalence in 2-10y  = sum(pos_M[2-10]) / sum(prop[2-10])
##   inc_py      all-age episodes per person per year = sum(inc) * 365
##   incU5_py    episodes per under-5 child per year  = sum(inc[<5]) / sum(prop[<5]) * 365
##   share_*     fraction of clinical cases in each age band (sum to 1)
summarise_eq <- function(EIR) {
  m <- human_equilibrium(EIR = EIR, ft = ft, p = pset, age = age)$states
  a <- m[, "age"]
  inc <- m[, "inc"]
  tot <- sum(inc)
  u5 <- sum(inc[a < 5]); mid <- sum(inc[a >= 5 & a < 15]); ad <- sum(inc[a >= 15])
  sel <- a >= 2 & a < 10
  propU5 <- sum(m[a < 5, "prop"])
  list(
    EIR = EIR,
    pfpr = sum(m[sel, "pos_M"]) / sum(m[sel, "prop"]),
    inc_py = tot * 365,
    incU5_py = if (propU5 > 0) u5 / propU5 * 365 else 0,
    share_u5 = u5 / tot, share_mid = mid / tot, share_ad = ad / tot
  )
}

## Vectorised helper: a data frame of the equilibrium over a vector of EIRs.
eq_grid <- function(eirs) {
  do.call(rbind, lapply(eirs, function(E) as.data.frame(summarise_eq(E))))
}

## ===========================================================================
## 3. NUMERIC CROSS-CHECK at the plotted scenario
## ===========================================================================
eq0 <- as.data.frame(summarise_eq(eir_base))         # baseline (no nets)
eqC <- as.data.frame(summarise_eq(eir_at(cov_mark))) # at the marker coverage
cat(sprintf(
  "\nScenario: baseline EIR = %.0f, coverage = %.0f%%\n", eir_base, 100 * cov_mark
))
cat(sprintf(
  "  EIR reduction at %.0f%% coverage : %.0f%%\n",
  100 * use_max, 100 * (1 - g_net(use_max))
))
cat(sprintf("  EIR with nets (%.0f%%)           : %.1f\n", 100 * cov_mark, eir_at(cov_mark)))
cat(sprintf("  baseline PfPR (2-10)            : %.0f%%\n", 100 * eq0$pfpr))
cat(sprintf("  PfPR with nets                  : %.0f%%\n", 100 * eqC$pfpr))
cat(sprintf("  prevalence reduction            : %.0f%%\n", 100 * (1 - eqC$pfpr / eq0$pfpr)))
cat(sprintf(
  "  under-5 clinical reduction      : %.0f%%\n",
  100 * (1 - eqC$incU5_py / eq0$incU5_py)
))
cat("  These should match the explainer's readouts (same solver, grid, ft = 0).\n")

## ===========================================================================
## 4. THE FIVE PLOTS
## ===========================================================================
base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = col_grey, size = 9.5)
  )

## Upper EIR bound, shared with the web page (nonlinearities.html LUT_EMAX): capped just past
## the incidence plateau, before the Griffin equilibrium's extreme-EIR upturn.
eir_top <- 300
cov_grid <- seq(0, use_max, length.out = 120)
eir_curve <- 10^seq(log10(0.02), log10(eir_top), length.out = 120)
inc_ymax <- NA_real_ # filled after the grids are computed

## Equilibrium over the EIR curve (general relationships) and the coverage grid.
eq_eir <- eq_grid(eir_curve)
eq_cov <- eq_grid(eir_at(cov_grid))
eq_cov$cov <- cov_grid
inc_ymax <- max(c(eq_eir$inc_py, eq_cov$inc_py)) * 1.1

## Field observations for the PfPR~EIR panel: Smith et al. (2005) annual EIR vs PfPR in
## children <15 across African sites [DATA]. Only sites within the plotted range are shown
## (as on the web page); a handful of high-transmission sites lie beyond the EIR 300 cap.
load("data/EIR_prev_hay2005.RData") # provides data frame EIR_prev_hay2005
smith <- as.data.frame(EIR_prev_hay2005)
names(smith) <- c("eir", "pfpr")
smith <- subset(smith, eir >= 0.02 & eir <= eir_top)

## Stack the age-band contributions (inc * share) into long form for geom_area.
stack_long <- function(df, xcol) {
  x <- df[[xcol]]
  rbind(
    data.frame(x = x, band = "under 5", y = df$inc_py * df$share_u5),
    data.frame(x = x, band = "5-15", y = df$inc_py * df$share_mid),
    data.frame(x = x, band = "15+", y = df$inc_py * df$share_ad)
  ) %>% mutate(band = factor(band, levels = c("under 5", "5-15", "15+")))
}

## --- Plot 1: PfPR vs EIR (Griffin equilibrium; log EIR) ---------------------
p1 <- ggplot(eq_eir, aes(EIR, pfpr)) +
  geom_point(data = smith, aes(eir, pfpr), colour = col_grey, alpha = 0.30, size = 1.3) +
  geom_line(colour = col_blue, linewidth = 1) +
  geom_point(data = eqC, aes(EIR, pfpr), colour = col_warm, size = 3) +
  scale_x_log10(
    limits = c(0.02, eir_top), oob = scales::squish,
    breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")
  ) +
  scale_y_continuous(labels = percent, limits = c(0, 1), oob = scales::squish) +
  labs(title = "1. PfPR vs EIR", subtitle = "line: Griffin equilibrium; points: Smith 2005 (PfPR<15)",
       x = "EIR (log scale)", y = "prevalence (PfPR 2-10)") +
  base_theme

## --- Plot 2: incidence vs PfPR (overall line + age ribbon; parametric) ------
ribbon_pf <- stack_long(eq_eir, "pfpr")
p2 <- ggplot() +
  geom_area(data = ribbon_pf, aes(x, y, fill = band), alpha = 0.34, position = "stack") +
  geom_line(data = eq_eir, aes(pfpr, inc_py), colour = col_ink, linewidth = 1) +
  geom_point(data = eqC, aes(pfpr, inc_py), colour = col_warm, size = 3) +
  scale_fill_manual(values = age_cols, name = NULL) +
  scale_x_continuous(labels = percent, limits = c(0, 0.80), oob = scales::squish) +
  coord_cartesian(ylim = c(0, inc_ymax)) +
  labs(title = "3. Incidence vs PfPR", subtitle = "all-age; bands = age groups",
       x = "PfPR (2-10)", y = "incidence (/person/yr)") +
  base_theme + theme(legend.position = c(0.20, 0.78),
                     legend.background = element_rect(fill = "white", colour = NA))

## --- Plot 3: incidence vs EIR (overall line + age ribbon; log EIR) ----------
ribbon_eir <- stack_long(eq_eir, "EIR")
p3 <- ggplot() +
  geom_area(data = ribbon_eir, aes(x, y, fill = band), alpha = 0.34, position = "stack") +
  geom_line(data = eq_eir, aes(EIR, inc_py), colour = col_ink, linewidth = 1) +
  geom_point(data = eqC, aes(EIR, inc_py), colour = col_warm, size = 3) +
  scale_fill_manual(values = age_cols, name = NULL, guide = "none") +
  scale_x_log10(
    limits = c(0.02, eir_top), oob = scales::squish,
    breaks = c(0.1, 1, 10, 100), labels = c("0.1", "1", "10", "100")
  ) +
  coord_cartesian(ylim = c(0, inc_ymax)) +
  labs(title = "2. Incidence vs EIR", subtitle = "all-age; rises then plateaus",
       x = "EIR (log scale)", y = "incidence (/person/yr)") +
  base_theme

## --- Plot 4: EIR vs coverage (mechanistic net effect) -----------------------
p4 <- ggplot(
  data.frame(cov = cov_grid, eir = eir_at(cov_grid)), aes(cov, eir)
) +
  geom_line(colour = col_accent, linewidth = 1) +
  geom_point(data = data.frame(cov = cov_mark, eir = eir_at(cov_mark)),
             aes(cov, eir), colour = col_warm, size = 3) +
  scale_x_continuous(labels = percent) +
  labs(title = "4. EIR vs coverage", subtitle = "barrier (a^2) + killing (p^n)",
       x = "bed-net coverage", y = "EIR") +
  base_theme

## --- Plot 5: incidence vs coverage (full cascade; overall line + ribbon) ----
ribbon_cov <- stack_long(eq_cov, "cov")
p5 <- ggplot() +
  geom_area(data = ribbon_cov, aes(x, y, fill = band), alpha = 0.34, position = "stack") +
  geom_line(data = eq_cov, aes(cov, inc_py), colour = col_ink, linewidth = 1) +
  geom_point(data = data.frame(cov = cov_mark, inc = eqC$inc_py),
             aes(cov, inc), colour = col_warm, size = 3) +
  scale_fill_manual(values = age_cols, name = NULL, guide = "none") +
  scale_x_continuous(labels = percent) +
  coord_cartesian(ylim = c(0, inc_ymax)) +
  labs(title = "5. Incidence vs coverage", subtitle = "full cascade",
       x = "bed-net coverage", y = "incidence (/person/yr)") +
  base_theme

## ===========================================================================
## 5. COMBINE & SAVE
## ===========================================================================
combined <- (p1 | p3 | p2) /
  (p4 | p5 | plot_spacer()) +
  plot_annotation(
    title = sprintf(
      "Bed-net cascade via the Griffin equilibrium (baseline EIR = %.0f, marker at %.0f%% coverage)",
      eir_base, 100 * cov_mark
    ),
    theme = theme(plot.title = element_text(face = "bold"))
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/cascade_validation.png", combined, width = 15, height = 8.5, dpi = 130)
cat("\nSaved figures/cascade_validation.png\n")
