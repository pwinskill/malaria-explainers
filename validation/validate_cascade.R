# =============================================================================
# validate_cascade.R
#
# Static R/ggplot reproduction of the four plots in index.html, as an
# independent cross-check of the interactive model. It re-implements the same
# cascade
#
#        bed-net coverage  ->  EIR  ->  { PfPR , clinical incidence }
#
# in R, reading the underlying data where the web page uses it, and refitting /
# documenting every parameter.
#
# Parameter provenance is labelled throughout as one of:
#   [DATA]      read from a file (Hay EIR-PfPR; Battle incidence-prevalence)
#   [FITTED]    estimated here from the data (the EIR->PfPR curve)
#   [DIGITISED] read off a published figure (Cameron age curves, Fig. 1)
#   [ASSUMED]   a fixed model choice (entomology, gonotrophic cycle, scenario)
#
# Run from the repository root (where the data files live):
#   Rscript validation/validate_cascade.R
# =============================================================================

## ---- libraries -------------------------------------------------------------
#.libPaths("C:/Users/pwinskil/Documents/r_packages_arm64")  # user's package lib
suppressPackageStartupMessages({
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

## --- EIR -> PfPR saturating curve -------------------------------------------
## Form: PfPR = 1 - exp(-a * EIR^b). a, b are [FITTED] to the Hay data below;
## these starting values are also the constants used in the web page.
hay_a_start <- 0.455
hay_b_start <- 0.216

## --- Cameron et al. (2015) PfPR2-10 -> clinical incidence -------------------
## [DIGITISED] from Fig. 1 (Griffin individual-based model, low seasonality),
## absolute incidence (cases person^-1 yr^-1) vs PfPR2-10, by age group.
cam_young <- data.frame( # 0-5y  (convex)
  pfpr = c(0, .05, .10, .15, .20, .25, .30, .35, .40, .45, .50, .55, .60),
  inc = c(0, .055, .133, .232, .354, .498, .653, .808, .973, 1.139, 1.305, 1.471, 1.637)
)
cam_old <- data.frame( # 5-15y (saturating)
  pfpr = c(0, .05, .10, .15, .20, .25, .30, .35, .45, .50, .55),
  inc = c(0, .077, .177, .277, .387, .476, .553, .608, .675, .697, .719)
)
cam_adult <- data.frame( # >15y  (low, peaks at intermediate PfPR)
  pfpr = c(0, .05, .10, .15, .20, .35, .40, .45, .50),
  inc = c(0, .055, .122, .166, .21, .232, .221, .221, .21)
)

## --- Colours (match the web page) -------------------------------------------
col_accent <- "#2f8f7e"
col_blue <- "#3d6fb4"
col_clin <- "#7b5bd6"
col_warm <- "#e8804b"
col_grey <- "#5b6573"
age_cols <- c("0-5y" = col_clin, "5-15y" = col_accent, ">15y" = col_blue)

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

## Monotone cubic (Fritsch-Carlson) interpolation of a digitised Cameron curve,
## matching the web page's PCHIP. Beyond the last digitised point we extrapolate:
##   tail_exp = a number -> power-law  last_y * (pr/last_x)^tail_exp  (young: 1.2)
##   tail_exp = NA        -> flat       last_y                        (older / adult)
make_cam <- function(pts, tail_exp = NA_real_) {
  sp <- splinefun(pts$pfpr, pts$inc, method = "monoH.FC")
  last_x <- tail(pts$pfpr, 1)
  last_y <- tail(pts$inc, 1)
  function(pr) {
    pr <- pmax(pr, 0)
    out <- sp(pmin(pr, last_x)) # monotone interp within range
    beyond <- pr > last_x
    if (any(beyond)) {
      out[beyond] <- if (!is.na(tail_exp)) {
        last_y * (pr[beyond] / last_x)^tail_exp
      } else {
        last_y
      }
    }
    out[pr <= 0] <- 0
    pmax(out, 0)
  }
}
inc_young <- make_cam(cam_young, tail_exp = 1.2)
inc_old <- make_cam(cam_old, tail_exp = NA)
inc_adult <- make_cam(cam_adult, tail_exp = NA)

## ===========================================================================
## 3. EIR -> PfPR : read the Hay data and REFIT the saturating curve  [DATA/FITTED]
## ===========================================================================
## Source: P. falciparum parasite rate in children (<15 y) vs annual EIR across
## African sites (Smith et al. 2005, Nature; data file named for senior author
## Hay). 130 site-level points.
load("data/EIR_prev_hay2005.RData") # provides data frame EIR_prev_hay2005
hay <- as.data.frame(EIR_prev_hay2005)
names(hay) <- c("eir", "pfpr")

## Refit PfPR = 1 - exp(-a * EIR^b) to the data (EIR > 0; the curve is 0 at EIR 0).
hay_fit <- nls(
  pfpr ~ 1 - exp(-a * eir^b),
  data = subset(hay, eir > 0),
  start = list(a = hay_a_start, b = hay_b_start)
)
hay_a <- coef(hay_fit)[["a"]]
hay_b <- coef(hay_fit)[["b"]]
cat(sprintf(
  "[FITTED] EIR->PfPR:  a = %.3f, b = %.3f  (web page uses %.3f, %.3f)\n",
  hay_a,
  hay_b,
  hay_a_start,
  hay_b_start
))

## Now that a, b are known, define the prevalence map and the cascade prevalence.
pfpr_from_eir <- function(eir) {
  ifelse(eir <= 0, 0, 1 - exp(-hay_a * eir^hay_b))
}
pf_at <- function(cov, eir0 = eir_base) {
  pfpr_from_eir(eir_at(cov, eir0))
}

## ===========================================================================
## 4. PfPR -> incidence : read & process the Battle field data  [DATA]
## ===========================================================================
## Source: matched incidence-prevalence records (Battle et al. 2015). We take the
## P. falciparum, sub-Saharan Africa, Cameron-calibration subset, one point per
## record, and bin each by its reported incidence age range.
battle <- read.csv(
  "data/PfPvAllData01042015_AgeStand.csv",
  fileEncoding = "latin1",
  stringsAsFactors = FALSE
)
names(battle) <- tolower(names(battle)) # lowercase all column names

battle_pts <- battle %>%
  filter(species == "Pf", region == "Africa+", cameron == "Yes") %>%
  transmute(
    pfpr = suppressWarnings(as.numeric(pfpr2_10)), # standardised PfPR2-10
    inc = suppressWarnings(as.numeric(inc)) / 1000, # raw inc is per 1000 PYO -> per person-yr
    lar = suppressWarnings(as.numeric(inc_lar)),
    uar = suppressWarnings(as.numeric(inc_uar))
  ) %>%
  filter(pfpr > 0, inc > 0) %>%
  mutate(
    age = case_when(
      uar <= 6 & lar < 6 ~ "0-5y", # young children
      lar >= 4 & uar <= 16 ~ "5-15y", # older children
      lar >= 15 ~ ">15y", # adults
      TRUE ~ "mixed" # spans groups -> dropped
    )
  ) %>%
  filter(age != "mixed") %>%
  mutate(age = factor(age, levels = c("0-5y", "5-15y", ">15y")))

cat(sprintf(
  "[DATA] Battle points: %d (0-5y), %d (5-15y), %d (>15y)\n",
  sum(battle_pts$age == "0-5y"),
  sum(battle_pts$age == "5-15y"),
  sum(battle_pts$age == ">15y")
))

## ===========================================================================
## 5. NUMERIC CROSS-CHECK at the plotted scenario
## ===========================================================================
pf0 <- pf_at(0, eir_base) # baseline PfPR (no nets)
pf_cov <- pf_at(cov_mark, eir_base) # PfPR at the marker coverage
cat(sprintf(
  "\nScenario: baseline EIR = %.0f, coverage = %.0f%%\n",
  eir_base,
  100 * cov_mark
))
cat(sprintf(
  "  EIR reduction at %.0f%% coverage : %.0f%%\n",
  100 * use_max,
  100 * (1 - g_net(use_max))
))
cat(sprintf("  baseline PfPR (<15)             : %.0f%%\n", 100 * pf0))
cat(sprintf("  PfPR with nets                  : %.0f%%\n", 100 * pf_cov))
cat(sprintf(
  "  under-5 clinical reduction      : %.0f%%\n",
  100 * (1 - inc_young(pf_cov) / inc_young(pf0))
))

## ===========================================================================
## 6. THE FIVE PLOTS
## ===========================================================================
base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(color = col_grey, size = 9.5)
  )

cov_grid <- seq(0, use_max, length.out = 200)
pfpr_grid <- seq(0.001, 0.80, length.out = 300)

## --- Plot 1: EIR vs coverage ------------------------------------------------
p1 <- ggplot(
  data.frame(cov = cov_grid, eir = eir_at(cov_grid)),
  aes(cov, eir)
) +
  geom_line(colour = col_accent, linewidth = 1) +
  geom_point(
    data = data.frame(cov = cov_mark, eir = eir_at(cov_mark)),
    aes(cov, eir),
    colour = col_warm,
    size = 3
  ) +
  scale_x_continuous(labels = percent) +
  labs(
    title = "1. EIR vs coverage",
    subtitle = "barrier (a^2) + killing (p^n)",
    x = "bed-net coverage",
    y = "EIR"
  ) +
  base_theme

## --- Plot 2: PfPR vs EIR (data + fitted curve) ------------------------------
eir_curve <- 10^seq(log10(0.02), log10(1000), length.out = 300)
p2 <- ggplot() +
  geom_point(
    data = subset(hay, eir >= 0.02),
    aes(eir, pfpr),
    colour = col_grey,
    alpha = 0.45,
    size = 1.3
  ) +
  geom_line(
    data = data.frame(eir = eir_curve, pfpr = pfpr_from_eir(eir_curve)),
    aes(eir, pfpr),
    colour = col_blue,
    linewidth = 1
  ) +
  geom_point(
    data = data.frame(eir = eir_at(cov_mark), pfpr = pf_at(cov_mark)),
    aes(eir, pfpr),
    colour = col_warm,
    size = 3
  ) +
  scale_x_log10(
    limits = c(0.02, 1000),
    oob = scales::squish,
    breaks = c(0.1, 1, 10, 100, 1000),
    labels = c("0.1", "1", "10", "100", "1000")
  ) +
  scale_y_continuous(labels = percent, limits = c(0, 1), oob = scales::squish) +
  labs(
    title = "2. PfPR vs EIR",
    subtitle = "dots: Smith et al. 2005 (children <15); curve: fitted saturating model",
    x = "EIR (log scale)",
    y = "prevalence (PfPR)"
  ) +
  base_theme

## --- Plot 4: incidence vs EIR (composed; log EIR, linear incidence) ---------
## incidence(EIR) = cameron_age( pfpr_from_eir(EIR) ): the fork's second arm.
curve_ie <- bind_rows(
  data.frame(eir = eir_curve, inc = inc_young(pfpr_from_eir(eir_curve)), age = "0-5y"),
  data.frame(eir = eir_curve, inc = inc_old(pfpr_from_eir(eir_curve)), age = "5-15y"),
  data.frame(eir = eir_curve, inc = inc_adult(pfpr_from_eir(eir_curve)), age = ">15y")
) %>%
  mutate(age = factor(age, levels = c("0-5y", "5-15y", ">15y")))

marker_ie <- data.frame(
  eir = eir_at(cov_mark),
  inc = c(inc_young(pf_cov), inc_old(pf_cov), inc_adult(pf_cov)),
  age = factor(c("0-5y", "5-15y", ">15y"), levels = c("0-5y", "5-15y", ">15y"))
)

p4 <- ggplot() +
  geom_line(data = curve_ie, aes(eir, inc, colour = age), linewidth = 1) +
  geom_point(data = marker_ie, aes(eir, inc), colour = col_warm, size = 2.6) +
  scale_colour_manual(values = age_cols, name = NULL) +
  scale_x_log10(
    limits = c(0.02, 1000),
    oob = scales::squish,
    breaks = c(0.1, 1, 10, 100, 1000),
    labels = c("0.1", "1", "10", "100", "1000")
  ) +
  scale_y_continuous(limits = c(0, 3), oob = scales::squish) +
  labs(
    title = "4. Incidence vs EIR",
    subtitle = "log EIR, linear incidence: older ages plateau, young keep rising",
    x = "EIR (log scale)",
    y = "incidence (/person/yr)"
  ) +
  base_theme +
  theme(legend.position = "none")

## --- Plot 3: incidence vs PfPR (3 age curves + field data; linear) ----------
curve_pfpr <- bind_rows(
  data.frame(pfpr = pfpr_grid, inc = inc_young(pfpr_grid), age = "0-5y"),
  data.frame(pfpr = pfpr_grid, inc = inc_old(pfpr_grid), age = "5-15y"),
  data.frame(pfpr = pfpr_grid, inc = inc_adult(pfpr_grid), age = ">15y")
) %>%
  mutate(age = factor(age, levels = c("0-5y", "5-15y", ">15y")))

marker_pfpr <- data.frame(
  pfpr = pf_cov,
  inc = c(inc_young(pf_cov), inc_old(pf_cov), inc_adult(pf_cov)),
  age = factor(c("0-5y", "5-15y", ">15y"), levels = c("0-5y", "5-15y", ">15y"))
)

## Linear incidence axis (0-3) to match plot 4; the highest field records
## (incidence > 3) are squished to the top edge, as on the web page.
p3 <- ggplot() +
  geom_point(
    data = battle_pts,
    aes(pfpr, inc, colour = age),
    alpha = 0.30,
    size = 1.1
  ) +
  geom_line(data = curve_pfpr, aes(pfpr, inc, colour = age), linewidth = 1) +
  geom_point(data = marker_pfpr, aes(pfpr, inc), colour = col_warm, size = 2.6) +
  scale_colour_manual(values = age_cols, name = NULL) +
  scale_x_continuous(labels = percent, limits = c(0, 0.80), oob = scales::squish) +
  scale_y_continuous(limits = c(0, 3), oob = scales::squish) +
  labs(
    title = "3. Incidence vs PfPR",
    subtitle = "curves: Cameron 2015 by age; dots: Battle 2015 (highest run off top)",
    x = "PfPR (2-10)",
    y = "incidence (/person/yr)"
  ) +
  base_theme +
  theme(
    legend.position = c(0.18, 0.80),
    legend.background = element_rect(fill = "white", colour = NA)
  )

## --- Plot 5: incidence vs coverage (3 age curves) ---------------------------
curve_cov <- bind_rows(
  data.frame(cov = cov_grid, inc = inc_young(pf_at(cov_grid)), age = "0-5y"),
  data.frame(cov = cov_grid, inc = inc_old(pf_at(cov_grid)), age = "5-15y"),
  data.frame(cov = cov_grid, inc = inc_adult(pf_at(cov_grid)), age = ">15y")
) %>%
  mutate(age = factor(age, levels = c("0-5y", "5-15y", ">15y")))

marker_cov <- data.frame(
  cov = cov_mark,
  inc = c(inc_young(pf_cov), inc_old(pf_cov), inc_adult(pf_cov)),
  age = factor(c("0-5y", "5-15y", ">15y"), levels = c("0-5y", "5-15y", ">15y"))
)

## fixed 0-3 incidence axis, matching plots 3 and 4 (as on the web page).
p5 <- ggplot() +
  geom_line(data = curve_cov, aes(cov, inc, colour = age), linewidth = 1) +
  geom_point(data = marker_cov, aes(cov, inc), colour = col_warm, size = 2.6) +
  scale_colour_manual(values = age_cols, name = NULL) +
  scale_x_continuous(labels = percent) +
  scale_y_continuous(limits = c(0, 3), oob = scales::squish) +
  labs(
    title = "5. Incidence vs coverage",
    subtitle = "absolute incidence by age group",
    x = "bed-net coverage",
    y = "incidence (/person/yr)"
  ) +
  base_theme +
  theme(legend.position = "none")

## ===========================================================================
## 7. COMBINE & SAVE
## ===========================================================================
combined <- (p1 | p2 | p3) /
  (p4 | p5 | plot_spacer()) +
  plot_annotation(
    title = sprintf(
      "Bed-net cascade (baseline EIR = %.0f, marker at %.0f%% coverage)",
      eir_base,
      100 * cov_mark
    ),
    theme = theme(plot.title = element_text(face = "bold"))
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/cascade_validation.png", combined, width = 15, height = 8.5, dpi = 130)
cat("\nSaved figures/cascade_validation.png\n")
