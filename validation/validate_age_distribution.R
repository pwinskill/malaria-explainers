# =============================================================================
# validate_age_distribution.R
#
# Cross-check the "age distribution of malaria burden" explainer against the
# canonical Griffin-model equilibrium (mrc-ide/malariaEquilibrium).
#
# The explainer's toy is a direct JavaScript port of malariaEquilibrium's
# human_equilibrium(), so this script runs the same solver in R over a range of
# transmission intensities (EIR) and checks the quantities the page reports:
#   1. Parasite prevalence (microscopy, 2-10y) rises with EIR.
#   2. The age distribution of CLINICAL cases shifts to younger children as EIR
#      rises: the under-5 share climbs from well below a fifth at low
#      transmission to above three-fifths at high transmission, while the 15+
#      share collapses. (Griffin et al. 2014, Nat Commun 5:3136.)
#   3. All-age clinical incidence (episodes per person-year) rises then saturates.
#
# The JS port uses a uniform 0.1-year age grid 0-80 and treatment ft = 0; this
# script matches that exactly, so the printed table should equal the page's
# readouts to rounding.
#
# Run (from the repository root, so figures/ resolves):
#   & 'C:/Program Files/R-aarch64/R-4.5.2/bin/Rscript' validation/validate_age_distribution.R
# =============================================================================

.libPaths('C:/Users/pwinskil/Documents/r_packages_arm64')
suppressMessages({
  library(malariaEquilibrium)
  library(ggplot2)
})

p   <- load_parameter_set()          # standard Griffin fitted parameters
ft  <- 0                             # treatment coverage off, matching the toy
age <- seq(0, 80, 0.1)               # uniform grid, matching the JS port

# ----------------------------------------------------------------------------
# Aggregate one equilibrium into the quantities the explainer shows.
#   case shares use inc directly (inc is cases per total-population; it already
#   embeds the demographic weighting), prevalence uses pos_M / prop over 2-10y.
# ----------------------------------------------------------------------------
summarise_eq <- function(EIR) {
  m   <- human_equilibrium(EIR = EIR, ft = ft, p = p, age = age)$states
  a   <- m[, "age"]
  inc <- m[, "inc"]
  tot <- sum(inc)
  u5  <- sum(inc[a < 5]) / tot
  mid <- sum(inc[a >= 5 & a < 15]) / tot
  ad  <- sum(inc[a >= 15]) / tot
  sel <- a >= 2 & a < 10
  pfpr <- sum(m[sel, "pos_M"]) / sum(m[sel, "prop"])
  list(EIR = EIR, pfpr = pfpr, u5 = u5, mid = mid, ad = ad,
       inc_py = tot * 365, m = m)
}

EIRs <- c(1, 3, 10, 30, 100, 300)
res  <- lapply(EIRs, summarise_eq)

cat("\n================ AGE DISTRIBUTION OF CLINICAL CASES vs EIR (ft = 0) ================\n")
tab <- data.frame(
  EIR          = EIRs,
  PfPR_2_10    = sprintf("%.0f%%", 100 * sapply(res, `[[`, "pfpr")),
  episodes_py  = round(sapply(res, `[[`, "inc_py"), 2),
  pct_U5       = sprintf("%.0f%%", 100 * sapply(res, `[[`, "u5")),
  pct_5_15     = sprintf("%.0f%%", 100 * sapply(res, `[[`, "mid")),
  pct_15plus   = sprintf("%.0f%%", 100 * sapply(res, `[[`, "ad"))
)
print(tab, row.names = FALSE)
cat("\nExpect: PfPR and the under-5 share RISE with EIR; the 15+ share FALLS.\n",
    "        These values should match the explainer's readouts (same solver, grid, ft).\n", sep = "")

# ----------------------------------------------------------------------------
# Figure: (A) age-group shares vs EIR, (B) case-by-age density at low/mid/high.
# ----------------------------------------------------------------------------
grid_EIR <- exp(seq(log(0.5), log(512), length.out = 60))
shares   <- lapply(grid_EIR, summarise_eq)
df_share <- rbind(
  data.frame(EIR = grid_EIR, share = sapply(shares, `[[`, "u5"),  band = "under 5"),
  data.frame(EIR = grid_EIR, share = sapply(shares, `[[`, "mid"), band = "5-15"),
  data.frame(EIR = grid_EIR, share = sapply(shares, `[[`, "ad"),  band = "15+")
)
df_share$band <- factor(df_share$band, levels = c("under 5", "5-15", "15+"))

gA <- ggplot(df_share, aes(EIR, 100 * share, colour = band)) +
  geom_line(linewidth = 1) +
  scale_x_log10() +
  scale_colour_manual(values = c("under 5" = "#e8804b", "5-15" = "#2f8f7e", "15+" = "#3d6fb4"),
                      name = NULL) +
  labs(title = "Age distribution of clinical malaria vs transmission intensity",
       subtitle = "Griffin equilibrium (malariaEquilibrium), ft = 0; the toy reproduces this",
       x = "EIR (infectious bites / person / year, log scale)", y = "% of clinical cases") +
  theme_bw(base_size = 11) + theme(legend.position = "top")

# density-by-age at three intensities (scaled to own peak, as the page draws it)
dens_at <- function(EIR, lab) {
  m <- summarise_eq(EIR)$m
  a <- m[, "age"]; d <- m[, "inc"]; d <- d / max(d[a <= 60])
  data.frame(age = a, dens = d, level = lab)
}
df_age <- rbind(dens_at(2, "low (EIR 2)"), dens_at(20, "moderate (EIR 20)"),
                dens_at(200, "high (EIR 200)"))
df_age$level <- factor(df_age$level, levels = c("low (EIR 2)", "moderate (EIR 20)", "high (EIR 200)"))

gB <- ggplot(subset(df_age, age <= 60), aes(age, dens, colour = level)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("low (EIR 2)" = "#9ec2e8",
                                 "moderate (EIR 20)" = "#7b5bd6",
                                 "high (EIR 200)" = "#3b1f7a"), name = NULL) +
  labs(title = "Where the cases are, by age",
       subtitle = "Each curve scaled to its own peak (as the explainer draws it)",
       x = "age (years)", y = "relative cases") +
  theme_bw(base_size = 11) + theme(legend.position = "top")

dir.create("figures", showWarnings = FALSE)
png_path <- file.path(getwd(), "figures", "age_distribution.png")
if (requireNamespace("gridExtra", quietly = TRUE)) {
  g <- gridExtra::arrangeGrob(gA, gB, ncol = 2)
  ggsave(png_path, g, width = 11, height = 4.5, dpi = 120)
} else {
  ggsave(png_path, gA, width = 7, height = 4.5, dpi = 120)   # fallback: shares panel only
}
cat("\nSaved figure to: ", png_path, "\n", sep = "")
