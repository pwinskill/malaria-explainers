# =============================================================================
# validate_rebounds.R
#
# Static R/ggplot reproduction of the interactive toy model in rebounds.html,
# as an independent cross-check of the JavaScript. It re-implements the same
# "rebound" dynamics
#
#        total protection(t)  =  acquired immunity I(t)  +  intervention J(t)
#
# where a transmission-reducing intervention suppresses exposure (so immunity
# wanes toward a lower setpoint), then is withdrawn (so its protection decays),
# and total protection can dip BELOW the pre-intervention baseline = a rebound.
#
# All parameters are [ASSUMED] illustrative choices (a toy model), sense-checked
# against the Griffin et al. / malariasimulation transmission model for
# structure and timescales. Not data-fitted.
#
# Run from the repository root:
#   Rscript validate_rebounds.R
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
## 1. PARAMETERS  (identical to the JavaScript in rebounds.html)
## ===========================================================================
k_half <- 28        # [ASSUMED] half-saturation of immunity-exposure curve
t_min  <- -3        # [ASSUMED] time axis start (years)
t_max  <- 20        # [ASSUMED] time axis end (years)
dt     <- 0.02      # [ASSUMED] Euler step (years)

## --- Colours (match the web page) -------------------------------------------
col_accent <- "#2f8f7e"  # immunity
col_warm   <- "#e8804b"  # intervention
col_danger <- "#d1495b"  # rebound
col_ink    <- "#26303b"
col_muted  <- "#6b7785"

## ===========================================================================
## 2. MODEL FUNCTIONS  (identical logic to the JavaScript)
## ===========================================================================

## Saturating (Hill) exposure -> immunity setpoint: g(x) = x^2 / (x^2 + K^2)
g_imm <- function(x) {
  x2 <- x * x
  x2 / (x2 + k_half * k_half)
}

## Immunity relaxes toward its setpoint at a rate that rises mildly with baseline
## transmission (faster rebuild where transmission is higher).
rate_of <- function(Tr) 0.08 + 0.20 * (Tr / 100)

## Speed-of-loss slider (0..100) -> intervention decay half-life (years)
##   0   = gradual waning -> ~3.0 yr ; 100 = abrupt stop -> ~0.17 yr (~2 mo)
loss_half_life <- function(v) 3.0 * (0.17 / 3.0) ^ (v / 100)

## ---- simulate one scenario -------------------------------------------------
## Tr    baseline transmission (5..100)
## S     control strength (0..0.95, fraction reduction in exposure while on)
## D     control duration (years)
## lossV speed-of-loss slider value (0..100)
simulate <- function(Tr, S, D, lossV) {
  B       <- g_imm(Tr)                  # baseline immunity floor
  r       <- rate_of(Tr)
  hl      <- loss_half_life(lossV)
  k_decay <- log(2) / hl

  ts <- seq(t_min, t_max, by = dt)
  I  <- numeric(length(ts))
  J  <- numeric(length(ts))
  I[1] <- B                             # start at baseline equilibrium

  for (i in seq_along(ts)) {
    t <- ts[i]
    Jt <- if (t >= 0 && t < D) S
          else if (t >= D)     S * exp(-k_decay * (t - D))
          else                 0
    J[i] <- Jt
    x   <- Tr * (1 - Jt)                # suppressed exposure
    Ieq <- g_imm(x)
    if (i > 1) I[i] <- I[i - 1] + r * (Ieq - I[i - 1]) * dt
  }
  I <- pmax(0, I)
  tot <- I + J

  ## rebound metric: how far total protection dips below baseline after withdrawal
  post     <- ts >= D
  min_tot  <- min(tot[post])
  depth    <- max(0, B - min_tot)

  list(
    df = data.frame(t = ts, I = I, J = J, tot = tot),
    B = B, depth = depth,
    rel_depth = if (B > 0) depth / B else 0,
    Tr = Tr, S = S, D = D, lossV = lossV
  )
}

## ===========================================================================
## 3. NUMERIC CROSS-CHECK  (compare these against the web-page readouts)
## ===========================================================================
scenarios <- list(
  "default            (Tr=40, S=0.70, D=6,  abrupt)"   = c(40, 0.70, 6,  75),
  "strong+long+abrupt (Tr=40, S=0.90, D=10, abrupt)"   = c(40, 0.90, 10, 95),
  "gradual loss       (Tr=40, S=0.90, D=10, gradual)"  = c(40, 0.90, 10, 10),
  "very high transm.  (Tr=95, S=0.90, D=10, abrupt)"   = c(95, 0.90, 10, 95),
  "very low transm.   (Tr=10, S=0.90, D=10, abrupt)"   = c(10, 0.90, 10, 95)
)

cat("\n=== Rebound toy model — numeric cross-check ===\n")
cat(sprintf("%-52s %8s %8s %8s\n", "scenario", "baseB", "depth", "rel.dep"))
sim_list <- lapply(scenarios, function(p) simulate(p[1], p[2], p[3], p[4]))
for (nm in names(sim_list)) {
  s <- sim_list[[nm]]
  cat(sprintf("%-52s %8.3f %8.3f %7.1f%%\n",
              nm, s$B, s$depth, 100 * s$rel_depth))
}

## --- Baseline-transmission check: rebound DIP vs baseline transmission -------
## Per Ghani et al. 2009, rebound risk grows with transmission. The protection
## dip should rise with baseline transmission and be largest at high transmission.
trange <- seq(5, 100, by = 1)
dip_curve <- sapply(trange, function(Tr) simulate(Tr, 0.90, 10, 95)$depth)
peak_Tr <- trange[which.max(dip_curve)]
cat(sprintf("\nBaseline check (S=0.90, D=10, abrupt): dip at Tr=5 is %.3f, ",
            dip_curve[1]))
cat(sprintf("at Tr=100 is %.3f, max at Tr=%d.\n",
            dip_curve[length(dip_curve)], peak_Tr))
cat(if (peak_Tr >= 75 && dip_curve[length(dip_curve)] > 3 * dip_curve[1])
      "  -> rebound dip largest at high baseline transmission, as expected.\n"
    else "  -> WARNING: expected the largest dip at high baseline transmission.\n")

## ===========================================================================
## 4. PLOTS
## ===========================================================================

## --- (a) the headline scenario: stacked protection-over-time ----------------
plot_scenario <- function(s, title) {
  df  <- s$df
  long <- df %>%
    transmute(t,
              Immunity     = I,
              Intervention = J) %>%
    pivot_longer(-t, names_to = "component", values_to = "v")
  long$component <- factor(long$component, levels = c("Intervention", "Immunity"))

  ## rebound (dip) ribbon: total below baseline, after withdrawal
  dip <- df %>%
    filter(t >= s$D, tot < s$B) %>%
    mutate(lo = tot, hi = s$B)

  ggplot() +
    geom_area(data = long, aes(t, v, fill = component),
              position = "stack", alpha = 0.42, colour = NA) +
    { if (nrow(dip) > 0)
        geom_ribbon(data = dip, aes(t, ymin = lo, ymax = hi),
                    fill = col_danger, alpha = 0.5) } +
    geom_hline(yintercept = s$B, linetype = "dashed", colour = col_muted) +
    geom_vline(xintercept = c(0, s$D), linetype = "dotted",
               colour = col_ink, alpha = 0.5) +
    geom_line(data = df, aes(t, tot), colour = col_ink, linewidth = 1) +
    annotate("text", x = 0,    y = max(df$tot) * 1.05, label = "control on",
             hjust = 0, size = 3, colour = col_ink) +
    annotate("text", x = s$D,  y = max(df$tot) * 1.05, label = "control off",
             hjust = 0, size = 3, colour = col_ink) +
    scale_fill_manual(values = c(Immunity = col_accent, Intervention = col_warm),
                      breaks = c("Immunity", "Intervention")) +
    labs(title = title,
         subtitle = sprintf("protection dips %.0f%% of baseline below it after withdrawal",
                            100 * s$rel_depth),
         x = "time (years, control starts at 0)", y = "total protection",
         fill = NULL) +
    coord_cartesian(xlim = c(t_min, t_max), ylim = c(0, max(df$tot) * 1.12)) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

p_main    <- plot_scenario(sim_list[[2]], "Strong, long, abrupt withdrawal — deep rebound")
p_gradual <- plot_scenario(sim_list[[3]], "Same programme, gradual loss — rebound softened")

## --- (b) rebound dip vs baseline transmission -------------------------------
dip_df <- data.frame(Tr = trange, dip = dip_curve)
p_hump <- ggplot(dip_df, aes(Tr, dip)) +
  geom_line(colour = col_danger, linewidth = 1) +
  labs(title = "Rebound deepens with baseline transmission (toy model)",
       subtitle = "strong, long, abruptly-withdrawn programme — dip of total protection below baseline",
       x = "baseline transmission", y = "rebound dip (protection units)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

## --- combine & save ---------------------------------------------------------
combined <- (p_main | p_gradual) / p_hump +
  plot_annotation(
    title = "Rebound toy model — R cross-check of rebounds.html",
    caption = "Illustrative toy model. Not for decision making."
  )

ggsave("rebounds_validation.png", combined,
       width = 11, height = 8, dpi = 130, bg = "white")
cat("\nSaved rebounds_validation.png\n")
