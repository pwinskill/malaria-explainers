# =============================================================================
# validate_attribution.R
#
# Cross-check the "attributing impact of malaria interventions"
# explainer (attribution.html) against the canonical Griffin-model equilibrium
# (mrc-ide/malariaEquilibrium).
#
# The explainer combines three unnamed controls with fixed, illustrative effects:
#   - A  reduces TRANSMISSION by a fraction t_A (moving the setting along the
#        saturating clinical-incidence curve f(EIR)) AND gives DIRECT protection
#        d_A (a multiplicative reduction of the residual clinical cases);
#   - B  gives DIRECT protection d_B;
#   - C  gives DIRECT protection d_C (< d_B): B and C are both direct-only, B the stronger.
# For a deployed set S:  cases(S) = N * f(EIR_S) * prod(1 - d_j over direct effects in S),
# with EIR_S = EIR*(1 - t_A) if A is in S and EIR otherwise.
#
# All three interventions are given the SAME cost, so cost per case averted
# depends only on the cases each is credited with.
#
# f(EIR) is all-age clinical incidence per person-year from human_equilibrium(),
# the same solver the explainer ports; this script matches its grid (seq(0,80,0.1))
# and ft = 0, so the numbers should equal the page's readouts to rounding.
#
# It checks the claims the page makes:
#   1. The COMBINED cases averted is identical for every order of introduction.
#   2. For any order, the per-control increments sum to that combined total.
#   3. The fair share (Shapley value) also sums to the combined total.
#   4. Going first is worth more than going last, for every control.
#   5. B and C are both direct-only, with B stronger than C (more credit assessed
#      alone and by fair share).
#   6. At equal cost, a control is dearer per case averted added last than assessed
#      alone (the per-control backdrop effect).
#
# Run (from the repository root, so figures/ resolves):
#   & 'C:/Program Files/R-aarch64/R-4.5.2/bin/Rscript' validation/validate_attribution.R
# =============================================================================

.libPaths('C:/Users/pwinskil/Documents/r_packages_arm64')
suppressMessages({
  library(malariaEquilibrium)
  library(ggplot2)
})

## ===========================================================================
## 1. PARAMETERS  (identical to the defaults in attribution.html)
## ===========================================================================
N     <- 100000                       # [ASSUMED] illustrative population
ft    <- 0                            # [MODEL]   treatment coverage off, matching the toy
age   <- seq(0, 80, 0.1)              # [MODEL]   uniform grid, matching the JS port

# the transmission slider default is 60; the page maps it to EIR on a log scale
# EIR = EIR_MIN * (EIR_MAX / EIR_MIN)^(v/100), with EIR_MIN = 1, EIR_MAX = 256
EIR0  <- 1 * (256 / 1) ^ (60 / 100)   # [MODEL]   ~= 27.86 infectious bites / person / year

t_A   <- 0.50                         # [ASSUMED] A: fractional reduction in transmission
d     <- c(A = 0.30, B = 0.50, C = 0.40)  # [ASSUMED] direct protective effects (B stronger than C)
cost_pp <- 1                          # [ASSUMED] equal cost per person for every intervention ($)

# web-page colours (so the figure matches the explainer)
col_A <- "#2f8f7e"; col_B <- "#b0501f"; col_C <- "#7b5bd6"   # match the page's A/B/C label hues
col_ink <- "#26303b"; col_muted <- "#6b7785"

## ===========================================================================
## 2. MODEL  (identical logic to the JavaScript)
## ===========================================================================
p <- load_parameter_set()

# f(EIR): all-age clinical incidence per person-year, the burden map
f_eir <- function(EIR) {
  m <- human_equilibrium(EIR = EIR, ft = ft, p = p, age = age)$states
  sum(m[, "inc"]) * 365
}
f0 <- f_eir(EIR0)                     # A not deployed (baseline transmission)
fA <- f_eir(EIR0 * (1 - t_A))         # A deployed (transmission reduced)

ints <- c("A", "B", "C")

# cases remaining for a deployed subset (named logical over A/B/C)
cases_of <- function(S) {
  E  <- if (S["A"]) fA else f0                                 # only A changes transmission
  df <- prod(ifelse(S[ints], 1 - d[ints], 1))                 # direct effects (A's own included when present)
  N * E * df
}
empty      <- c(A = FALSE, B = FALSE, C = FALSE)
base_cases <- cases_of(empty)

# cases averted credited to each control as it is added along one order
increments_for <- function(order) {
  S <- empty; out <- setNames(numeric(3), ints); before <- cases_of(S)
  for (i in order) {
    S[i]  <- TRUE
    after <- cases_of(S)
    out[i] <- before - after
    before <- after
  }
  out
}

# the six orders in which three controls can be introduced
perms <- list(
  c("A","B","C"), c("A","C","B"), c("B","A","C"),
  c("B","C","A"), c("C","A","B"), c("C","B","A")
)

shapley <- Reduce(`+`, lapply(perms, increments_for)) / length(perms)   # fair share
alone   <- setNames(sapply(ints, function(i) { S <- empty; S[i] <- TRUE; base_cases - cases_of(S) }), ints)

# cost per case averted (equal cost per person for all three)
ce <- function(averted) ifelse(averted > 1e-6, (cost_pp * N) / averted, Inf)   # guard zero, mirroring the JS perCase()

## ===========================================================================
## 3. NUMERIC CROSS-CHECK  (compare against the web-page readouts)
## ===========================================================================
cat(sprintf("\nEIR0 = %.3f   f(EIR0) = %.4f   f(EIR0*(1-t_A)) = %.4f   baseline cases = %.0f\n",
            EIR0, f0, fA, base_cases))

total_all <- base_cases - cases_of(c(A = TRUE, B = TRUE, C = TRUE))
cat(sprintf("Combined cases averted (all three) = %.0f   (web: 58,884)\n\n", total_all))

cat("Cases averted credited by order (each row sums to the same combined total):\n")
cat(sprintf("  %-14s %8s %8s %8s %9s\n", "order", "A", "B", "C", "total"))
for (o in perms) {
  v <- increments_for(o)
  cat(sprintf("  %-14s %8.0f %8.0f %8.0f %9.0f\n", paste(o, collapse = ">"), v["A"], v["B"], v["C"], sum(v)))
}
inc_default <- increments_for(c("A", "B", "C"))
cat(sprintf("\nDefault order A>B>C credits  A %.0f (web 27,631), B %.0f (web 22,324), C %.0f (web 8,930)\n",
            inc_default["A"], inc_default["B"], inc_default["C"]))

cat(sprintf("\nFair share (Shapley): A %.0f, B %.0f, C %.0f   (sum %.0f)\n",
            shapley["A"], shapley["B"], shapley["C"], sum(shapley)))
cat(sprintf("Assessed alone:       A %.0f, B %.0f, C %.0f\n", alone["A"], alone["B"], alone["C"]))

cat(sprintf("\nCost per case averted, order A>B>C ($, equal cost): A $%.1f (web $3.6), B $%.1f (web $4.5), C $%.1f (web $11.2)\n",
            ce(inc_default["A"]), ce(inc_default["B"]), ce(inc_default["C"])))

## ---- assertions ------------------------------------------------------------
totals  <- sapply(perms, function(o) sum(increments_for(o)))
sums_ok <- sapply(perms, function(o) abs(sum(increments_for(o)) - total_all) < 1e-6)
first_of <- function(i) max(sapply(perms, function(o) increments_for(o)[i]))
last_of  <- function(i) min(sapply(perms, function(o) increments_for(o)[i]))

stopifnot(
  # 1. combined total identical for every order
  max(totals) - min(totals) < 1e-6,
  # 2. each order's increments sum to the combined total
  all(sums_ok),
  # 3. the fair share sums to the combined total
  abs(sum(shapley) - total_all) < 1e-6,
  # 4. going first is worth more than going last, for every control. first_of = max increment
  #    (occurs at first position) and last_of = min (last position) by submodularity; the STRICT
  #    inequality is non-vacuous (an additive / order-invariant model would fail it).
  all(sapply(ints, function(i) first_of(i) > last_of(i))),
  # 5. B and C are both direct-only, with B stronger than C: more credit assessed alone and by fair share
  alone["B"] > alone["C"],
  shapley["B"] > shapley["C"],
  # 6. at equal cost, a control is dearer per case averted added last (min credit) than first (max credit)
  ce(last_of("C")) > ce(first_of("C")),
  # cross-implementation: combined total matches the web readout to ~1%
  abs(total_all - 58884) < 0.01 * 58884,
  # per-control web values (guard JS<->R per-control drift, not just the total)
  round(inc_default["A"]) == 27631, round(inc_default["B"]) == 22324, round(inc_default["C"]) == 8930,
  abs(ce(inc_default["A"]) - 3.6) < 0.05, abs(ce(inc_default["B"]) - 4.5) < 0.05, abs(ce(inc_default["C"]) - 11.2) < 0.05
)
cat("\nAll assertions passed.\n\n")

## ===========================================================================
## 4. PLOTS
## ===========================================================================
pal <- c(A = col_A, B = col_B, C = col_C)
ord <- c("A", "B", "C")
inc0 <- increments_for(ord)
cumv <- cumsum(inc0[ord])

# Panel A: cumulative cases-averted staircase for the default order
dfA <- data.frame(step = factor(ord, levels = ord), x = seq_along(ord),
                  ymin = c(0, head(cumv, -1)), ymax = cumv, credit = inc0[ord])
gA <- ggplot(dfA) +
  geom_rect(aes(xmin = x - 0.34, xmax = x + 0.34, ymin = ymin, ymax = ymax, fill = step)) +
  geom_hline(yintercept = total_all, linetype = "dashed", colour = col_ink) +
  geom_text(aes(x = x, y = (ymin + ymax) / 2, label = format(round(credit), big.mark = ",")),
            colour = "white", fontface = "bold", size = 3.4) +
  annotate("text", x = 0.55, y = total_all, label = "combined total (fixed)",
           hjust = 0, vjust = -0.6, colour = col_muted, size = 3) +
  scale_x_continuous(breaks = dfA$x, labels = paste0(ord, "\n(", c("1st","2nd","3rd"), ")")) +
  scale_fill_manual(values = pal, guide = "none") +
  labs(title = "Cases averted credited as each control is added",
       subtitle = "Default order A > B > C; the top of the staircase is the fixed combined total",
       x = NULL, y = "cases averted / 100,000") +
  theme_bw(base_size = 11)

# Panel B: cost per case averted for this order (equal cost)
dfB <- data.frame(intervention = factor(ord, levels = ord), value = ce(inc0[ord]),
                  pos = paste0("(", c("1st","2nd","3rd"), ")"))
gB <- ggplot(dfB, aes(intervention, value, fill = intervention)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("$%.1f", value)), vjust = -0.5, fontface = "bold", size = 3.4) +
  scale_x_discrete(labels = paste0(ord, "\n", dfB$pos)) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Cost per case averted, order A > B > C",
       subtitle = "Equal cost, so the control credited with fewer cases costs more per case averted",
       x = NULL, y = "cost per case averted ($)") +
  theme_bw(base_size = 11)

dir.create("figures", showWarnings = FALSE)
png_path <- file.path(getwd(), "figures", "attribution_validation.png")
if (requireNamespace("patchwork", quietly = TRUE)) {
  ggsave(png_path, patchwork::wrap_plots(gA, gB, ncol = 2), width = 11, height = 4.6, dpi = 120, bg = "white")
} else {
  ggsave(png_path, gA, width = 6, height = 4.6, dpi = 120, bg = "white")   # fallback: staircase only
}
cat("Saved figure to: ", png_path, "\n", sep = "")
