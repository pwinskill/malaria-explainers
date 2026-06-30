# Malaria modelling explainers

A small, self-contained collection of interactive explainers on how malaria transmission behaves
and how control interventions reshape it. Each explainer pairs a plain-language walkthrough with a
"toy" model you can play with, and is cross-checked in R against published model structures.

Everything is plain HTML, CSS and vanilla JavaScript (canvas + SVG), with no build step or
dependencies. `index.html` is a hub linking to the explainers; shared styling lives in `styles.css`.

**Live page:** https://pwinskill.github.io/malaria-explainers/

## Explainers

1. **Non-linearities in malaria epidemiology and control** ([`nonlinearities.html`](nonlinearities.html)).
   Why a modest change in mosquito survival or bed-net coverage can produce a much larger (or
   surprisingly small) change in prevalence and clinical cases. Traces a bed-net effect through the
   EIR → prevalence → incidence cascade. Cross-check: `validate_cascade.R`.

2. **Rebounds: why protection can dip below baseline** ([`rebounds.html`](rebounds.html)).
   When transmission-reducing control is withdrawn, naturally-acquired immunity has waned underneath
   it, so total protection can fall *below* where it started. A triangle of three drivers (baseline
   transmission, strength and duration of control, speed of withdrawal and loss of protection) plus
   a toy protection-over-time model. Structure and timescales sense-checked against the Griffin
   et al. / [malariasimulation](https://github.com/mrc-ide/malariasimulation) model.
   Cross-check: `validate_rebounds.R`.

## R cross-checks

The `validate_*.R` scripts independently reproduce each explainer's model in R/ggplot and print
numeric checks, as a guard against bugs in the JavaScript. Run from the repository root, e.g.:

```
Rscript validate_rebounds.R
```

> The numbers are illustrative. They show the shape of the relationships, not real-world impact,
> and are not to be used for decision making.
