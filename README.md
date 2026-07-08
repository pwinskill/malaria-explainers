# Malaria modelling explainers

A small, self-contained collection of interactive explainers on how malaria transmission behaves
and how control interventions reshape it. Each explainer pairs a plain-language walkthrough with a
"toy" model you can play with, and is cross-checked in R against published model structures.

Everything is plain HTML, CSS and vanilla JavaScript (canvas + SVG), with no build step or
dependencies. `index.html` is a hub linking to the explainers; shared styling lives in `styles.css`.
The R cross-check scripts live in `validation/`; their input data (`data/`), reference papers
(`references/`), and generated figures (`figures/`) are kept out of version control.

**Live page:** https://pwinskill.github.io/malaria-explainers/

Visits to the live page are counted anonymously with [GoatCounter](https://www.goatcounter.com/) — a cookieless, privacy-friendly analytics service that collects no personal data.

## Explainers

1. **Non-linearities in malaria epidemiology and control** ([`nonlinearities.html`](nonlinearities.html)).
   Why a modest change in mosquito survival or bed-net coverage can produce a much larger (or
   surprisingly small) change in prevalence and clinical cases. Traces a bed-net effect through the
   coverage → EIR → prevalence → incidence cascade. The coverage → EIR step is a mechanistic
   experimental-hut net model; the EIR → PfPR and EIR → clinical-incidence relationships both come
   from the Griffin-model equilibrium ([malariaEquilibrium](https://github.com/mrc-ide/malariaEquilibrium)),
   the same solver used by the age-distribution explainer, with overall incidence split into age-group
   contributions. Cross-check: `validation/validate_cascade.R`.

2. **Rebounds: why protection can dip below baseline** ([`rebounds.html`](rebounds.html)).
   When transmission-reducing control is withdrawn, naturally-acquired immunity has waned underneath
   it, so total protection can fall *below* where it started. A triangle of three drivers (baseline
   transmission, strength and duration of control, speed of withdrawal and loss of protection) plus
   a toy protection-over-time model. Structure and timescales sense-checked against the Griffin
   et al. / [malariasimulation](https://github.com/mrc-ide/malariasimulation) model.
   Cross-check: `validation/validate_rebounds.R`.

3. **Seasonality and the timing of control** ([`seasonality.html`](seasonality.html)).
   Where transmission is seasonal, cases arrive in a short window and each dose of seasonal malaria
   chemoprevention (SMC) protects for only a few weeks, so a programme has to tile its rounds across
   the season. A toy model of childhood cases through the year, with sliders for seasonality, peak
   timing, the number, spacing and start of the SMC rounds, coverage per round, and whether the same,
   random or different children are reached across rounds (inter-round correlation), plus an
   "optimise timing" button that searches for the round schedule averting the most cases. Shows how
   impact depends on when the rounds land and on who is reached. The drug-prophylaxis curve is taken directly
   from [malariasimulation](https://github.com/mrc-ide/malariasimulation), and the seasonal SMC and
   correlation behaviour is cross-checked against it. Cross-check: `validation/validate_seasonality.R`.

4. **The age distribution of malaria burden** ([`age-distribution.html`](age-distribution.html)).
   Malaria is a disease of young children only where transmission is intense: as transmission rises,
   acquired immunity builds faster and clinical cases concentrate in the very young, while older
   children and adults are increasingly protected. A live, in-browser port of the Griffin-model
   equilibrium ([malariaEquilibrium](https://github.com/mrc-ide/malariaEquilibrium)) with a single
   transmission-intensity slider. The chart shows clinical incidence per person by age (rising into
   adulthood at low transmission, concentrating in infancy at high transmission) alongside the
   case distribution, with the under-5 / 5–15 / 15+ shares beside it. Cross-check:
   `validation/validate_age_distribution.R`.

5. **Non-malaria fevers and inflated malaria case counts** ([`nmf.html`](nmf.html)).
   Routine surveillance records a malaria case whenever someone with a fever tests
   positive, but where malaria is endemic many infections are asymptomatic, so a
   non-malaria fever that coincides with one is recorded as malaria all the same.
   Because non-malaria fevers occur at a roughly constant rate while true malaria is
   seasonal, this over-counting is differential: it is worst in the low season, when
   true clinical malaria has collapsed but the asymptomatic reservoir has not. A toy
   model of true versus routine-recorded childhood cases through the year, with a
   slider for the non-malaria fever rate alongside the transmission setting. The
   coincidence structure follows the non-malarial-fever process in
   [malariasimulation](https://github.com/mrc-ide/malariasimulation/pull/372), and
   the over-count fractions are anchored to the estimates of Dalrymple et al. (2017).
   Cross-check: `validation/validate_nmf.R`.

## R cross-checks

The `validation/validate_*.R` scripts independently reproduce each explainer's model in R/ggplot
and print numeric checks, as a guard against bugs in the JavaScript. Run from the repository root,
e.g.:

```
Rscript validation/validate_rebounds.R
```

> The numbers are illustrative. They show the shape of the relationships, not real-world impact,
> and are not to be used for decision making.
