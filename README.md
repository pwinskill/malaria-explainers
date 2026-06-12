# Non-linear impacts of bed nets

A small, self-contained interactive explainer on why insecticide-treated bed nets (ITNs)
have a *non-linear* effect on malaria transmission.

It walks through the Ross–Macdonald model of the basic reproduction number (R₀), shows why
mosquito survival matters so much (the parasite must survive the ~10-day extrinsic incubation
period, so survival enters as an exponential term), and lets you drag an ITN-coverage slider to
watch R₀ fall.

**Live page:** https://pwinskill.github.io/nonlinear-bednets/

Everything lives in a single `index.html` — plain HTML, CSS and vanilla JavaScript (canvas charts),
with no build step or dependencies.

> The numbers are illustrative — they show the shape of the relationship, not real-world impact.
> Not to be used for decision making.
