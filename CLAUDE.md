# CLAUDE.md

Project overview, file layout and the R cross-check workflow live in [README.md](README.md).
This file only notes things specific to working on the repo with Claude.

## Conventions

- **No build step, no dependencies.** Plain HTML, CSS and vanilla JavaScript
  (canvas + SVG). Don't add a bundler, framework or npm packages.
- **Each explainer is one self-contained `.html` file** with its model logic in
  an inline `<script>`. The only shared asset is `styles.css`. `index.html` is
  the hub and carries no script.
- **Numbers are illustrative.** They show the shape of relationships, not
  real-world impact. Never frame them as fit for decision-making.

## Verifying changes

- There is no test/lint runner. Verify the web side by opening the `.html` in a
  browser.
- Every explainer is mirrored by an R cross-check in `validation/validate_*.R`.
  **When you change a model in an explainer, keep its `validate_*.R` in sync (and
  vice versa)** — the two are meant to reproduce the same numbers.
- Run the R scripts from the repository root (they read `data/`), e.g.
  `Rscript validation/validate_rebounds.R`. R toolchain paths (Windows arm64) are
  in the user's global `~/.claude/CLAUDE.md`.
