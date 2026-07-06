# Expression-matched bulk scoring

This folder contains the unified bookdown source for the expression-matched
bulk scoring tutorial.

View the knitted report at:

```text
report/index.html
```

Render it again from the repository root with:

```bash
Rscript -e 'bookdown::render_book("report")'
```

Run the analysis workflow that generates the figures first:

```bash
Rscript scripts/01_run_GSE122380_module_score_tutorial.R
```

The scoring functions live in `functions/` and are loaded with
`source("load_functions.R")` from the repository root.
