# Notum vGEP

R Shiny application for visualizing and exploring spatial gene expression data in the *Drosophila melanogaster* pupal notum epithelium.

Developed by the [Bellaiche Lab](https://science.curie.fr/recherche/biologie-interactive-des-tumeurs-immunologie-environnement/equipe-bellaiche/) (Institut Curie).

**Live app:** https://ybellaichelab.shinyapps.io/NotumvGEP/

The data come from [vgep](https://github.com/BellaicheLab/vgep); the paper's analyses are in [morphodomains](https://github.com/BellaicheLab/morphodomains).

## Features

- **Atlas** — Browse protein expression patterns (24 genes)
- **Validation** — Compare FISH and protein validation data (7 protein + 13 RNA smiFISH patterns, plus 3 additional protein patterns: Tollo, yrt, stan)
- **vGEP** — Explore predicted spatial expression (4,312 genes)
- **Correlate** — Find genes correlated with a gene of interest
- **Draw** — Lasso-select a spatial region to discover enriched genes
- **Compare** — Side-by-side comparison across datasets
- **Overlay** — Two-color (red/green) overlay of gene pairs

## Running locally

```r
install.packages(c("shiny", "shinyjs", "shinyWidgets", "shinythemes", "ggplot2", "ggiraph", "DT"))
shiny::runApp("app.R")
```

## Reference

Please acknowledge the original publication when using this application or its data:
Leroy, Van Leen *et al.*, Cell Reports (2026). Code and data: [10.5281/zenodo.19633285](https://doi.org/10.5281/zenodo.19633285)