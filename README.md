# Phylogenetic Tree: *Yersinia pestis* Metalloproteome

Circular phylogenetic tree figure for the metalloproteomic study of
*Yersinia pestis* CO92 (UniProt reference proteome UP000000815).

**Author:** Divyam Urja

## Contents

| Path | Contents |
| --- | --- |
| `data/raw/` | Original MEGA exports, unmodified |
| `data/alignment/` | Cleaned protein alignment, 92 sequences x 1412 columns |
| `data/tree/` | Newick tree, 92 tips, bootstrap support as node labels |
| `data/metadata/` | Per-tip table for annotation |
| `data/annotation/` | eggNOG-mapper functional annotation |
| `reference/` | Style reference figures |
| `R/` | Pipeline scripts |
| `figures/output/` | Generated figures |

## Dataset

92 protein-metal pairs across 78 distinct proteins. Fourteen proteins bind
more than one metal and appear once per metal; these are recoded with a
`_copyN` suffix because parentheses are structural characters in Newick.

The tree relates different proteins of a single organism, so it is a
clustering of sequence similarity across the metalloproteome rather than a
species phylogeny.

## Build

See [`IMPLEMENTATION_PLAN.md`](IMPLEMENTATION_PLAN.md).

```bash
python3 R/00_prepare_data.py     # regenerates data/ from data/raw/
Rscript R/01_annotate.R          # join eggNOG + COG categories
Rscript R/03_tree_figure.R       # Figure 1, circular clustergram
Rscript R/04_function_figure.R   # Figure 2, functional composition
Rscript R/05_verify.R            # build gate, fails on any violation
```

Requires `ape`, `ggplot2`, `dplyr`, `readr` and `colorspace`. No Bioconductor
dependency.

## Figure standards

- Colours: Paul Tol `light` qualitative scheme, pastel and colour-blind safe
  under deuteranopia, protanopia and tritanopia
- Tip labels verified non-overlapping by construction at the specified
  canvas size, and asserted in `R/05_verify.R`
- Export: vector PDF plus 600 dpi PNG
- Figure 1 is a cladogram: branch lengths are not represented, because they
  span three orders of magnitude and collapse the topology when drawn to scale

## Figures

| File | Contents |
| --- | --- |
| `figures/output/fig1_tree.pdf` | Figure 1, circular clustergram, vector |
| `figures/output/fig1_tree_600dpi.png` | Figure 1, 600 dpi raster |
| `figures/output/fig2_function.pdf` | Figure 2, functional composition, vector |
| `figures/output/fig2_function.png` | Figure 2, raster |

Regenerate either with the commands above. `R/05_verify.R` gates the build and
fails on tip-set mismatch, missing annotation, an unsafe palette or any risk of
overlapping tip labels.
