# Phylogenetic Tree Figure: Implementation Plan

**Project:** Metalloproteome of *Yersinia pestis* CO92 (UP000000815)
**Author:** Divyam Urja
**Deliverable:** Publication-quality circular phylogenetic tree figure, built in R
**Status:** Built. Both figures generate and pass the verification gate.

---

## 1. What the data actually is

Verified against the raw files, not assumed:

| Property | Value |
| --- | --- |
| Alignment | 92 sequences x 1412 columns, protein, MEGA export |
| Tree | 92 tips, 89 internal nodes carrying bootstrap support |
| Distinct proteins | 78 |
| Duplicated tips | 14 (proteins binding more than one metal; MEGA wrote them as `(2)`/`(3)`) |
| Reviewed (Swiss-Prot) | 36 of 92 tips |
| eggNOG annotation | 92 rows, covering **all 78 proteins, 100%** |
| Gene names available | 92 of 92, no gaps |

### Two defects found in the raw export, both fixed

1. `Divyam_tree.meg` contained a sequence named `Sequence_1` that is **entirely gaps**. It is an artefact of the MEGA export and is absent from the tree. It has been dropped, which is what brings the alignment from 93 to 92 and makes it agree with the tree.
2. One record carried its full UniProt description inside the sequence name
   (`sp|Q8ZJ91|RL361_YERPE_"Large ribosomal subunit protein bL36A OS=..."`), so it did not match its own tip in the tree. The name is now normalised.

After both fixes, the tip sets of the tree, the alignment and the metadata table are **identical, all 92**. This is asserted in the prep script and must stay true.

### Alignment quality, stated honestly

The alignment is sparse: mean column occupancy 26.6%, and 756 of 1412 columns are more than 90% gaps. That is expected, because these are 78 different proteins rather than orthologs of one protein.

It is **not**, however, the pathological "block diagonal" case where sequences share no columns. Measured: median pairwise overlap is 246 ungapped columns, and only 2.2% of sequence pairs share no columns at all. The alignment carries real signal.

What follows from this, and belongs in the methods section:

> This tree relates different proteins of one organism, not one protein across organisms. It should be described as a clustering of sequence similarity across the metalloproteome, not as a species phylogeny, and no divergence times or ancestry are inferred from it.

### Bootstrap support is weak, and the figure must not hide it

| Support band | Nodes |
| --- | --- |
| 90-100 | 14 |
| 70-89 | 1 |
| 50-69 | 1 |
| below 50 | **73** |

Only 15 of 89 internal nodes reach 70%. **The backbone is essentially unresolved.** Three consequences that drive the design:

- Show support **only at nodes >= 70**, and say in the caption that unlabelled nodes fall below 70. Printing 89 numbers, most of them single digits, would be both ugly and misleading.
- Do **not** cut the tree into clades algorithmically to drive the colours. With an unsupported backbone, any such grouping is an artefact of an arbitrary cut height, and a reviewer can say so. This supersedes the earlier "let the tree define the clades" decision.
- Colour by **eggNOG COG functional category** instead. It is externally derived, biologically meaningful, independent of the weak topology, and it is already complete for every protein.

---

## 2. Does the eggNOG data need its own tree? No.

The eggNOG output annotates the **same 78 proteins** already in the tree. A second tree would require a second alignment, and there is only one. A "tree of COG categories" is not a phylogeny, it is a classification.

So the functional data becomes an **annotation layer on the one tree**, which is exactly what the outer rings in reference image 1 are for:

- **Figure 1** (main): the tree, branches coloured by COG functional category, with outer rings carrying category and, later, metal.
- **Figure 2** (companion, not a tree): COG category composition. A horizontal bar chart of proteins per category, optionally split by metal once that mapping exists. This is the honest "functionality" figure and it is far more readable than a second tree.

---

## 3. Repository layout

```
data/
  raw/          Divyam_tree.meg, mega_tree_raw.nwk, tree_trial_2.mas   (untouched originals)
  alignment/    metalloproteome_aln.fasta      92 x 1412, cleaned
  tree/         mega_tree.nwk                  92 tips, support as 0-100 node labels
  metadata/     tips.csv                       one row per tip, ready to join
  annotation/   query.emapper.*                eggNOG-mapper output
reference/      ref1..ref4                     style references
R/              00_prepare_data.py + the R scripts below
figures/output/ generated figures (git-ignored except finals)
```

### Tip label recoding, important for downstream tools

MEGA wrote repeated accessions as `ACC_YERPE(2)`. Parentheses are **structural characters in Newick**, so leaving them in bare labels produces a file that silently mis-parses. They are recoded to `ACC_YERPE_copy2`. Keep this convention everywhere; `tips.csv` carries the `copy` number in its own column.

---

## 4. Pipeline

Tree inference is **not** repeated. MEGA did the inference and its bootstrap values are used as-is. R renders. The methods section says exactly that:

> The maximum-likelihood tree was inferred in MEGA with bootstrap resampling; the figure was rendered in R (ggtree) from the exported Newick file.

| Script | Does |
| --- | --- |
| `R/00_prepare_data.py` | Already written and run. Cleans the MEGA exports, drops the all-gap sequence, normalises names, converts support to percent, emits the three `data/` files, asserts the tip sets match. |
| `R/01_annotate.R` | Join eggNOG onto `tips.csv`. Resolve COG ID to functional category letter. Build `display_label`. Writes `data/metadata/tips_annotated.csv`. |
| `R/02_palette.R` | Define and **verify** the palette. Fails the build if it does not pass. |
| `R/03_tree_figure.R` | Figure 1. |
| `R/04_function_figure.R` | Figure 2. |
| `R/05_verify.R` | Overlap and colour checks. Fails the build on violation. |

### 4.1 `01_annotate.R`

```r
ann <- readr::read_tsv("data/annotation/query.emapper.annotations", comment = "##")
tips <- readr::read_csv("data/metadata/tips.csv")
# tips.csv already carries `copy`; strip the suffix to get the join key
tips$protein <- sub("_copy\\d+$", "", tips$tip_label)
tips <- dplyr::left_join(tips, ann, by = c("protein" = "#query"))
stopifnot(!any(is.na(tips$Preferred_name)))   # must be 92/92
```

**COG category resolution.** This emapper build puts the **COG ID** (`COG0086`) in the `COG_category` column, not the single-letter functional category. It is resolved against a COG definitions table committed at `data/annotation/cog_categories.tsv` (`COGID <tab> category <tab> name`, 4631 rows). 69 of the 70 COG IDs present resolve directly. `COG6021` (*bepA*) is absent and is assigned category **O** in `01_annotate.R`, with the reason stated inline: it is a periplasmic M48 zinc metalloprotease acting in outer-membrane protein quality control. Then:

- Join on COG ID to get the category letter(s).
- A COG may carry several letters (e.g. `EG`). Take the **first** as primary, and keep the full string in a `cog_all` column.
- 3 proteins have no COG and 1 is category `S` (function unknown). Map all 4 to a single **"Unassigned"** group and give it Tol's pale grey, `#DDDDDD`, which is reserved for exactly this.
- Collapse the 20-odd COG letters into **6 to 8 display groups**. With 78 proteins, 20 categories would give singleton colours and an unreadable legend. Group by the standard COG super-categories and merge thin ones:
  - Translation and ribosome (J)
  - Transcription and replication/repair (K, L)
  - Energy and central metabolism (C, G)
  - Amino acid, nucleotide and coenzyme metabolism (E, F, H)
  - Inorganic ion transport (P)
  - Cell wall, membrane and envelope (M, N, U)
  - Other / poorly characterised (R, S, and any remainder)
  - Unassigned

  Confirm the final grouping against the real counts after the join and merge any group that ends up below ~4 proteins. **Target 7 groups plus Unassigned.**

**Display label:** `paste0(Preferred_name, " (", accession, ")")`, giving `rpoC (Q8D1H3)`. For a `_copy2` tip append the metal once known, or a thin superscript, so duplicates are distinguishable rather than looking like a mistake.

### 4.2 `02_palette.R`, colours

**Paul Tol `light`** qualitative scheme. Pastel by design rather than by tinting a saturated palette, and the only soft scheme that survives all three CVD simulations at this many groups while keeping dark text readable on it.

```r
tol_light <- c("#77AADD","#EE8866","#EEDD88","#FFAABB",
               "#99DDFF","#44BB99","#BBCC33","#DDDDDD")
```

Source: Paul Tol, *Colour Schemes*, SRON/EPS/TN/09-002, figure 4.

Rules, each with a reason:

- **Reserve `#DDDDDD` for "Unassigned".** Tol intends the grey for exactly this.
- **Do not place `#BBCC33` (pear) next to `#AAAA00` (olive).** They are the weakest pair under deuteranopia. Dropping olive entirely, as above, removes the problem; use it only if an 8th real group is unavoidable.
- **Two tiers.** Pale fills for sector backgrounds and rings; the same colours darkened 40% for branch strokes. A 0.3pt line in a pastel reads as grey, so strokes must be darker than fills.
  ```r
  fill_pal   <- setNames(tol_light[seq_along(groups)], groups)
  stroke_pal <- colorspace::darken(fill_pal, 0.40)
  ```
- **`geom_hilight(alpha = 1)`.** The default 0.5 composites the fill against what is behind it and invalidates the contrast guarantees. Use an already-pale colour rather than a transparent saturated one.
- **Every ring gets a hairline stroke** (`colour = "grey70", linewidth = 0.1`). Pale fills on white have only 1.4-2.5:1 edge contrast and otherwise vanish in print.

All 7 working colours clear **8.3:1 contrast against black text**, so tip labels sit safely on any sector.

### 4.3 `03_tree_figure.R`, Figure 1

**Design direction:** restrained and minimal, in the spirit of the LiveMap work. The tree is the only ink that matters; everything else recedes.

- One sans-serif family throughout. Three sizes only: tip labels 2.2mm, legend 8pt, title 10pt.
- No panel border, no grid, no background fill, no shadows or gradients.
- Hairlines, not rules. Guide lines `grey85` at 0.15pt.
- Legend sits **inside the empty centre of the circle**. No box, no title where the meaning is obvious. This uses space the layout already wastes and keeps the figure a single object.
- Scale bar bottom-left, thin, with a plain `substitutions/site` label.
- Negative space is deliberate. Do not fill the centre with anything but the legend.

**Geometry, computed for 92 tips rather than guessed:**

```r
ggtree(tr, layout = "fan", open.angle = 12, linewidth = 0.4)
ggsave(width = 12, height = 12, units = "in", dpi = 600)
```

At a 12 inch canvas, `open.angle = 12` and tip labels on a ring at 60% of the radius, each tip gets **6.04 mm of arc** while a 2.2 mm label needs **3.19 mm** of line height. That is **1.89x clearance**, and the longest label (30.4 mm of radial run) still leaves 30 mm of margin inside the canvas. Overlap is therefore impossible by construction, not by inspection.

Do not shrink the canvas below 12 inches without re-running the check in `05_verify.R`.

**Layers, in order:**

1. `geom_hilight()` on each functional group's MRCA, pale fill, `alpha = 1`, `to.bottom = TRUE`.
2. Tree branches, `aes(colour = group)`, `stroke_pal`, `linewidth = 0.4`.
3. `geom_tiplab()` with `align = TRUE`, `linetype = "dotted"`, `linesize = 0.15`, `size = 2.2`, `colour = "grey10"`, `offset` tuned so text clears the sector edge. ggtree computes per-tip rotation automatically in fan layout; **do not** use the legacy `geom_tiplab2()`.
4. `geom_nodepoint()` for support **>= 70 only**, size by band (70-89, 90-94, 95-100), single muted colour, not a rainbow.
5. `ggtreeExtra::geom_fruit(geom = geom_tile)` ring for functional category.
6. A second ring left as a hook for metal, once that mapping exists.
7. `geom_treescale()`, then legends.

**Export:** PDF via `cairo_pdf` (vector, for the journal) and TIFF at 600 dpi (for submission systems that demand raster). Embed fonts.

### 4.4 `04_function_figure.R`, Figure 2

Horizontal bar chart, proteins per functional group, same palette and same group order as Figure 1 so the two read as one system. Sorted by count. Direct value labels at the bar ends rather than an axis, no gridlines. Add a metal facet once the metal mapping lands.

### 4.5 `05_verify.R`, the build gate

This is what makes "no overlapping text" a guarantee rather than a hope. It must **fail the build**, not warn:

1. **Tip label overlap.** Compute each label's angular position and its rendered extent, assert the minimum gap between consecutive labels is positive. Report the worst pair and its margin.
2. **Tip set integrity.** Tree, alignment and metadata tip sets identical, all 92.
3. **Annotation completeness.** Every tip has a group and a display label; no `NA` reaches the plot.
4. **CVD check.** Simulate deuteranopia, protanopia and tritanopia on the final palette; assert minimum pairwise CIE Lab distance >= 10.
5. **Contrast check.** Assert every fill clears 4.5:1 against the label colour.
6. **Visual check.** Render at final size, then inspect the PNG. Mechanical checks cannot catch everything; look at it.

---

## 5. Build order

1. Download `cog-20.def.tab` into `data/annotation/`, commit it.
2. `01_annotate.R`, then **stop and read the group counts** before going further. The grouping in 4.1 is a starting proposal and must be reconciled with the real distribution.
3. `02_palette.R` and its verification.
4. `03_tree_figure.R`, iterating on offsets at full export size, never on a small preview.
5. `05_verify.R`. Fix anything it rejects.
6. `04_function_figure.R`.
7. Final visual pass at 100% zoom and in greyscale.

## 6. Environment

The figures are built on **ape + ggplot2 only**, with no Bioconductor dependency:

```r
install.packages(c("ape", "ggplot2", "dplyr", "readr", "colorspace"))
```

The polar layout is computed directly rather than delegated to `ggtree`. That is
deliberate: the brief required guaranteed non-overlapping labels, and computing
the tip angles explicitly is what makes the clearance check in `05_verify.R`
possible. It also removes a heavy dependency. `ggtree` remains a valid
alternative if the methods section would rather cite it.

Record `sessionInfo()` alongside the figures.

## 7. Open item

**Metal mapping is still missing.** The second ring and the Figure 2 facet stay as hooks until a CSV of `accession,metal` (Zn/Fe/Mg/Mn/Cu/Co) is added to `data/metadata/`. The 14 duplicated tips are exactly the multi-metal proteins, so that file also resolves how the copies should be labelled. Everything else builds without it.

---

## 8. Build record

Both figures are built and the gate passes.

### Decisions taken during the build

**Branch lengths are not drawn.** The branch lengths in this tree span three
orders of magnitude, from 0 to 1.21. Drawn to scale they collapse the entire
topology into the centre of the circle and the figure reads as a starburst
rather than a tree. The figure is therefore a **cladogram**: radius is
topological depth, so every tip lands on one ring and the annotation ring lines
up with it. This must be stated in the caption. It is also consistent with the
support situation: with 73 of 89 nodes below 50% bootstrap, precise branch
lengths would imply a precision the data does not carry.

**Tip labels are dark grey, not coloured by group.** Colouring the text as well
as the ring is redundant, and pale yellows and cyans fail as text on white. The
ring carries the colour, the text carries the name.

**Labels use ASCII only.** A middle dot for copy number rendered as `<c2><b7>`
under the build locale. Duplicated tips are marked `#2` / `#3` instead, and
`05_verify.R` fails if any non-renderable character reaches a label.

### Measured results

| Check | Result |
| --- | --- |
| Tip sets, tree vs alignment vs metadata | identical, 92 |
| Annotation coverage | 92 of 92, 7 groups |
| Palette, worst CVD separation | 10.9 Lab dE (deuteranopia) |
| Palette, text contrast | 8.35:1 minimum against black |
| Tip label clearance | **2.28x** (7.43 mm arc vs 3.26 mm needed) |
| Longest label vs canvas | clears the edge by 11 mm |
| Support shown | 15 of 90 nodes at >= 70 |

### Group sizes, 78 proteins

Amino acid, nucleotide and coenzyme metabolism 17; general and unknown function
14; energy and carbohydrate metabolism 13; cell envelope, division and defense
10; transcription, replication and repair 9; translation and ribosome 8;
inorganic ion transport 7.

### Outputs

`figures/output/fig1_tree.pdf`, `fig1_tree_600dpi.png`,
`fig2_function.pdf`, `fig2_function.png`.

### Caption text for Figure 1

> Circular clustergram of the 92 candidate metal-associated protein entries
> (78 distinct proteins) of *Yersinia pestis* CO92. Branches and the outer ring
> are coloured by eggNOG COG functional category. Tip labels give the gene name
> and UniProt accession; proteins associated with more than one metal appear
> once per association and are marked `#2` or `#3`. Nodes with bootstrap
> support of at least 70% are marked, scaled by support band; unmarked nodes
> fall below 70%. The tree is drawn as a cladogram, so branch lengths are not
> represented. The maximum-likelihood tree was inferred in MEGA and rendered in
> R.

## 9. Still open

**Metal mapping.** A CSV of `accession,metal` for the six metals would add a
second ring and a metal split on Figure 2, and would let the `#2` / `#3` tips
be labelled by metal instead of by number. Everything else is complete.

**Subcellular localization.** The PSORTb output would add the innermost ring,
matching the reference paper's layout and feeding the vaccine-candidate section.
