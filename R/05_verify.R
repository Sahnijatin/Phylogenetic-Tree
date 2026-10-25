# Build gate. Fails loudly rather than warning.
suppressPackageStartupMessages({library(ape); library(readr); library(dplyr)})
source("R/02_palette.R")

fail <- function(...) stop(sprintf(...), call. = FALSE)
ok   <- function(...) cat(sprintf("  PASS  %s\n", sprintf(...)))

cat("1. tip set integrity\n")
tr   <- read.tree("data/tree/mega_tree.nwk")
meta <- read_csv("data/metadata/tips_annotated.csv", show_col_types = FALSE)
fa   <- sub("^>", "", grep("^>", readLines("data/alignment/metalloproteome_aln.fasta"),
                           value = TRUE))
if (!setequal(tr$tip.label, meta$tip_label)) fail("tree and metadata tips differ")
if (!setequal(tr$tip.label, fa))             fail("tree and alignment tips differ")
if (Ntip(tr) != 92L)                         fail("expected 92 tips, got %d", Ntip(tr))
ok("tree, alignment and metadata all carry the same 92 tips")

cat("2. annotation completeness\n")
if (any(is.na(meta$fun_group)))     fail("%d tips lack a functional group",
                                         sum(is.na(meta$fun_group)))
if (any(is.na(meta$display_label))) fail("tips lack a display label")
if (any(!nzchar(meta$gene)))        fail("tips lack a gene name")
bad <- meta$display_label[grepl("<[0-9a-f]{2}>", meta$display_label)]
if (length(bad)) fail("non-renderable characters in labels: %s", bad[1])
ok("all 92 tips carry a gene name, label and group; %d groups in use",
   n_distinct(meta$fun_group))

cat("3. palette\n")
verify_palette(fill_pal)

cat("4. tip label overlap\n")
# Geometry must match R/03_tree_figure.R
OPEN_ANGLE <- 14; R_LAB <- 1.085; LAB_SIZE <- 2.25
CANVAS_IN  <- 12;  LIM <- 1.46
n       <- Ntip(tr)
mm_unit <- (CANVAS_IN * 25.4) / (2 * LIM)      # mm per plot unit
step    <- (360 - OPEN_ANGLE) / n              # degrees between adjacent tips
arc_mm  <- 2 * pi * R_LAB * mm_unit * step / 360
need_mm <- LAB_SIZE * 1.45                     # glyph height plus leading
ratio   <- arc_mm / need_mm
cat(sprintf("    arc per tip %.2f mm, label needs %.2f mm, clearance %.2fx\n",
            arc_mm, need_mm, ratio))
if (ratio < 1) fail("tip labels would overlap (clearance %.2fx)", ratio)
ok("tangential clearance %.2fx at a %d inch canvas", ratio, CANVAS_IN)

cat("5. radial fit\n")
longest_mm <- max(nchar(meta$display_label)) * LAB_SIZE * 0.60
outer      <- R_LAB * mm_unit + longest_mm
half       <- CANVAS_IN * 25.4 / 2
cat(sprintf("    label ring %.0f mm + longest label %.0f mm = %.0f mm vs %.0f mm\n",
            R_LAB * mm_unit, longest_mm, outer, half))
if (outer > half) fail("longest label overflows the canvas")
ok("longest label clears the canvas edge by %.0f mm", half - outer)

cat("\nall checks passed\n")
