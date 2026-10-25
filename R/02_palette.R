# Pastel, colour-blind-safe categorical palette, with verification.
# Paul Tol "light" qualitative scheme, SRON/EPS/TN/09-002, figure 4.
suppressPackageStartupMessages(library(colorspace))

TOL_LIGHT <- c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB",
               "#99DDFF", "#44BB99", "#DDDDDD")

FUN_LEVELS <- c("Translation and ribosome",
                "Transcription, replication and repair",
                "Energy and carbohydrate metabolism",
                "Amino acid, nucleotide and coenzyme metab.",
                "Inorganic ion transport",
                "Cell envelope, division and defense",
                "General and unknown function")   # takes the reserved grey

fill_pal   <- setNames(TOL_LIGHT, FUN_LEVELS)
stroke_pal <- setNames(colorspace::darken(TOL_LIGHT, 0.45), FUN_LEVELS)

# --- verification -----------------------------------------------------------
rel_lum <- function(hex) {
  v <- t(grDevices::col2rgb(hex)) / 255
  lin <- ifelse(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055)^2.4)
  as.vector(lin %*% c(0.2126, 0.7152, 0.0722))
}
contrast_vs_black <- function(hex) (rel_lum(hex) + 0.05) / 0.05

min_dE <- function(cols) {
  lab <- as(hex2RGB(cols), "LAB")@coords
  min(dist(lab))
}

verify_palette <- function(pal, min_sep = 10, min_contrast = 4.5) {
  sims <- list(normal = pal, deutan = deutan(pal),
               protan = protan(pal), tritan = tritan(pal))
  seps <- vapply(sims, min_dE, numeric(1))
  cr   <- contrast_vs_black(pal)
  cat("minimum pairwise Lab dE per vision type:\n"); print(round(seps, 1))
  cat(sprintf("contrast vs black text: min %.2f:1, max %.2f:1\n",
              min(cr), max(cr)))
  if (any(seps < min_sep))
    stop(sprintf("palette fails CVD separation: worst %.1f < %.1f",
                 min(seps), min_sep))
  if (any(cr < min_contrast))
    stop(sprintf("palette fails text contrast: worst %.2f < %.2f",
                 min(cr), min_contrast))
  cat("palette PASSES\n")
  invisible(TRUE)
}

if (sys.nframe() == 0L) verify_palette(fill_pal)
