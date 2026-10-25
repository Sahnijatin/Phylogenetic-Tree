# Figure 1: circular clustergram of the Yersinia pestis metalloproteome.
suppressPackageStartupMessages({
  library(ape); library(ggplot2); library(dplyr); library(readr)
})
source("R/02_palette.R")

# ---------------------------------------------------------------- parameters
OPEN_ANGLE <- 14      # degrees of gap at the top, for the scale bar
START_DEG  <- 90      # first tip sits just clockwise of 12 o'clock
R_INNER    <- 0.30    # empty centre, holds the legend
R_TIP      <- 1.00    # radius the tip ring sits at
R_RING     <- 1.035   # functional-category ring
RING_W     <- 0.032
R_LAB      <- 1.085   # tip text starts here
LAB_SIZE   <- 2.25    # mm
LAB_COL    <- "grey15"

# ---------------------------------------------------------------- data
tr <- read.tree("data/tree/mega_tree.nwk")
tr <- ladderize(tr)
meta <- read_csv("data/metadata/tips_annotated.csv", show_col_types = FALSE)
stopifnot(setequal(tr$tip.label, meta$tip_label))

n  <- Ntip(tr); N <- n + tr$Nnode
meta <- meta[match(tr$tip.label, meta$tip_label), ]

# ---------------------------------------------------------------- layout
# tip order as drawn: cladewise traversal
cw      <- reorder(tr, "cladewise")
tiporder <- cw$edge[cw$edge[, 2] <= n, 2]
ypos <- numeric(N); ypos[tiporder] <- seq_len(n)

# internal node y = mean of its children, computed leaves-first
po <- reorder(tr, "postorder")
for (i in seq_len(nrow(po$edge))) {
  p <- po$edge[i, 1]
  kids <- tr$edge[tr$edge[, 1] == p, 2]
  ypos[p] <- mean(ypos[kids])
}

# Radius is topological depth, not branch length. The branch lengths in this
# tree span three orders of magnitude, which collapses the whole topology into
# the centre and renders the figure unreadable. Drawing it as a cladogram puts
# every tip on one ring, which is also what lets the annotation rings line up.
# Branch lengths are therefore not represented; the caption must say so.
# node height = longest path, in edges, down to a descendant tip.
# Tips are 0 and therefore all land on one ring; the root is deepest.
height <- numeric(N)
for (i in seq_len(nrow(po$edge))) {
  p_ <- po$edge[i, 1]
  kids <- tr$edge[tr$edge[, 1] == p_, 2]
  height[p_] <- max(height[kids]) + 1
}
rad <- R_TIP - height / max(height) * (R_TIP - R_INNER)

step  <- (360 - OPEN_ANGLE) / n
ang   <- START_DEG + (ypos - 1) * step       # degrees
rang  <- ang * pi / 180
px <- function(r, a) r * cos(a)
py <- function(r, a) r * sin(a)

# ---------------------------------------------------------------- grouping
# a branch inherits a colour only if every tip below it shares one group
Descendants_ <- function(tr, node, n) {
  out <- c(); stack <- node
  while (length(stack)) {
    cur <- stack[1]; stack <- stack[-1]
    kids <- tr$edge[tr$edge[, 1] == cur, 2]
    if (!length(kids)) out <- c(out, cur) else stack <- c(stack, kids)
  }
  out
}
desc <- lapply(seq_len(N), function(i) if (i <= n) i else Descendants_(tr, i, n))
grp  <- vapply(desc, function(d) {
  g <- unique(as.character(meta$fun_group[d]))
  if (length(g) == 1L) g else NA_character_
}, character(1))

# ---------------------------------------------------------------- edges
arcs <- do.call(rbind, lapply(unique(tr$edge[, 1]), function(p) {
  kids <- tr$edge[tr$edge[, 1] == p, 2]
  a <- range(rang[kids])
  s <- seq(a[1], a[2], length.out = max(3, ceiling(diff(a) / (pi / 720))))
  data.frame(id = p, x = px(rad[p], s), y = py(rad[p], s),
             grp = grp[p], stringsAsFactors = FALSE)
}))
spokes <- data.frame(
  x    = px(rad[tr$edge[, 1]], rang[tr$edge[, 2]]),
  y    = py(rad[tr$edge[, 1]], rang[tr$edge[, 2]]),
  xend = px(rad[tr$edge[, 2]], rang[tr$edge[, 2]]),
  yend = py(rad[tr$edge[, 2]], rang[tr$edge[, 2]]),
  grp  = grp[tr$edge[, 2]], stringsAsFactors = FALSE)

# ---------------------------------------------------------------- tips
deg <- ang[seq_len(n)] %% 360
flip <- deg > 90 & deg < 270
tips <- data.frame(
  label = meta$display_label,
  grp   = as.character(meta$fun_group),
  a     = rang[seq_len(n)],
  r     = rad[seq_len(n)],
  angle = ifelse(flip, deg + 180, deg),
  hjust = ifelse(flip, 1, 0), stringsAsFactors = FALSE)
tips$x <- px(R_LAB, tips$a); tips$y <- py(R_LAB, tips$a)

guides <- data.frame(x = px(tips$r, tips$a),    y = py(tips$r, tips$a),
                     xend = px(R_LAB - 0.012, tips$a),
                     yend = py(R_LAB - 0.012, tips$a))

ring <- do.call(rbind, lapply(seq_len(n), function(i) {
  a <- rang[i] + c(-1, 1) * (step * pi / 180) / 2
  s <- seq(a[1], a[2], length.out = 6)
  data.frame(id = i,
             x = c(px(R_RING, s), rev(px(R_RING + RING_W, s))),
             y = c(py(R_RING, s), rev(py(R_RING + RING_W, s))),
             grp = as.character(meta$fun_group[i]), stringsAsFactors = FALSE)
}))

# ---------------------------------------------------------------- support
sup <- suppressWarnings(as.numeric(tr$node.label))
sidx <- which(!is.na(sup) & sup >= 70)
nodes <- data.frame(
  x = px(rad[n + sidx], rang[n + sidx]),
  y = py(rad[n + sidx], rang[n + sidx]),
  band = cut(sup[sidx], c(70, 90, 95, 100), right = FALSE,
             labels = c("70-89", "90-94", "95-100"), include.lowest = TRUE))

# ---------------------------------------------------------------- plot
lim <- 1.46
p <- ggplot() +
  geom_polygon(data = ring, aes(x, y, group = id, fill = grp),
               colour = "white", linewidth = 0.12) +
  geom_segment(data = guides, aes(x, y, xend = xend, yend = yend),
               colour = "grey82", linewidth = 0.15, linetype = "12") +
  geom_path(data = arcs, aes(x, y, group = id, colour = grp),
            linewidth = 0.42, lineend = "round") +
  geom_segment(data = spokes, aes(x, y, xend = xend, yend = yend,
                                  colour = grp),
               linewidth = 0.42, lineend = "round") +
  geom_point(data = nodes, aes(x, y, size = band),
             colour = "grey25", stroke = 0) +
  geom_text(data = tips, aes(x, y, label = label, angle = angle,
                             hjust = hjust),
            size = LAB_SIZE, colour = LAB_COL, show.legend = FALSE) +
  scale_fill_manual(values = fill_pal, breaks = FUN_LEVELS, name = NULL) +
  scale_colour_manual(values = stroke_pal, na.value = "grey72",
                      breaks = FUN_LEVELS, guide = "none") +
  scale_size_manual(values = c(0.7, 1.15, 1.6), name = "Bootstrap support") +
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
  theme_void(base_size = 9) +
  theme(legend.position = c(0.5, 0.5),
        legend.justification = c(0.5, 0.5),
        legend.key.size = unit(3.0, "mm"),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.margin = margin(3, 5, 3, 3),
        legend.spacing.y = unit(1, "mm"),
        legend.text = element_text(size = 7, colour = "grey15"),
        legend.title = element_text(size = 7, colour = "grey35"),
        plot.margin = margin(2, 2, 2, 2))

ggsave("figures/output/fig1_tree.png", p, width = 12, height = 12,
       dpi = 200, bg = "white")                       # quick review copy
if (!nzchar(Sys.getenv("TREE_FAST"))) {
  ggsave("figures/output/fig1_tree.pdf", p, width = 12, height = 12,
         device = cairo_pdf, bg = "white")            # vector, for the journal
  ggsave("figures/output/fig1_tree_600dpi.png", p, width = 12, height = 12,
         dpi = 600, bg = "white")
}
cat("wrote figures/output/fig1_tree.{png,pdf}\n")
cat(sprintf("tips %d | support shown %d of %d | tree height %d levels\n",
            n, length(sidx), tr$Nnode, max(height)))
