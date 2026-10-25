# Figure 2: functional composition of the metalloproteome.
# Not a second tree. The eggNOG annotation describes the same 78 proteins,
# so it is shown as composition rather than as another topology.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(readr)})
source("R/02_palette.R")

meta <- read_csv("data/metadata/tips_annotated.csv", show_col_types = FALSE)

d <- meta |>
  distinct(protein, .keep_all = TRUE) |>        # one row per protein
  count(fun_group, name = "n") |>
  mutate(fun_group = factor(fun_group, levels = FUN_LEVELS)) |>
  arrange(n) |>
  mutate(fun_group = factor(as.character(fun_group),
                            levels = as.character(fun_group)),
         pct = n / sum(n) * 100)

p <- ggplot(d, aes(n, fun_group, fill = as.character(fun_group))) +
  geom_col(width = 0.62, colour = "grey70", linewidth = 0.15) +
  geom_text(aes(label = sprintf("%d  (%.0f%%)", n, pct)),
            hjust = -0.18, size = 2.9, colour = "grey20") +
  scale_fill_manual(values = fill_pal, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = NULL,
       title = "Functional composition of the metalloproteome",
       subtitle = sprintf("%d proteins, eggNOG-mapper COG categories",
                          sum(d$n))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "grey15", size = 9),
        plot.title = element_text(size = 11, colour = "grey10"),
        plot.subtitle = element_text(size = 8.5, colour = "grey45"),
        plot.margin = margin(10, 14, 10, 10))

ggsave("figures/output/fig2_function.png", p, width = 7.2, height = 3.4,
       dpi = 400, bg = "white")
ggsave("figures/output/fig2_function.pdf", p, width = 7.2, height = 3.4,
       device = cairo_pdf, bg = "white")
cat("wrote figures/output/fig2_function.{png,pdf}\n"); print(as.data.frame(d))
