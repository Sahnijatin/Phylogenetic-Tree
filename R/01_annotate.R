# Join eggNOG-mapper output and COG functional categories onto the tip table.
# Output: data/metadata/tips_annotated.csv
suppressPackageStartupMessages({library(readr); library(dplyr)})

tips <- read_csv("data/metadata/tips.csv", show_col_types = FALSE) |>
  select(tip_label, db, accession, entry_name, copy, reviewed)

# emapper was run on the 92 protein-metal pairs, so multi-metal proteins
# appear more than once. The functional content of those duplicate rows is
# identical; only the alignment e-value/score differs for one query. Keep the
# best-scoring row per protein so the join is one-to-one.
ann <- read_tsv("data/annotation/query.emapper.annotations",
                comment = "##", show_col_types = FALSE) |>
  rename(query = `#query`) |>
  group_by(query) |>
  slice_max(score, n = 1, with_ties = FALSE) |>
  ungroup()
stopifnot(nrow(ann) == 78L)

cogdef <- read_tsv("data/annotation/cog_categories.tsv",
                   col_names = c("cog_id", "cog_letters", "cog_name"),
                   show_col_types = FALSE)

# COG6021 (bepA) is absent from the COG 2014 table. It is a periplasmic
# M48 zinc metalloprotease acting in outer-membrane protein quality control,
# which places it in category O (post-translational modification, protein
# turnover, chaperones).
cogdef <- bind_rows(cogdef, tibble(
  cog_id = "COG6021", cog_letters = "O",
  cog_name = "Periplasmic M48 metalloprotease BepA, OMP quality control"))

tips <- tips |>
  mutate(protein = sub("_copy\\d+$", "", tip_label)) |>
  left_join(ann, by = c("protein" = "query")) |>
  left_join(cogdef, by = c("COG_category" = "cog_id"))

stopifnot(nrow(tips) == 92L, !any(is.na(tips$Preferred_name)))

# Primary category letter. A COG may list several; the first is primary.
tips <- tips |>
  mutate(cog_primary = ifelse(is.na(cog_letters), NA_character_,
                              substr(cog_letters, 1, 1)))

# Collapse ~21 COG letters into 7 groups. Twenty-one colours would be an
# unreadable legend; these are the standard COG super-categories, merged
# until every group holds enough proteins to be worth a colour.
GROUPS <- list(
  "Translation and ribosome"                   = c("J"),
  "Transcription, replication and repair"      = c("K", "L"),
  "Energy and carbohydrate metabolism"         = c("C", "G"),
  "Amino acid, nucleotide and coenzyme metab." = c("E", "F", "H", "I"),
  "Inorganic ion transport"                    = c("P"),
  "Cell envelope, division and defense"        = c("D", "M", "N", "O", "T", "U", "V")
)
lookup <- stack(setNames(GROUPS, names(GROUPS)))
map <- setNames(as.character(lookup$ind), lookup$values)

tips <- tips |>
  mutate(fun_group = unname(map[cog_primary]),
         fun_group = ifelse(is.na(fun_group), "General and unknown function",
                            fun_group),
         gene = Preferred_name,
         display_label = sprintf("%s (%s)", gene, accession),
         # distinguish the duplicated multi-metal tips
         display_label = ifelse(copy > 1,
                                sprintf("%s (%s) #%d", gene, accession, copy),
                                display_label))

LEVELS <- c(names(GROUPS), "General and unknown function")
tips$fun_group <- factor(tips$fun_group, levels = LEVELS)
stopifnot(!any(is.na(tips$fun_group)))

out <- tips |>
  select(tip_label, protein, accession, gene, display_label, reviewed, copy,
         cog_id = COG_category, cog_primary, cog_letters, cog_name,
         fun_group, pfams = PFAMs, kegg_ko = KEGG_ko, ec = EC)
write_csv(out, "data/metadata/tips_annotated.csv")

cat(sprintf("%d tips annotated, %d distinct proteins\n",
            nrow(out), dplyr::n_distinct(out$protein)))
print(as.data.frame(table(out$fun_group)))
cat(sprintf("\nlongest display label: %d chars (%s)\n",
            max(nchar(out$display_label)),
            out$display_label[which.max(nchar(out$display_label))]))
