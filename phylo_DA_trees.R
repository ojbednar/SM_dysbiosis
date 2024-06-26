
bb_models_asv_0_cc <- asv_cc %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("Group")
  )
bb_stats_asv_0_cc <- taxatree_models2stats(bb_models_asv_0_cc, param = "mu")
bb_stats_asv_0_cc
bb_stats_asv_0_cc %>% taxatree_stats_get()
tree_cc <- bb_stats_asv_0_cc %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  # keep only first 4 plots
  .[1-4] %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) #+ labs(title = "death in children with inpat abx")

tree_cc

key <- taxatree_plotkey(
  data = bb_stats_asv_0_cc,
  taxon_renamer = function(x) stringr::str_remove(x, "[PCFGS]: "),
  # 2 lines of conditions below, for filtering taxa to be labelled
  rank == "Family" | rank == "Genus" |rank == "Species", #& prevalence > 0.2,
  p.value < 0.05, 
  !grepl("Kingdom", taxon)
) +
  # add a bit more space for the longer labels by expanding the x axis
  scale_x_continuous(expand = expansion(mult = 0.2))

ggsave("16S_corncob_asv_0_pfpos_v_pfneg_circle_tree_genus.png", tree_group, width = 13, height = 5.5, dpi = 1200, device = "png")