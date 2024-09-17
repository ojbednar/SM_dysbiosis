library(phyloseq)
library(tidyverse)
library(microViz)
library(corncob)

asv <- read.csv("ndi_asv.rds")
asv_0 <- asv %>% ps_filter(month == 0)
asv_cc <- asv %>% ps_filter(month == 0)%>% 
  ps_filter(sick == FALSE)
asv_pfpos <- asv %>% ps_filter(Group == "PfPos")
asv_sm <- asv %>% ps_filter(sick == TRUE)

#SM vs CC
bb_models_asv_0_group <- asv_0 %>% 
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("sick")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
bb_stats_asv_0_group %>% taxatree_stats_get()
tree_group <- bb_stats_asv_0_group %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 
tree_group 

#PfPos vs PfNeg
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
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 

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

#pfpos month 12 vs month 0
bb_models_asv_pfpos <- asv_pfpos %>% 
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("month")
  )
bb_stats_asv_pfpos<- taxatree_models2stats(bb_models_asv_pfpos, param = "mu")
bb_stats_asv_pfpos
bb_stats_asv_pfpos %>% taxatree_stats_get()
tree_group <- bb_models_asv_pfpos %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 
tree_group 

#SM month 12 vs month 0
bb_models_asv_sm <- asv_sm %>% 
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("month")
  )
bb_stats_asv_sm <- taxatree_models2stats(bb_models_asv_sm, param = "mu")
bb_stats_asv_sm
bb_stats_asv_sm %>% taxatree_stats_get()
tree_group <- bb_models_asv_sm %>%
  taxatree_plots(
    node_size_range = c(1, 4), palette = "Blue-Red 3"
  ) %>%
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  ) 
tree_group 
