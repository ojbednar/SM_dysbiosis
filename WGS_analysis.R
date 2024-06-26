setwd("/Users/ojbednar/Library/CloudStorage/OneDrive-IndianaUniversity/Schmidt lab/microbiome/16S_till_OneDrive_works")
library(phyloseq)
library(corncob)
library(ggthemes)
# Load your phyloseq objects (replace 'your_phyloseq_object' with your actual object names)
krak <- readRDS("wgs_kraken_phyloseq.rds")
metaphlan <- readRDS("wgs_metaphlan_phyloseq.rds")

#barokits
k_barplot <- krak_0 %>%
  ps_mutate(Outcome = cut(outcome, breaks=2, labels=c("survived", "died"))) %>%
  ps_select(Outcome, Group, sick) %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 7,tax_order = c("Enterobacteriaceae", "Bacteroidaceae",
                                                                    "Enterococcaceae", "Prevotellaceae", 
                                                                    "Lachnospiraceae", 
                                                                    "Bifidobacteriaceae", "Streptococcaceae")) +
  labs(x= NULL, y = NULL, title= "Kraken")
ggsave("ccvsm_family_kraken.png", k_barplot, width = 8.3, height = 10, dpi = 1200, device = "png")

m_barplot <- metaphlan_0 %>%
  ps_mutate(Outcome = cut(outcome, breaks=2, labels=c("survived", "died"))) %>%
  ps_select(Outcome, Group, sick) %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 7, tax_order = c("Enterobacteriaceae", "Bacteroidaceae",
                                                                                "Enterococcaceae", "Prevotellaceae", 
                                                                                "Lachnospiraceae", 
                                                                                "Bifidobacteriaceae", "Streptococcaceae")) +
  labs(x= NULL, y = NULL, title= "Metaphlan sick")
ggsave("ccvsm_family_metaphlan.png", m_barplot, width = 8.3, height = 10, dpi = 1200, device = "png")

##beta diversity
beta_2 <- metaphlan_0 %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "sick", alpha = 0.6, size = 2) +
  theme_classic(12) +
  coord_fixed(0.7) + #scale_color_brewer(palette = "Set1") +
  labs(title = "Metaphlan") + scale_color_calc()
metaphlan_0 %>% tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("Group", "site", "age", "stool_day"), n_perms = 999, seed = 123) %>%
  perm_get()
ggsave("ccvsm_beta_aitchison_metaphlan_sick.png", beta_2, width = 5, height = 7, dpi = 1200, device = "png")
beta_3 <- krak_0 %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "sick", alpha = 0.6, size = 2) +
  theme_classic(12) +
  coord_fixed(0.7) + scale_color_calc() +labs(title = "Kraken")
krak_0 %>% tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("sick", "site", "age", "date_to_stool"), n_perms = 999, seed = 123) %>%
  perm_get()
ggsave("ccvsm_beta_aitchison_kraken_sick.png", beta_3, width = 6, height = 8, dpi = 1200, device = "png")

##alpha diversity table
krak <- readRDS("wgs_kraken_phyloseq.rds")
alpha_stat_krak <-  krak_0 %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_krak_export <- alpha_stat_krak %>% select(shannon_Species, observed_Species, sick, Group) %>% group_by(sick)
write_csv(alpha_krak_export, "kraken_alpha.csv")

alpha_stat_met <-  metaphlan_0 %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_met_export <- alpha_stat_met %>% select(shannon_Species, observed_Species, sick, Group) %>% group_by(sick)
write_csv(alpha_met_export, "met_alpha.csv")

#output tables for species level
met_t <- metaphlan_0 %>% ps_select(c("studyid", "sick"))
plot_data_species_met <- met_t %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Species") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>%
  samdat_tbl()
write_csv(plot_data_species_met, "species_from_met_0.csv")

krak_t <- krak_0 %>% ps_select(c("studyid", "sick"))
plot_data_species_krak <- krak_t %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Species") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% # adds Parabacteroides as sample data!
  samdat_tbl()
write_csv(plot_data_species_krak, "species_from_krak_0.csv")
