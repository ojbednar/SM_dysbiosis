library(microViz)

#16S SM v CC
beta_1 <- asv_0 %>% 
  ps_mutate(group_order = case_when(group %in% c(1,2,3,4,5) ~ "A",
                                    group == 6 ~ "B")) %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "sick", alpha = 0.6, size = 2) +
  theme_classic(12) +
  coord_fixed(0.7) + scale_color_calc()
asv_0 %>% tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("sick", "site", "age", "antibio2", "stool_day"), n_perms = 999, seed = 123) %>%
  perm_get()
ggsave("ccvsm_beta_aitchison_asv_0_sm_v_cc.png", beta_1, width = 6, height = 8, dpi = 1200, device = "png")
#Metaphlan WGS SM v CC
beta_2 <- metaphlan %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "sick", alpha = 0.6, size = 2) +
  theme_classic(12) +
  coord_fixed(0.7) + #scale_color_brewer(palette = "Set1") +
  labs(title = "Metaphlan")

#Kraken WGS SM v CC
beta_3 <- krak %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "sick", alpha = 0.6, size = 2) +
  theme_classic(12) +
  coord_fixed(0.7) + scale_color_calc()
labs(title = "Kraken")
metaphlan %>% tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("sick", "site", "age", "antibio2", "stool_day"), n_perms = 999, seed = 123) %>%
  perm_get()

#16S PfPos vs PfNeg
beta_cc <- asv_0 %>% ps_filter(X12mo_sample == 0, sick == FALSE) %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(color = "Group", alpha = 0.85, size = 3) +
  theme_classic(12) +
  coord_fixed(0.7) + scale_color_brewer(palette = "Paired")
asv_0 %>% ps_filter(sick == FALSE) %>%
  tax_filter(min_prevalence = 1/ 100, verbose = FALSE) %>%
  tax_agg(rank = "Species") %>%
  dist_calc(dist = "aitchison") %>%
  dist_permanova(variables = c("Group", "site", "sex", "age", "antibio2", "stool_day"), n_perms = 999, seed = 123) %>%
  perm_get()

##alpha
#asv_0
alpha_stat <-  asv_0 %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_export <- alpha_stat %>% select(sick, shannon_Species, observed_Species)
write_csv(alpha_export, "asv_0_alpha.csv")
#pfposvpfneg
alpha_stat_cc <-  asv_cc %>%
  ps_calc_diversity(rank = "Species", index = 'shannon') %>%
  ps_calc_richness(rank = "Species", index = 'observed') %>%
  samdat_tbl()
alpha_export_cc <- alpha_stat_cc %>% select(Group, shannon_Species, observed_Species)
write_csv(alpha_export_cc, "asv_0_cc_alpha.csv")


