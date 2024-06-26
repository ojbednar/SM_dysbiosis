

##barplots###
#family
ccvsm_barplot <- asv_0 %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 8, tax_order = c("Enterobacteriaceae","Enterococcaceae",  
                                                                    "Bacteroidaceae", "Tannerellaceae","Lachnospiraceae",
                                                                    "Prevotellaceae","Bifidobacteriaceae","Streptococcaceae")
               ) +
  labs(x= NULL, y = NULL, title = "sick")
ggsave("ccvsm_family_16S.png", ccvsm_barplot, width = 8.3, height = 10, dpi = 1200, device = "png") 

k_barplot <- krak_0 %>%
  ps_mutate(Outcome = cut(outcome, breaks=2, labels=c("survived", "died"))) %>%
  ps_select(Outcome, Group, sick) %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 8, tax_order = c("Enterobacteriaceae","Enterococcaceae",  
                                                                                "Bacteroidaceae", "Tannerellaceae","Lachnospiraceae",
                                                                                "Prevotellaceae","Bifidobacteriaceae","Streptococcaceae")) +
  labs(x= NULL, y = NULL, title= "Kraken")
ggsave("ccvsm_family_kraken.png", k_barplot, width = 8.3, height = 10, dpi = 1200, device = "png")

m_barplot <- metaphlan_0 %>%
  ps_mutate(Outcome = cut(outcome, breaks=2, labels=c("survived", "died"))) %>%
  ps_select(Outcome, Group, sick) %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Family", bar_width = 0.8, n_taxa = 8, tax_order = c("Enterobacteriaceae","Enterococcaceae",  
                                                                                "Bacteroidaceae", "Tannerellaceae","Lachnospiraceae",
                                                                                "Prevotellaceae","Bifidobacteriaceae","Streptococcaceae")) +
  labs(x= NULL, y = NULL, title= "Metaphlan sick")
ggsave("ccvsm_family_metaphlan.png", m_barplot, width = 8.3, height = 10, dpi = 1200, device = "png") 

####
#timeline
ccvsm_barplot_12 <- asv_12 %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Genus", bar_width = 0.8, n_taxa = 12, tax_order = c("Prevotella", "Faecalibacterium", "Bifidobacterium", 
                                                                                "Blautia", "Ruminococcus", "Roseburia","Escherichia", 
                                                                                "Enterococcus", "Bacteroides",
                                                                                "Parabacteroides", "Klebsiella", "Shigella")) +
  labs(x= NULL, y = NULL, title = "Month 12")
ccvsm_barplot <- asv_0 %>%
  phyloseq::merge_samples(group = "sick") %>%
  comp_barplot(tax_level = "Genus", bar_width = 0.8, n_taxa = 12, tax_order = c("Prevotella", "Faecalibacterium", "Bifidobacterium", 
                                                                                "Blautia", "Ruminococcus", "Roseburia","Escherichia", 
                                                                                "Enterococcus", "Bacteroides",
                                                                                "Parabacteroides", "Klebsiella", "Shigella")) +
  labs(x= NULL, y = NULL, title = "Enrollement")
library(gridExtra)
both <- grid.arrange(ccvsm_barplot, ccvsm_barplot_12, ncol = 2)
ggsave("ccvsm_genus_enroll_v_12m.png", both, width = 8, height =6 , dpi = 1200, device = "png")

##deathplot
death_barplot <- asv_sick %>%
  ps_select(sick, inpathemosvno12m, hyperpa0, Group, thrombocytopenia, enr_cm, enr_sma, haki,
            tooweak, respdist, acidosis, jaundice_adm, convulhx, lowph, hyperbili, hypoglu, sick, 
            stool_dt0, doe, antibio2, antibiot1, antibiot2, sickvmal6m, inpatmal6m, feedhour, tff3_inj, hyper_uri, death) %>%
  phyloseq::merge_samples(group = "death") %>%
  comp_barplot(tax_level = "Genus", bar_width = 0.8, tax_order = c("Escherichia", "Shigella", "Enterobacter", "Bacteroides")) +
  labs(x= NULL, y = NULL, title = "death")
ggsave("genus_sick_death.png", death_barplot, width = 8.3, height = 10, dpi = 1200, device = "png")


