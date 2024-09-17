setwd("/Users/ojbednar/Library/CloudStorage/OneDrive-IndianaUniversity/Schmidt lab/WGS_pathway/humann")
library(readr)
library(tidyverse)
library(Maaslin2)

#pathways analysis
parsed_input_data <- data.frame(data.table::fread("all_pathabundance_231207_cpm.tsv", header = TRUE, sep = "\t"), row.names = 1)
parsed_genes <- data.frame(data.table::fread("all_genefamilies_03112024_cpm.tsv", header = TRUE, sep = "\t"), row.names = 1)

path_0 <- parsed_input_data %>% select(!c("CM_8_merge_Abundance","SMA_12_merge_Abundance"))
gene_0 <- parsed_genes %>% select(!c("CM_8_merge_Abundance.RPKs","SMA_12_merge_Abundance.RPKs"))

wgs_meta <- read.csv("NDI_wgs_meta_fix.csv")
wgs_meta$Sample.Name <- str_replace(wgs_meta$Sample.Name, "^(\\D+_\\d+).*", "\\1_merge_Abundance")
row.names(wgs_meta) <- wgs_meta$Sample.Name

wgs_meta <- wgs_meta %>% mutate(sick = Group %in% c("CM", "SMA", "RDS", "M/S", "Prostration")) %>%
  mutate(hyper_uri = case_when(uric_se1 > 7 ~ "yes",
                               uric_se1 <= 7 ~ "no"),
         tff3_inj = case_when(tff3_pl1 >= 4.078 ~ "yes",
                              tff3_pl1 < 4.078 ~ "no"),
         month = case_when(studyid == 1745 ~ 1,
                           studyid == 1703 ~ 1,
                           .default = 0)) %>%
  select(Sample.Name, studyid, sick, month)
wgs_0 <- wgs_meta %>% filter(month == 0)

#pathways
wgs_sm <- Maaslin2(
  input_data = path_0,
  input_metadata = wgs_0,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  output = "wgs_0_maaslin2_SM_v_CC_LM",
  fixed_effects = "sick"
)
#gene families
wgs_meta <- read.csv("NDI_wgs_meta_fix.csv")
wgs_meta$Sample.Name <- str_replace(wgs_meta$Sample.Name, "^(\\D+_\\d+).*", "\\1_merge_Abundance.RPKs")
row.names(wgs_meta) <- wgs_meta$Sample.Name

wgs_meta_gene <- wgs_meta %>% mutate(sick = Group %in% c("CM", "SMA", "RDS", "M/S", "Prostration")) %>%
  mutate(hyper_uri = case_when(uric_se1 > 7 ~ "yes",
                               uric_se1 <= 7 ~ "no"),
         tff3_inj = case_when(tff3_pl1 >= 4.078 ~ "yes",
                              tff3_pl1 < 4.078 ~ "no"),
         month = case_when(studyid == 1745 ~ 1,
                           studyid == 1703 ~ 1,
                           .default = 0)) %>%
  select(Sample.Name, studyid, sick, month)
wgs_gene_0 <- wgs_meta_gene %>% filter(month == 0)

wgs_genes <- Maaslin2(
  input_data = gene_0,
  input_metadata = wgs_gene_0,
  normalization = "NONE",
  transform = "NONE",
  analysis_method = "LM",
  output = "wgs_maaslin2_genes_0_CC_vs_SM_LM",
  fixed_effects = "sick"
)

#graphing
gene_sig <- read_csv("wgs_maaslin2_genes_0_CC_vs_SM_LM/gene_sig_sick_0.csv")
gene_sig <- gene_sig %>% mutate(exp = exp(coef)) %>% mutate(fc = log2(exp)) ##converting coef to log2 fold change will help interpret the data
gene_sig_unclass_only <- gene_sig %>% filter(!str_detect(feature, "g__"))
gene_sig_classified <- gene_sig %>% filter(str_detect(feature, "g__")) %>% 
  mutate(bug = str_extract(feature, "s__[A-Za-z_]+")) %>%
  mutate(bug = str_replace_all(bug, "s__", "")) %>%
  mutate(bug = str_replace_all(bug, "_", " ")) 
gene_all <- read_tsv("wgs_maaslin2_genes_0_CC_vs_SM_LM/significant_results.tsv")
gene_all <- gene_all %>% filter(!feature %in% c("UNMAPPED", "UniRef90_A0A174IWH6.g__Prevotella.s__Prevotella_copri", "UniRef90_A0A174IWH6"))%>%
  mutate(diff_reg= case_when(coef>0 ~ "UP",
                             coef<0 ~ "DOWN"))
gene_qsig_unclass <- gene_all %>% filter(qval < 0.05) %>% filter(!str_detect(feature, "g__"))
table(gene_qsig_unclass$diff_reg)

path_sig <- read_csv("wgs_0_maaslin2_SM_v_CC_LM/sig_results_path_SM.csv")
path_sig <- path_sig %>% mutate(exp = exp(coef)) %>% mutate(fc = log2(exp)) %>% 
  mutate(bug = str_extract(feature, "s__[A-Za-z_]+")) %>%
  mutate(bug = str_replace_all(bug, "s__", "")) %>%
  mutate(bug = str_replace_all(bug, "_", " ")) 
path_all <- read_tsv("wgs_0_maaslin2_SM_v_CC_LM/significant_results.tsv")
path_all <- path_all %>% 
  mutate(diff_reg= case_when(coef>0 ~ "UP",
                             coef<0 ~ "DOWN"))
path_qsig <- path_all %>% filter(qval < 0.05)
table(path_qsig$diff_reg)

library(viridis)
library(ggpubr)
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))
gene_v <- ggplot(gene_all, aes(x=coef, y=-log10(qval), col = diff_reg)) + 
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "darkred"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Negative assoc", "Positive assoc")) +
  labs(title = "Significant Gene Families")
ggsave("volacno_genefamilies.png", plot = gene_v, height = 5, width = 8)


path_v <- ggplot(path_all, aes(x=coef, y=-log10(qval), col = diff_reg)) + 
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("#00AFBB", "darkred"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Negative assoc", "Positive assoc")) +
  labs(title="Significant Pathways") #+ theme_pubclean()
ggsave("volacno_pathway.png", plot = path_v, height = 5, width = 8)

#tables for categories
gene_sig_unclass_only %>% group_by(Category) %>% summarize(n=n())
view(gene_sig_classified %>% group_by(organism) %>% summarize(n=n()))
path_sig %>% group_by(category) %>% summarize(n=n())
path_sig %>% group_by(bug) %>% summarize(n=n())

#table for bugs
gene_sig_classified %>% group_by(bug) %>% summarize(n=n())
path_sig %>% group_by(bug) %>% summarize(n=n())



