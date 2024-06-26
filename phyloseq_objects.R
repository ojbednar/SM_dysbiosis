setwd("/Users/ojbednar/Library/CloudStorage/OneDrive-IndianaUniversity/NDI_16S")

library(phyloseq)
library(microViz)
library(tidyverse)

asv <- import_biom("fixed_ndi.biom") #otu and phylogeny
meta <- read_csv("NDI_samples_enroll_ext_20240619_clean_meta.csv", col_select = 2:5) #cleaned metadata from stool samples
new <- read_csv("NDI_Analytic_20240405.csv") #study metadata

meta_new <- left_join(meta, new, by = "studyid") %>% as.data.frame #merge stool and study metadata
rownames(meta_new) <- meta_new$SampleID 
meta_new <- meta_new %>% select(-SampleID) %>% mutate(Group = case_when(group == 1 ~ "CM",
                                                                        group == 2 ~ "RDS",
                                                                        group == 3 ~ "M/S",
                                                                        group == 4 ~ "SMA",
                                                                        group == 5 ~ "Prostration",
                                                                        group == 6 & pcr_fal0 == 0 ~ "PfNeg",
                                                                        group == 6 & pcr_fal0 == 1 ~ "PfPos")) %>%
  select(studyid, uric_se1, month, Group, group, urinehx, acidosis, ldh_se1, ifabp_se1,
        uremia, mod_kidney, mod_brain, mod_liver, mod_heart,mod_blood, mod_lung, mod_multi, mod_systems, sysbp_low, acidotic, sepsis, crea_se1, krt_indic, ngal_pl1, ngal_pl3,
        ceftri, ceftri1, ceftri2, antibiot2, antibiot1, antibio2, chloram, cipro, xpen, gentam, ceftri,flagyl, cotrimoxazole, il10_se1, il2_se1,
        il17a_se1, il1b_se1, plt, ifabp_se1, sevaki0, akd_overall, ngal_highrisk0, akd, hrp2_pl1, doe, stool_dt0,
        chloram1, cipro1, xpen1, gentam1, flagyl1, amoxyl, lacidosis, scd14_pl1, lbp_pl1, tff3_pl1, outcome, sm_death12m,
        akistage0, lactate0hr, death, inpathemosv12m, inpathemosvno12m, site, sex, age, doe, hbsgt, haki, sickvmal6m,
        inpathemosvno12m, hyperpa0, thrombocytopenia, coma, sev_anemia, jaundice_adm, inpatany6m, inpatmal6m, shock,
        enr_hb, coldperi, sirscriteria, inpatmalno12m, inpatmalno6m, inpatmal12m, inpatmalno12m, inpatany12m, tbili_se1, enr_sma, enr_cm, hypoglycemia, iglucose,
        respdist, aniongap, ph, tooweak, convulhx, diarr, cough, fever, sickvhemosv12m, sickvhemosv6m, leukocytosis,
        tce_mort_com, akd, inpatsepsisno6m, hgb1m, il6_se1, tnfa_se1, il17a_se1, il17e_se1, paracetamol, paracet1, hgbbase, 
        stool_dt0, doe, feedhour, drinhour, outpatmal6m, neu, neuabs, hm_se1, ho1_pl1, blcx_pos, 
        inpatmal1stdaysto6m, DNA, wbc, icc1q_se1, hrp2_pl1, follow12date, pi0, ang1_pl1, ang2_pl1, psel_pl1, 
        icam_pl1, vcam_pl1, crp_pl1, esel_pl1, il1ra_se1, pct_pl1, stool0, dod, trem_pl1, inpatsepsis12m, deathdaysto12m)
sample_meta <- sample_data(meta_new) 
asv <- merge_phyloseq(asv, sample_meta) #create phyloseq object

#clean up tax table
tax_table(asv)[, colnames(tax_table(asv))] <- gsub(tax_table(asv)[, colnames(tax_table(asv))],     pattern = "[a-z]__", replacement = "")
colnames(tax_table(asv))[1:8] <- c("Domain",  "Phylum",  "Class",   "Order",   "Family",  "Genus",   "Species", "Strain")
colnames(tax_table(asv))
asv <- tax_fix(asv, verbose = FALSE)
asv <- tax_rename(asv, rank = "Species") %>% ps_mutate(sick = Group %in% c("CM", "SMA", "RDS", "M/S", "Prostration"))

#create necessary phyloseq objects by filtering based on metadata
asv_0 <- asv %>% ps_filter(month == 0) %>%
  ps_mutate(hyper_uri = case_when(uric_se1 > 7 ~ 1,
                                  uric_se1 <= 7 ~ 0),
            high_ldh = case_when(ldh_se1 > 400 ~ 1,
                                 ldh_se1 <=400 ~ 0),
            ifabp_inj = case_when(ifabp_se1 > 3.54 ~ 1,
                                  ifabp_se1 <= 3.54 ~ 0),
            sCD14_inj = case_when(scd14_pl1 > 8.2 ~ 1, 
                                  scd14_pl1 <= 8.2 ~ 0),
            lbp_inj = case_when(lbp_pl1 > 50 ~ 1,
                                lbp_pl1 <= 50 ~ 0),
            tff3_inj = case_when(tff3_pl1 >= 4.078 ~ 1,
                                 tff3_pl1 < 4.078 ~ 0),
            int_inj = case_when(tff3_pl1 >= 4.087 | ifabp_se1 >= 15.433 ~ 1,
                                tff3_pl1 < 4.087 | ifabp_se1 < 15.433 ~ 0),
            research_mal = case_when(plt <= 150 & hrp2_pl1 >= 1000 ~ "RM",
                                     #plt > 150 & hrp2_pl1 < 1000 ~ "CM",
                                     .default = "other"),
            h_uric_and_aki = case_when(hyper_uri == "yes" & haki ==1 ~ "both",
                                       hyper_uri == "no" & haki == 0 ~ "neither"),
            hus = case_when(plt <= 150 & haki == 1 & ldh_se1 > 450 & hgbbase < 10 ~ 1,
                            .default = 0),
            hyperbili = case_when(tbili_se1 > 3 ~ "elevated",
                                  tbili_se1 <= 3 ~ "not elevated"),
            hypoglu = case_when(iglucose < 3.9 ~ "hypoglu",
                                iglucose >= 3.9 ~ "not hypoglu"),
            elev_anion = case_when(aniongap > 12 ~ "yes",
                                   aniongap <= 12 ~ "no"),
            lowph = case_when(ph < 7.35 ~ "yes",
                              ph >= 7.35 ~ "no"),
            date_to_stool = (as.Date(stool_dt0, format="%d/%m/%Y") - as.Date(doe, format="%d/%m/%Y")),
            stool_day = case_when(date_to_stool == 0 ~ "day 0",
                                  date_to_stool == 1 ~ "day 1",
                                  date_to_stool == 2 ~ "day 2",
                                  date_to_stool == 3 ~ "day 3",
                                  date_to_stool == 4 ~ "day 4",
                                  date_to_stool == 5 ~ "day 5",
                                  date_to_stool == 6 ~ "day 6",
                                  date_to_stool == 7 ~ "day 7",
                                  date_to_stool > 7 ~ "past one week")
  )

asv_sick <- asv_0 %>% ps_filter(sick == TRUE)
asv_cc <- asv_0 %>% ps_filter(sick == FALSE)
asv_12 <- asv %>% ps_filter(month == 12)
asv_12_cc <- asv_12 %>% ps_filter(sick == FALSE)
asv_pfpos <- asv %>% ps_filter(Group == "PfPos")
asv_sm <- asv %>% ps_filter(sick == TRUE)


#output plot for exporting and other graphing
plot_data <- asv_sick %>%
  tax_fix() %>% 
  tax_transform("compositional", rank = "Genus") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% 
  samdat_tbl()

#for family analysis
plot_data_family <- asv_sick %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Family") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat(c("Enterobacteriaceae")) %>% # adds Parabacteroides as sample data!
  samdat_tbl()
plot_data_family <- plot_data_family %>% select(Enterobacteriaceae, haki, enr_sma, enr_cm,  tff3_inj, hyper_uri,acidotic)
write.csv(plot_data_family, "sick_0_enterobacteriaceae_clinical_symptoms_20240619.csv")

#for species plotting
asv_0_t <- asv_0 %>% ps_select(c("studyid", "sick"))
plot_data_species <- asv_0_t %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Species") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% # adds Parabacteroides as sample data!
  samdat_tbl()
write_csv(plot_data_species, "species_16S.csv")



##################
#metaphlan
setwd("/Users/ojbednar/Library/CloudStorage/OneDrive-IndianaUniversity/Schmidt lab/microbiome/metaphlan_NDI")
species <- read.csv("abs_merged_species_only.csv")
mphln.species <- species %>% filter(grepl('s__', clade_name)) %>%
  filter(!grepl('t__', clade_name)) #strain is repeat of reads and needs to be removed, could be analyzed separately

#otu table
otu <- subset(mphln.species, select = -clade_name) %>% as.matrix()
rownames(otu) <- paste0("OTU", 1:nrow(otu))
otu_species = otu_table(otu, taxa_are_rows = TRUE)

#tax table
tax <- mphln.species[,1] %>% as.data.frame()
split_data <- strsplit(tax$., "\\|") %>% as.matrix()
split_matrix <- do.call(rbind, split_data)
col_names <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
expanded_df <- as.data.frame(split_matrix, stringsAsFactors = FALSE)
colnames(expanded_df) <- col_names
rownames(expanded_df) <- rownames(otu)
expanded_df <- expanded_df %>% as.matrix()
tax_species <- tax_table(expanded_df)
#add metameta
metadata <- read.csv("NDI_wgs_meta_abs.csv")
metadata <- metadata %>% select(sample, Sample.Name, studyid)
meta_join <- left_join(metadata, new, by = "studyid") #merge stool and study metadata
sample_df <- meta_join %>% mutate(Group = case_when(group == 1 ~ "CM",
                                                    group == 2 ~ "RDS",
                                                    group == 3 ~ "M/S",
                                                    group == 4 ~ "SMA",
                                                    group == 5 ~ "Prostration",
                                                    group == 6 & pcr_fal0 == 0 ~ "PfNeg",
                                                    group == 6 & pcr_fal0 == 1 ~ "PfPos")) %>%
  select(sample, Sample.Name,studyid, uric_se1, Group, group, urinehx, acidosis, ldh_se1, ifabp_se1,
         uremia, mod_kidney, mod_brain, mod_liver, mod_heart,mod_blood, mod_lung, mod_multi, mod_systems, sysbp_low, acidotic, sepsis, crea_se1, krt_indic, ngal_pl1, ngal_pl3,
         ceftri, ceftri1, ceftri2, antibiot2, antibiot1, antibio2, chloram, cipro, xpen, gentam, ceftri,flagyl, cotrimoxazole, il10_se1, il2_se1,
         il17a_se1, il1b_se1, plt, ifabp_se1, sevaki0, akd_overall, ngal_highrisk0, akd, hrp2_pl1, doe, stool_dt0,
         chloram1, cipro1, xpen1, gentam1, flagyl1, amoxyl, lacidosis, scd14_pl1, lbp_pl1, tff3_pl1, outcome, sm_death12m,
         akistage0, lactate0hr, death, inpathemosv12m, inpathemosvno12m, site, sex, age, doe, hbsgt, haki, sickvmal6m,
         inpathemosvno12m, hyperpa0, thrombocytopenia, coma, sev_anemia, jaundice_adm, inpatany6m, inpatmal6m, shock,
         enr_hb, coldperi, sirscriteria, inpatmalno12m, inpatmalno6m, inpatmal12m, inpatmalno12m, inpatany12m, tbili_se1, enr_sma, enr_cm, hypoglycemia, iglucose,
         respdist, aniongap, ph, tooweak, convulhx, diarr, cough, fever, sickvhemosv12m, sickvhemosv6m, leukocytosis,
         tce_mort_com, akd, inpatsepsisno6m, hgb1m, il6_se1, tnfa_se1, il17a_se1, il17e_se1, paracetamol, paracet1, hgbbase, 
         stool_dt0, doe, feedhour, drinhour, outpatmal6m, neu, neuabs, hm_se1, ho1_pl1, blcx_pos, 
         inpatmal1stdaysto6m, wbc, icc1q_se1, hrp2_pl1, follow12date, pi0, ang1_pl1, ang2_pl1, psel_pl1, 
         icam_pl1, vcam_pl1, crp_pl1, esel_pl1, il1ra_se1, pct_pl1, stool0, dod, trem_pl1, inpatsepsis12m, deathdaysto12m) %>%
  mutate(month = case_when(studyid == 1745 ~ 1,
                           studyid == 1703 ~ 1,
                           .default = 0))
sample_df[26,50] <- "11/04/2016" #missing stool_dt0
row.names(sample_df) <- sample_df$sample
samples_df <- sample_df %>% select (-sample)
sample <- sample_data(sample_df)

#now to remove strains since they still complicate the data with repeat reads,
#in theory I could only keep them and agg on species if i wanted...
phy_spec <- merge_phyloseq(otu_species, tax_species, metadat = sample)
tax_table(phy_spec)[, colnames(tax_table(phy_spec))] <- gsub(tax_table(phy_spec)[, colnames(tax_table(phy_spec))],     pattern = "[a-z]__", replacement = "")
tax_table(phy_spec)[, colnames(tax_table(phy_spec))] <- gsub(tax_table(phy_spec)[, colnames(tax_table(phy_spec))],     pattern = "_", replacement = " ")

bacteria_physeq <- phy_spec %>% tax_select(tax_list = c("Archaea", "Eukaryota"), n_typos = 0, ranks_searched = "Domain", deselect = TRUE) %>%
  tax_rename(rank = "Species")

metaphlan_0 <- bacteria_physeq %>% ps_filter(month == 0) %>%
  ps_mutate(sick = Group %in% c("CM", "SMA", "RDS", "M/S", "Prostration")) %>%
  ps_mutate(date_to_stool = (as.Date(stool_dt0, format="%d/%m/%Y") - as.Date(doe, format="%d/%m/%Y"))) %>%
  ps_mutate(stool_day = case_when(date_to_stool == 0 ~ "day 0",
                                  date_to_stool == 1 ~ "day 1",
                                  date_to_stool == 2 ~ "day 2",
                                  date_to_stool == 3 ~ "day 3",
                                  date_to_stool == 4 ~ "day 4",
                                  date_to_stool == 5 ~ "day 5",
                                  date_to_stool == 6 ~ "day 6",
                                  date_to_stool >= 7 ~ "day 7"))

#############
#kraken
setwd("/Users/ojbednar/Library/CloudStorage/OneDrive-IndianaUniversity/Schmidt lab/microbiome/kraken_bracken_output")
kraken_wgs <- import_biom("breports_combined.biom")
#metadata
metadata <- read.csv("NDI_wgs_meta.csv")
metadata <- metadata %>% select(sample, Sample.Name, studyid)
meta_join <- left_join(metadata, new, by = "studyid") #merge stool and study metadata
sample_df <- meta_join %>% mutate(Group = case_when(group == 1 ~ "CM",
                                                    group == 2 ~ "RDS",
                                                    group == 3 ~ "M/S",
                                                    group == 4 ~ "SMA",
                                                    group == 5 ~ "Prostration",
                                                    group == 6 & pcr_fal0 == 0 ~ "PfNeg",
                                                    group == 6 & pcr_fal0 == 1 ~ "PfPos")) %>%
  select(sample, Sample.Name,studyid, uric_se1, Group, group, urinehx, acidosis, ldh_se1, ifabp_se1,
         uremia, mod_kidney, mod_brain, mod_liver, mod_heart,mod_blood, mod_lung, mod_multi, mod_systems, sysbp_low, acidotic, sepsis, crea_se1, krt_indic, ngal_pl1, ngal_pl3,
         ceftri, ceftri1, ceftri2, antibiot2, antibiot1, antibio2, chloram, cipro, xpen, gentam, ceftri,flagyl, cotrimoxazole, il10_se1, il2_se1,
         il17a_se1, il1b_se1, plt, ifabp_se1, sevaki0, akd_overall, ngal_highrisk0, akd, hrp2_pl1, doe, stool_dt0,
         chloram1, cipro1, xpen1, gentam1, flagyl1, amoxyl, lacidosis, scd14_pl1, lbp_pl1, tff3_pl1, outcome, sm_death12m,
         akistage0, lactate0hr, death, inpathemosv12m, inpathemosvno12m, site, sex, age, doe, hbsgt, haki, sickvmal6m,
         inpathemosvno12m, hyperpa0, thrombocytopenia, coma, sev_anemia, jaundice_adm, inpatany6m, inpatmal6m, shock,
         enr_hb, coldperi, sirscriteria, inpatmalno12m, inpatmalno6m, inpatmal12m, inpatmalno12m, inpatany12m, tbili_se1, enr_sma, enr_cm, hypoglycemia, iglucose,
         respdist, aniongap, ph, tooweak, convulhx, diarr, cough, fever, sickvhemosv12m, sickvhemosv6m, leukocytosis,
         tce_mort_com, akd, inpatsepsisno6m, hgb1m, il6_se1, tnfa_se1, il17a_se1, il17e_se1, paracetamol, paracet1, hgbbase, 
         stool_dt0, doe, feedhour, drinhour, outpatmal6m, neu, neuabs, hm_se1, ho1_pl1, blcx_pos, 
         inpatmal1stdaysto6m, wbc, icc1q_se1, hrp2_pl1, follow12date, pi0, ang1_pl1, ang2_pl1, psel_pl1, 
         icam_pl1, vcam_pl1, crp_pl1, esel_pl1, il1ra_se1, pct_pl1, stool0, dod, trem_pl1, inpatsepsis12m, deathdaysto12m) %>%
  mutate(month = case_when(studyid == 1745 ~ 1,
                           studyid == 1703 ~ 1,
                           .default = 0))
sample_df[26,50] <- "11/04/2016" #missing stool_dt0
row.names(sample_df) <- sample_df$sample
samples_df <- sample_df %>% select (-sample)
sample <- sample_data(sample_df)

krak = merge_phyloseq(kraken_wgs, metadat = sample)
krak
##now clean up the data, first fix up the tax table to have the correct names
tax_table(krak)[, colnames(tax_table(krak))] <- gsub(tax_table(krak)[, colnames(tax_table(krak))],     pattern = "[a-z]__", replacement = "")
colnames(tax_table(krak))[1:7] <- c("Domain",  "Phylum",  "Class",   "Order",   "Family",  "Genus", "Species")

krak <- krak %>% 
  tax_select(tax_list = c("Archaea", "Eukaryota", "Viruses"), n_typos = 0, ranks_searched = "Domain", deselect = TRUE) 
  

tax_table <- tax_table(krak)
tax_table <- as.data.frame(tax_table)
genus <- tax_table$Genus
species <- tax_table$Species

# Combine Genus and Species with a space in between
new_species <- paste(genus, species, sep = " ")

# Assign the modified species names back to the taxonomy table
tax_table$Species <- new_species
tax_table <- as.matrix(tax_table)
# Update the taxonomy table in the phyloseq object
tax_table(krak) <- tax_table

krak <- krak %>% tax_fix() %>% tax_rename(rank = "Species")

krak_0 <- krak %>% ps_filter(month == 0) %>%
  ps_mutate(sick = Group %in% c("CM", "SMA", "RDS", "M/S", "Prostration")) %>%
  ps_mutate(date_to_stool = (as.Date(stool_dt0, format="%d/%m/%Y") - as.Date(doe, format="%d/%m/%Y"))) %>%
  ps_mutate(stool_day = case_when(date_to_stool == 0 ~ "day 0",
                                  date_to_stool == 1 ~ "day 1",
                                  date_to_stool == 2 ~ "day 2",
                                  date_to_stool == 3 ~ "day 3",
                                  date_to_stool == 4 ~ "day 4",
                                  date_to_stool == 5 ~ "day 5",
                                  date_to_stool == 6 ~ "day 6",
                                  date_to_stool >= 7 ~ "day 7+"))
