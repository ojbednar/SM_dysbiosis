library(corncob)
library(tidyverse)
library(microViz)
library(gt)

#demographics table
plot_data <- asv_0 %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Family") %>%
  tax_transform("log", zero_replace = "halfmin", chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat() %>% 
  samdat_tbl()
dem_data <- plot_data %>% select(age, sex, site, sick, antibio2, stool_day) %>%
  mutate(sex = case_when(sex == 1 ~ "Male",
                         sex == 2 ~ "Female"),
         site = case_when(site == 1 ~ "Mulago",
                          site == 2 ~ "Jinja"),
         sick = case_when(sick == TRUE ~ "SM",
                          sick == FALSE ~ "CC")) %>%
  gt(rows_label(age = md("**Age**"),
               site = md("**Site**"),
               antibio2 = md("**Prior Antibiotics**"),
               stool_day = md("Day of stool collection post admission")),
    groupname_col = "category"
  ) |> 
  tab_header(
    title = "x.x: Demographic Characteristics",
    subtitle = "x.x.x: Demographic Characteristics - ITT Analysis Set"
  )



ccvsm <- asv_0 %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
    variables = c("sick")
  )
ccvsm <- taxatree_models2stats(ccvsm, param = "mu")
ccvsm
ccvsm_stats <- ccvsm %>% taxatree_stats_get()

#genus
ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_CCvSM_genus.docx")

#speices
ccvsm_stats_out_spec <- ccvsm_stats %>%
  filter(rank == "Species") %>%
  mutate(taxon = str_replace(taxon, "^S:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC 16S Differential Abundance Anaylsis: Species Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out_spec, "DA_16S_CCvSM_species.docx")

#16S family
ccvsm_stats_out_family <- ccvsm_stats %>%
  filter(rank == "Family") %>%
  mutate(taxon = str_replace(taxon, "^F:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC 16S Differential Abundance Anaylsis: Family Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out_family, "DA_16S_CCvSM_family.docx")

#16S 12-month
ccvsm <- asv_12 %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
    variables = c("sick")
  )
ccvsm <- taxatree_models2stats(ccvsm, param = "mu")
ccvsm
ccvsm_stats <- ccvsm %>% taxatree_stats_get()

#genus
ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC 16S Differential Abundance Analysis at 12 Month: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_CCvSM_month12_genus.docx")


#metaphlan family
ccvsm_met <- metaphlan_0 %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
    variables = c("sick")
  )
ccvsm_met <- taxatree_models2stats(ccvsm_met, param = "mu")
ccvsm_met_stats <- ccvsm_met %>% taxatree_stats_get()
ccvsm_met_stats_out_family <- ccvsm_met_stats %>%
  filter(rank == "Family") %>%
  mutate(taxon = str_replace(taxon, "^F:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC WGS (Metaphlan4) Differential Abundance Anaylsis: Family Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_met_stats_out_family, "DA_wgs_ccvsm_met_family.docx")
#metaphlan species
ccvsm_met_stats_out_spec <- ccvsm_met_stats %>%
  filter(rank == "Species") %>%
  mutate(taxon = str_replace(taxon, "^S:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC WGS (Metaphlan) Differential Abundance Anaylsis: Species Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_met_stats_out_spec, "DA_wgs_ccvsm_met_species.docx")
#metaphlan genus
ccvsm_met_stats_out_genus <- ccvsm_met_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC WGS (Metaphlan) Differential Abundance Anaylsis: Genus Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_met_stats_out_genus, "DA_wgs_ccvsm_met_genus.docx")


#kraken family
ccvsm_krak <- krak_0 %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
    variables = c("sick")
  )
ccvsm_krak <- taxatree_models2stats(ccvsm_krak, param = "mu")
ccvsm_krak_stats <- ccvsm_krak %>% taxatree_stats_get()
ccvsm_krak_stats_out_family <- ccvsm_krak_stats %>%
  filter(rank == "Family") %>%
  mutate(taxon = str_replace(taxon, "^F:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC WGS (Kraken2) Differential Abundance Anaylsis: Family Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_krak_stats_out_family, "DA_wgs_ccvsm_krak_family.docx")
#kraken species
ccvsm_krak_stats_out_spec <- ccvsm_krak_stats %>%
  filter(rank == "Species") %>%
  mutate(taxon = str_replace(taxon, "^S:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC WGS (Kraken2) Differential Abundance Anaylsis: Species Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_krak_stats_out_spec, "DA_wgs_ccvsm_krak_species.docx")

#kraken genus
ccvsm_krak_stats_out_genus <- ccvsm_krak_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM vs CC WGS (Kraken2) Differential Abundance Anaylsis: Genus Level**"))%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05. (p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_krak_stats_out_genus, "DA_wgs_ccvsm_krak_genus.docx")

####CC for figure 2
cc <- asv_0 %>% ps_filter(sick == FALSE) %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("Group")
  )
cc <- taxatree_models2stats(cc, param = "mu")
cc
ccvsm_stats <- cc %>% taxatree_stats_get()

#genus- PfPos v PfNeg
ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**PfPos vs PfNeg 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_PfPos_v_PfNeg_genus.docx")

#
pfpos <- asv %>% ps_filter(Group == "PfPos") %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("month")
  )
pfpos <- taxatree_models2stats(pfpos, param = "mu")
pfpos
ccvsm_stats <- pfpos %>% taxatree_stats_get()

#genus- PfPos v PfNeg
ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**PfPos Month 12 vs PfPos Enrollment 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_PfPos_month12_v_month0_genus.docx")

#SM 0 vs 12m
#
SM <- asv %>% ps_filter(sick == TRUE) %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("month")
  )
SM <- taxatree_models2stats(SM, param = "mu")
SM
ccvsm_stats <- SM %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Month 12 vs SM Enrollment 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_month12_v_month0_genus.docx")

#month 12 SM vs month 12 CC
#SM 0 vs 12m
month12 <- asv %>% ps_filter(month == 12) %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("sick")
  )
month12 <- taxatree_models2stats(month12, param = "mu")
month12
ccvsm_stats <- month12 %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Month 12 vs CC Month 12 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM month 12, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_month12_v_CC_month12_genus.docx")

#sick 
#prior ABX
sick <- asv_sick %>% ps_mutate(prior_abx_oral = case_when(amoxyl == 1 ~ 1,
                                                        cipro == 1 ~ 1,
                                                        cotrimoxazole == 1 ~ 1,
                                                        flagyl == 1 ~ 1,
                                                        othantibi == "Ampiclox" ~ 1,
                                                        othantibi == "Azithromycin syrup" ~ 1,
                                                        othantibi == "Cefladriox" ~ 1,
                                                        othantibi == "Cephalexin" ~ 1,
                                                        othantibi == "Doxycline" ~ 1,
                                                        othantibi == "Erythromycin" ~ 1,
                                                        .default = 0
)) %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("prior_abx_oral")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Prior Oral Antibiotics 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_prior_abx_genus.docx")

#feedhour
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("feedhour")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Hours Since Eating 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_feedhour_genus.docx")

#neuabs
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("neuabs")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Neutrophil Absolute Number 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_neuabs_genus.docx")

#wbc
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("wbc")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM White Blood Cell Count 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_wbc_genus.docx")

#HO-1
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("ho1_pl1")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Heme Oxygenase 1 Levels 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_ho1_pl1_genus.docx")

#uric
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("uric_se1")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Uric Acid Levels 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_uric_se1_genus.docx")

#enr_sma
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("enr_sma")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Severe Malaria Anemia 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_enr_sma_genus.docx")

#enr_cm
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("enr_cm")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Coma 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_enr_cm_genus.docx")

#lacidosis
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("lacidosis")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Lactic Acidosis 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_lacidosis_genus.docx")

#haki
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("haki")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM AKI 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_AKI_genus.docx")

#tff3_pl1
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("tff3_pl1")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Trefoil Factor 3 Levels 16S Differential Abundance Anaylsis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_tff3_pl1_genus.docx")

#death
sick <- asv_sick %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("death")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Family") %>%
  mutate(taxon = str_replace(taxon, "^F:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM All Mortality 16S Differential Abundance Analysis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_SM_death_family.docx")

#month 1 death
sick <- asv_1 %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
    variables = c("death")
  )
sick <- taxatree_models2stats(sick, param = "mu")
sick
ccvsm_stats <- sick %>% taxatree_stats_get()

ccvsm_stats_out <- ccvsm_stats %>%
  filter(rank == "Genus") %>%
  mutate(taxon = str_replace(taxon, "^G:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**SM Month 1: Post-Discharge Death 16S Differential Abundance Analysis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_month_1_SM_death_genus.docx")


###asv_cc combined with asv_1
# first make sure the rank_names are the same format (e.g. both lowercase)
colnames(asv_cc@tax_table) <- rank_names(asv_cc) %>% tolower()
colnames(asv_1@tax_table) <- rank_names(asv_1) %>% tolower()

# create a dataset ID variable to distinguish the two datasets after merging
asv_cc <- asv_cc %>% ps_mutate(dataset = "asv_0_cc")
asv_1 <- asv_1 %>% ps_mutate(dataset = "asv_1")

combined <- phyloseq::merge_phyloseq(
  asv_cc %>% tax_agg("genus") %>% ps_get(),
  asv_1 %>% tax_agg("genus") %>% ps_get()
)
combined

bb_models_asv_0_group <- combined %>%
  tax_fix() %>%
  tax_prepend_ranks() %>%
  tax_filter(min_prevalence = 0.1) %>%
  taxatree_models(
    type = corncob::bbdml,
    ranks = c("phylum", "class", "order", "family", "genus"),
    variables = c("sick")
  )
bb_stats_asv_0_group <- taxatree_models2stats(bb_models_asv_0_group, param = "mu")
bb_stats_asv_0_group
combined_stats<- bb_stats_asv_0_group %>% taxatree_stats_get()

combined_stats_out <- combined_stats %>%
  filter(rank == "genus") %>%
  mutate(taxon = str_replace(taxon, "^g:", "")) %>%
  select(taxon, estimate, std.error, t.statistic, p.value) %>%
  arrange((p.value)) %>%
  mutate(sig = case_when(
    between(p.value, 0, 0.001) ~ "***",
    between(p.value, 0, 0.01) ~ "**",
    between(p.value, 0, 0.05) ~ "*",
    TRUE ~ "."
  ))%>%
  gt() %>% tab_header(
    title = md("**CC enrollment vs Month 1 SM 16S Differential Abundance Analysis: Genus Level**")
  )%>%
  cols_label(taxon = md("**Taxon**"),
             estimate = md("**Estimate**"),
             std.error = md("**Std Error**"),
             t.statistic = md("**t.statistic**"),
             p.value = md("**P-value**")) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkblue"),
                  targets = "row", fn=function(x) x<0) %>% 
  tab_style_body( columns = estimate,
                  style = cell_text(color = "darkred"),
                  targets = "row", fn=function(x) x>0) %>% 
  tab_style_body( columns = p.value,
                  style = cell_text(style="italic"),
                  targets = "row", fn=function(x) x<0.05) %>% 
  tab_footnote(footnote = "Red indicates positive association with SM, blue indicates a negative association.") %>% 
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.(p<0.001=***, p<0.01=**, p<0.05=*)") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(combined_stats_out, "DA_16S_CC_0_SM_month1_genus.docx")
