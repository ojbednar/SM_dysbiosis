library(corncob)
library(tidyverse)
library(microViz)
library(gt)

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

gtsave(ccvsm_stats_out, "DA_16S_CCvSM.docx")

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
