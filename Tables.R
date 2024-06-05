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
    ranks = c("Phylum", "Class", "Order", "Family", "Genus"),
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
    title = md("**SM vs CC 16S Differential Abundance Anaylsis**"), subtitle = md("*Genus Level*")
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
  tab_footnote(footnote = "Italicized rows indicate p-value < 0.05.") %>% 
  fmt_number(decimals = 3) %>% 
  fmt_scientific(columns = p.value, decimals = 3, rows = p.value <= 0.001)

gtsave(ccvsm_stats_out, "DA_16S_CCvSM.docx")