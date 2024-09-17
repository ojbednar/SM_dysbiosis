
library(tidyverse)
library(ggthemes)
library(ggpubr)
library(gt)
library(gtsummary)

#re-admission analysis
inpat <- new %>%
  drop_na(inpatmal12m) %>%
  count(sick, inpatmal12m) %>%       
  group_by(sick) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatmal12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(sick = factor(sick, levels = c("SM", "CC"))) %>%  # Reorder factor levels
  ggplot(aes(sick, pct, fill = sick)) +
  geom_bar(stat = "identity") +
  ylab("% of children with inpatient malaria over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("red", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) + 
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatmal12m_sick_new.png", inpat, width = 3, height = 6, dpi = 1200, device = "png")

inpat<-new %>%
  drop_na(inpatmal12m) %>%
  count(sick, inpatmal12m) %>%       
  group_by(sick) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(sick, pct, fill=sick) +
  geom_bar(stat="identity") +
  ylab("% of children with inpatient malaria over 12 months") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  theme_pubclean() + scale_color_calc()
ggsave("inpatmal12m_sick.png", inpat, width = 4, height = 6, dpi = 1200, device = "png")

inpat_sep<-new %>%
  drop_na(inpatsepsis12m) %>%
  count(sick, inpatsepsis12m) %>%       
  group_by(sick) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(sick, pct, fill=factor(inpatsepsis12m)) +
  geom_bar(stat="identity") +
  ylab("% of children with inpatient sepsis over 12 months") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  theme_pubclean() + scale_fill_brewer(palette = "Reds")
ggsave("inpatsepsis12m_sick.png", inpat_sep, width = 4, height = 6, dpi = 1200, device = "png")

inpat_sep <- new %>%
  drop_na(inpatsepsis12m) %>%
  count(sick, inpatsepsis12m) %>%       
  group_by(sick) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatsepsis12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(sick = factor(sick, levels = c("SM", "CC"))) %>%  # Reorder factor levels
  ggplot(aes(sick, pct, fill = sick)) +
  geom_bar(stat = "identity") +
  ylab("% of children with inpatient sepsis over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("red", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatsepsis12m_sick_new.png", inpat_sep, width = 3, height = 6, dpi = 1200, device = "png")

outpat<-new %>%
  drop_na(outpatmal12m)%>%
  count(sick, outpatmal12m) %>%       
  group_by(sick) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(sick, pct, fill=factor(outpatmal12m)) +
  geom_bar(stat="identity") +
  ylab("% of children with malaria outpatient visits over 12 months") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  theme_pubclean() + scale_fill_brewer(palette = "Reds")
ggsave("outpatmal12m_sick.png", outpat, width = 4, height = 6, dpi = 1200, device = "png")

outpat <- new %>%
  drop_na(outpatmal12m) %>%
  count(sick, outpatmal12m) %>%       
  group_by(sick) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(outpatmal12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(sick = factor(sick, levels = c("SM", "CC"))) %>%  # Reorder factor levels
  ggplot(aes(sick, pct, fill = sick)) +
  geom_bar(stat = "identity") +
  ylab("% of children with outpatient malaria over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("red", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("outpatmal12m_sick_new.png", outpat, width = 3, height = 6, dpi = 1200, device = "png")

inpat_any <- new %>%
  drop_na(inpatany12m) %>%
  count(sick, inpatany12m) %>%       
  group_by(sick) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(sick, pct, fill=factor(inpatany12m)) +
  geom_bar(stat="identity") +
  ylab("% of children with inpatient for any reason over 12 months") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  theme_pubclean() + scale_fill_brewer(palette = "Reds")
ggsave("inpatany12m_sick.png", inpat_any, width = 4, height = 6, dpi = 1200, device = "png")

inpat_any <- new %>%
  drop_na(inpatany12m) %>%
  count(sick, inpatany12m) %>%       
  group_by(sick) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatany12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(sick = factor(sick, levels = c("SM", "CC"))) %>%  # Reorder factor levels
  ggplot(aes(sick, pct, fill = sick)) +
  geom_bar(stat = "identity") +
  ylab("% of children with any admission over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("red", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatanyl12m_sick_new.png", inpat_any,width = 3, height = 6, dpi = 1200, device = "png")


inpat<-new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(inpatmal12m) %>%
  count(Pf, inpatmal12m) %>%       
  group_by(Pf) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(Pf, pct, fill=factor(inpatmal12m)) +
  geom_bar(stat="identity") +
  ylab("% of children with inpatient malaria over 12 months") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  theme_pubclean() + scale_fill_brewer(palette = "Purples")
ggsave("inpatmal12m_cc.png", inpat, width = 4, height = 6, dpi = 1200, device = "png")

inpat <- new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(inpatmal12m) %>%
  count(Pf, inpatmal12m) %>%       
  group_by(Pf) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatmal12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(Pf = factor(Pf, levels = c("PfPos", "PfNeg"))) %>%  # Reorder factor levels
  ggplot(aes(Pf, pct, fill = Pf)) +
  geom_bar(stat = "identity") +
  ylab("% of children with inpatient malaria over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("lightskyblue2", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatmal12m_cc_new.png", inpat, width = 3, height = 6, dpi = 1200, device = "png")


outpat<- new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>%
  drop_na(outpatmal12m) %>%
  count(Pf, outpatmal12m) %>%       
  group_by(Pf) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(Pf, pct, fill=factor(outpatmal12m)) +
  geom_bar(stat="identity") +
  ylab("% of children with outpatient malaria over 12 months-CC") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  theme_pubclean() + scale_fill_brewer(palette = "Purples")
ggsave("outpatmal12m_cc.png", outpat, width = 4, height = 6, dpi = 1200, device = "png")

outpat <- new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(outpatmal12m) %>%
  count(Pf, outpatmal12m) %>%       
  group_by(Pf) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(outpatmal12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(Pf = factor(Pf, levels = c("PfPos", "PfNeg"))) %>%  # Reorder factor levels
  ggplot(aes(Pf, pct, fill = Pf)) +
  geom_bar(stat = "identity") +
  ylab("% of children with outpatient malaria over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("lightskyblue2", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("outpatmal12m_cc_new.png", outpat, width = 3, height = 6, dpi = 1200, device = "png")


inpat_sep<-new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(inpatsepsis12m) %>%
  count(Pf, inpatsepsis12m) %>%       
  group_by(Pf) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(Pf, pct, fill=factor(inpatsepsis12m)) +
  geom_bar(stat="identity") +
  ylab("% of children with sepsis over 12 months-CC") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  theme_pubclean() + scale_fill_brewer(palette = "Purples")
ggsave("inpatsepsis12m_cc.png", inpat_sep, width = 4, height = 6, dpi = 1200, device = "png")

inpat_sep <- new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(inpatsepsis12m) %>%
  count(Pf, inpatsepsis12m) %>%       
  group_by(Pf) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatsepsis12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(Pf = factor(Pf, levels = c("PfPos", "PfNeg"))) %>%  # Reorder factor levels
  ggplot(aes(Pf, pct, fill = Pf)) +
  geom_bar(stat = "identity") +
  ylab("% of children with outpatient malaria over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 2.5),
            size = 7,
            fontface = "bold", 
            color = "black") +
  scale_fill_manual(values = c("lightskyblue2", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatsepsis12m_cc_new.png", inpat_sep, width = 3, height = 6, dpi = 1200, device = "png")

inpat_any_cc<-new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(inpatany12m) %>%
  count(Pf, inpatany12m) %>%       
  group_by(Pf) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(Pf, pct, fill=factor(inpatany12m)) +
  geom_bar(stat="identity") +
  ylab("% of children with hospital admission over 12 months-CC") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  theme_pubclean() + scale_fill_brewer(palette = "Purples")
ggsave("inpatany12m_cc.png", inpat_any_cc, width = 4, height = 6, dpi = 1200, device = "png")

inpat_any <- new %>%
  filter(group == 6) %>%
  mutate(Pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                        pcr_fal0 == 0 ~ "PfNeg")) %>% 
  drop_na(inpatany12m) %>%
  count(Pf, inpatany12m) %>%       
  group_by(Pf) %>%
  mutate(total_pct = sum(n)) %>%   # Calculate total within each group
  mutate(pct = (n / sum(n)) * 100) %>%  # Calculate percentage before filtering
  filter(inpatany12m == 1) %>%  # Keep only rows where inpatmal12m is 1
  ungroup() %>%
  mutate(Pf = factor(Pf, levels = c("PfPos", "PfNeg"))) %>%  # Reorder factor levels
  ggplot(aes(Pf, pct, fill = Pf)) +
  geom_bar(stat = "identity") +
  ylab("% of children with any admission over 12 months") +
  geom_text(aes(label = n),  # Display count instead of percentage
            position = position_stack(vjust = 0.5),
            size = 7,
            fontface = "bold", 
            color = "white") +
  scale_fill_manual(values = c("lightskyblue2", "navyblue"))+
  scale_y_continuous(limits = c(0, 100)) +
  theme(legend.position = "none") +
  theme_pubclean()
ggsave("inpatanyl12m_cc_new.png", inpat_any, width = 3, height = 6, dpi = 1200, device = "png")

ane <- new %>% filter( hgbbase < 10.5)
##for fischer's exact test, one tailed\
table(new$sick, new$inpatany12m)
table(new$sick, new$inpatsepsis12m)
table(new$sick, new$outpatmal12m)
table(new$sick, new$inpatmal12m)

table(ane$sick, ane$inpatany12m)
table(ane$sick, ane$inpatsepsis12m)
table(ane$sick, ane$outpatmal12m)
table(ane$sick, ane$inpatmal12m)

inpat_sep<-ane %>%
  drop_na(inpatsepsis12m) %>%
  count(sick, inpatsepsis12m) %>%       
  group_by(sick) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(sick, pct, fill=factor(inpatsepsis12m)) +
  geom_bar(stat="identity") +
  ylab("% of children with sepsis over 12 months-CC") +
  geom_text(aes(label=paste0(sprintf("%1.1f", pct),"%")),
            position=position_stack(vjust=0.5)) +
  theme_pubclean() + scale_fill_brewer(palette = "Reds")
ggsave("inpatsepsis12m_sick_anemic.png", inpat_sep, width = 4, height = 6, dpi = 1200, device = "png")
cc <- new %>% filter(group == 6) %>% mutate(pf = case_when(pcr_fal0 == 1 ~ "PfPos",
                                                           pcr_fal0 == 0 ~ "PfNeg"))

table(cc$pf, cc$inpatany12m)
table(cc$pf, cc$inpatsepsis12m)
table(cc$pf, cc$outpatmal12m)
table(cc$pf, cc$inpatmal12m)

#for adjusted odds ratios
new %>%
  select(sick, outpatmal12m,inpatmal12m, inpatsepsis12m, inpatany12m) %>%
  tbl_uvregression(
    method = glm,
    y = sick,
    method.args = list(family = binomial),
    exponentiate = TRUE,
    pvalue_fun = ~ style_pvalue(.x, digits = 2)
  ) %>%
  add_global_p() %>% # add global p-value
  add_nevent() %>% # add number of events of the outcome
  add_q() %>% # adjusts global p-values for multiple testing
  bold_p() %>% # bold p-values under a given threshold (default 0.05)
  bold_p(t = 0.10, q = TRUE) %>% # now bold q-values under the threshold of 0.10
  bold_labels() %>%
  as_gt() %>%
  gt::gtsave(filename = "OR_cc_v_sm_hospitalization.pdf")
cc %>%
  select(Pf, outpatmal12m,inpatmal12m, inpatsepsis12m, inpatany12m) %>%
  tbl_uvregression(
    method = glm,
    y = Pf,
    method.args = list(family = binomial),
    exponentiate = TRUE,
    pvalue_fun = ~ style_pvalue(.x, digits = 2)
  ) %>%
  add_global_p() %>% # add global p-value
  add_nevent() %>% # add number of events of the outcome
  add_q() %>% # adjusts global p-values for multiple testing
  bold_p() %>% # bold p-values under a given threshold (default 0.05)
  bold_p(t = 0.10, q = TRUE) %>% # now bold q-values under the threshold of 0.10
  bold_labels() %>%
  as_gt() %>%
  gt::gtsave(filename = "OR_pfpos_v_pfneg_hospitalization.pdf")