# Author: Grant Hopkins
# Date: January 21, 2024
#
# Description: this script produces basic Exploratory Data Analysis (EDA )for the
#              studies that contain at least 5 doses, excluding three studies in
#              which the null model performed comparably
#
# Input: all functions included in the GitHub folder
#
# Output: plots in the EDA appendix

####################
# Library packages #
####################

# general analysis package
library(tidyverse) # data wrangling

# MCP-Mod packages
library(clinDR) # gather metaData
library(DoseFinding) # implement MCP-Mod curves

# plot packages
library(pheatmap) # hierarchical clustered heatmap
library(ggh4x) # scale_y_dendrogram on ggplot
library(dendextend) # flip y-axis of dendrogram
library(ggplotify) # convert pheatmap to ggplot

# table packages
library(matrixStats) # produce row ranks

####################
# Attach function #
####################

# general functions
source("pool_sum_stat.R")
source("wrangle_data.R")

################
# Wrangle data #
################

data(metaData)
md_list <- wrangle_data(metaData, min_num_dose = 5)$md_list
md_proc <- wrangle_data(metaData, min_num_dose = 5)$md_proc
rm(list = c("metaData"))

################################
# Manually correct/adjust data #
################################

# manually edit id1012_1 with non-placebo-adjusted raw data found via
# https://www.tga.gov.au/sites/default/files/auspar-macitentan-140428-cer.docx

md_list[["id1012_1"]] <- md_list[["id1012_1"]] %>%
  arrange(dose) %>%
  mutate(dose     = c(  0.0,   0.3,    1.0,    3.0,   10.0)) %>%
  mutate(rslt     = c( -7.9,  -8.8,  -10.0,  -10.8,  -11.8)) %>%
  mutate(sd       = c(  8.0,   8.6,    6.7,    7.3,    6.5)) %>%
  mutate(se       = c(  1.1,   1.2,    0.9,    1.0,    0.9)) %>%
  mutate(sampsize = c( 54  ,  54  ,   60  ,   57  ,   56  )) %>%
  mutate(lcl      = c(-10.1, -11.2,  -11.7,  -12.7,  -13.5)) %>%
  mutate(ucl      = c( -5.7,  -6.4,   -8.3,   -8.8,  -10.0))

md_proc[md_proc$trial_id == "id1012_1", ] <- md_list[["id1012_1"]]

##################################
# Manually remove good null fits #
##################################

# identify data to be removed
data_rm <- c("id1015_1", "id1041_1", "id2_1")

# filter accordingly
md_proc <- md_proc %>%
  filter(!(trial_id %in% data_rm))
md_list <- md_list[which(!(names(md_list) %in% data_rm))]

# validate it
md_proc %>% pull(trial_id) %>% unique %>% length
length(md_list)

############################
# Wrangle the data for EDA #
############################

dataEDA <- md_proc %>%
  # calculate summary statistics within each trial
  group_by(trial_id) %>%
  
    # total sample size
    mutate(totalss = sum(sampsize),
    
    # average uncertainty in parameter mean
    mean_se = mean(se),
    
    # expansiveness of doses
    dose_spread_1 = var(dose),
    dose_spread_2 = diff(range(dose)),
    
    # direction of general effect
    trend = ifelse(rslt[which.max(dose)] - rslt[which.min(dose)] > 0,
                  "UPWARD", "DOWNWARD"),
    
    # average of missing data proportion
    totpmiss = sum(sampsize * pmiss) / totalss,
    
    # which regimens were used
    regimens = paste0(sort(unique(regimen)), collapse = ",")) %>%
  ungroup() %>%
  
  # now identify if the trial was FDA or Pfizer
  mutate(org = case_when(metasource %in% c("FDA914", "FDA1417") ~ "FDA",
                         metasource %in% c("PF09", "PFIZERUPDATE18") ~ "Pfizer")) %>%
  
  # remove dose-level information
  select(-c(dose, rslt, se, sd, sampsize, ittsize, pmiss, lcl, ucl, alpha, regimen)) %>%
  
  # keep only unique rows
  distinct()


##############
# Make plots #
##############

# establish plot colors
my_cols <- c("#000000", "#56B4E9", "#CC79A7", "#F0E442",
             "#009E73", "#D55E00", "#0072B2", "#E69F00")


# make sure save directories exist
dir.create(file.path(getwd(), file.path("plots")), showWarnings = FALSE)
dir.create(file.path(getwd(), file.path("plots/EDA")), showWarnings = FALSE)


# Was an active comparator used?
tiff(filename = "./plots/EDA/EDA1.tiff",
     height = 4, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(data = dataEDA %>%
         group_by(actcomp) %>%
         summarize(count = n()),
       aes(x = reorder(actcomp, -count), y = count)) +
  geom_bar(stat = "identity") +
  geom_text(aes(y = count, label = count), vjust = -1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 40)) +
  labs(title = "Active Comparator",
       x = "Control Group", y = "Number of Studies"))
dev.off()

# show included therapeutic areas
tiff(filename = "./plots/EDA/EDA2.tiff",
     height = 4, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(dataEDA %>%
         group_by(broadta) %>%
         mutate(count = n()),
       aes(x = reorder(broadta, -count))) +
  geom_bar() +
  theme_bw() +
  labs(title = "Theraputic Areas",
       x = "Broad Theraputic Area", y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        legend.position = "none"))
dev.off()

# drugs delivery method
tiff(filename = "./plots/EDA/EDA3.tiff",
     height = 4, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(data = dataEDA %>%
         group_by(route) %>%
         summarize(count = n()) %>%
         mutate(route = fct_recode(route,
                                   "SUBCUTANEOUS\nINJECTION" =
                                     "SUBCUTANEOUS INJECTION")),
       aes(x = reorder(route, -count), y = count)) +
  geom_bar(stat = "identity") +
  geom_text(data = dataEDA %>%
              group_by(route) %>%
              summarize(count = n()) %>%
              mutate(route = fct_recode(route,
                                        "SUBCUTANEOUS\nINJECTION" =
                                          "SUBCUTANEOUS INJECTION")),
            aes(y = count, label = count), vjust = -1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 40)) +
  labs(title = "Drug Delivery",
       x = "Method", y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        legend.position = "none"))
dev.off()

# drugs delivery method
tiff(filename = "./plots/EDA/EDA4.tiff",
     height = 4, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(data = dataEDA %>%
         group_by(oralForm) %>%
         summarize(count = n()) %>%
         filter(!is.na(oralForm)),
       aes(x = reorder(oralForm, -count), y = count)) +
  geom_bar(stat = "identity") +
  geom_text(data = dataEDA %>%
              group_by(oralForm) %>%
              summarize(count = n()) %>%
              filter(!is.na(oralForm)),
            aes(y = count, label = count), vjust = -1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 30)) +
  labs(title = "Oral Drug Delivery",
       x = "Method", y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        legend.position = "none"))
dev.off()

# how many times was the same drug looked at?
tiff(filename = "./plots/EDA/EDA5.tiff",
     height = 4, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(data = dataEDA %>%
         group_by(drugid) %>%
         summarize(ntrial = n()) %>%
         group_by(ntrial) %>%
         summarize(count = n()),
       aes(x = reorder(ntrial, -count), y = count)) +
  geom_bar(stat = "identity") +
  geom_text(data = dataEDA %>%
              group_by(drugid) %>%
              summarize(ntrial = n()) %>%
              group_by(ntrial) %>%
              summarize(count = n()),
            aes(y = count, label = count), vjust = -1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 50)) +
  labs(title = "Repeated Drugs",
       x = "Number of Studies", y = "Number of Drugs") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.position = "none"))
dev.off()


# how were repeated doses clustered together?
write.csv(x = dataEDA %>%
            group_by(drugid) %>%
            mutate(ntrial = n()) %>%
            filter(ntrial > 1) %>%
            select(drugid) %>%
            group_by(drugid) %>%
            summarize(count = n()) %>%
            arrange(-count),
          row.names = FALSE,
          file = "./plots/EDA/repeated_drugs.csv")


# parallel versus crossover
tiff(filename = "./plots/EDA/EDA7.tiff",
     height = 4, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(data = dataEDA %>%
         group_by(design) %>%
         summarize(count = n()),
       aes(x = reorder(design, -count), y = count)) +
  geom_bar(stat = "identity") +
  geom_text(data = dataEDA %>%
              group_by(design) %>%
              summarize(count = n()),
            aes(y = count, label = count), vjust = -1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 50)) +
  labs(title = "Study Design",
       x = "Design", y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.position = "none"))
dev.off()

# understand how endpoints were calculated and type
tiff(filename = "./plots/EDA/EDA8.tiff",
     height = 4, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(data = dataEDA %>%
         group_by(primsource, primtype) %>%
         summarize(count = n(), .groups="keep") %>%
         ungroup() %>%
         mutate(primtype = str_to_title(primtype))) +
  geom_bar(aes(x = reorder(primsource, -count), y = count),
           stat = "identity") +
  facet_wrap(~ factor(primtype, levels = c("Continuous", "Binary"))) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 30)) +
  labs(title = "Endpoint Characteristics",
       x = "Reporting Method", y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.position = "none"))
dev.off()

# endpoint description
write.csv(x = dataEDA %>%
              group_by(endpointLong) %>%
              summarize(count = n()) %>%
              arrange(-count),
          row.names = FALSE,
          file = "./plots/EDA/endpoints.csv")

# how many doses were used in each trial
tiff(filename = "./plots/EDA/EDA10.tiff",
     height = 4, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(data = dataEDA %>%
         group_by(ndose) %>%
         summarize(count = n()),
       aes(x = reorder(ndose, -count), y = count)) +
  geom_bar(stat = "identity") +
  geom_text(data = dataEDA %>%
              group_by(ndose) %>%
              summarize(count = n()),
            aes(y = count, label = count), vjust = -1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 30)) +
  labs(title = "Doses Per Trial",
       x = "Number of Doses", y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.position = "none"))
dev.off()

# what was the direction of the trend?
tiff(filename = "./plots/EDA/EDA11.tiff",
     height = 4, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(data = dataEDA %>%
         group_by(trend) %>%
         summarize(count = n()),
       aes(x = reorder(trend, -count), y = count)) +
  geom_bar(stat = "identity") +
  geom_text(data = dataEDA %>%
              group_by(trend) %>%
              summarize(count = n()),
            aes(y = count, label = count), vjust = -1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 30)) +
  labs(title = "Effect Trend",
       x = "Direction", y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.position = "none"))
dev.off()

# CONTINUOUS PLOTS

# size of dose range
tiff(filename = "./plots/EDA/EDA12.tiff",
     height = 2, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(dataEDA, aes(x = dose_spread_2)) +
  geom_histogram(bins = 10) +
  theme_bw() +
  labs(title = "Range of Doses",
       x = "Range of Dose",
       y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.position = "none"))
dev.off()

# sample size
tiff(filename = "./plots/EDA/EDA13.tiff",
     height = 2, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(dataEDA, aes(x = totalss)) +
  geom_histogram(bins = 10) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 15)) +
  labs(title = "Sample Size",
       x = "Total Sample Size",
       y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.position = "none"))
dev.off()

# standard error in means
tiff(filename = "./plots/EDA/EDA14.tiff",
     height = 2, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(dataEDA, aes(x = mean_se)) +
  geom_histogram(bins = 10) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 23)) +
  labs(title = "Standard Error of Estimates",
       x = "Average SE Across Doses",
       y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.position = "none"))
dev.off()

# missing info
tiff(filename = "./plots/EDA/EDA15.tiff",
     height = 2, width = 3, units = 'in', res = 300, compression = "lzw")
print(ggplot(dataEDA %>% filter(!is.na(totpmiss)), aes(x = totpmiss)) +
  geom_histogram(bins = 10) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  labs(title = "Missing Data (if info available)",
       x = "Percentage of Missing Data",
       y = "Number of Studies") +
  guides(fill = guide_legend(title = "Cluster")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.position = "none"))
dev.off()


