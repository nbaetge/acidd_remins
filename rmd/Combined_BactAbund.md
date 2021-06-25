Combined\_BactAbund
================
Nicholas Baetge
6/21/2021

# Intro

This document shows how **individual bottle** bacterial abundance data
from ACIDD experiments were processed, QC’d, and analyzed.

``` r
library(tidyverse)
library(readxl)
library(lubridate)
library(ggpubr)
```

# Import data

``` r
excel_sheets("~/GITHUB/acidd_remins/input/ACIDD_Exp_BactAbund.xlsx")
```

    ## [1] "Metadata" "Data"

``` r
acidd_metadata <- read_excel("~/GITHUB/acidd_remins/input/ACIDD_Exp_BactAbund.xlsx", sheet = "Metadata")
glimpse(acidd_metadata)
```

    ## Rows: 84
    ## Columns: 18
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH1…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San D…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, …
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", …
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3,…
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, …
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, …
    ## $ Datetime                <chr> "2017-12-16T21:30", "2017-12-17T10:00", "2017…
    ## $ TOC_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALS…
    ## $ DOC_Sample              <lgl> TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE,…
    ## $ Parallel_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALS…
    ## $ Cell_Sample             <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FAL…
    ## $ DNA_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE…
    ## $ Nutrient_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE…
    ## $ DNA_SampleID            <chr> "ASH171-A0_S293", NA, NA, NA, NA, NA, "ASH171…

``` r
# unique(metadata$Experiment)
# unique(metadata$Location)
# unique(metadata$Bottle)
# unique(metadata$Treatment)

acidd_data <- read_excel("~/GITHUB/acidd_remins/input/ACIDD_Exp_BactAbund.xlsx", sheet = "Data")
glimpse(acidd_data)
```

    ## Rows: 52
    ## Columns: 5
    ## $ Experiment  <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH171", "ASH171…
    ## $ Bottle      <chr> "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B…
    ## $ Timepoint   <dbl> 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, …
    ## $ Cells_ml    <dbl> 130000, 134000, 128000, 155000, 155000, 200000, 377000, 1…
    ## $ Cells_ml_sd <dbl> 20900, 27600, 22200, 25200, 31900, 49100, 59700, 18400, 3…

``` r
acidd_joined <- left_join(acidd_metadata, acidd_data)
```

    ## Joining, by = c("Experiment", "Bottle", "Timepoint")

``` r
# names(joined)
# summary(joined)
glimpse(acidd_joined)
```

    ## Rows: 84
    ## Columns: 20
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH1…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San D…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, …
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", …
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3,…
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, …
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, …
    ## $ Datetime                <chr> "2017-12-16T21:30", "2017-12-17T10:00", "2017…
    ## $ TOC_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALS…
    ## $ DOC_Sample              <lgl> TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE,…
    ## $ Parallel_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALS…
    ## $ Cell_Sample             <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FAL…
    ## $ DNA_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE…
    ## $ Nutrient_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE…
    ## $ DNA_SampleID            <chr> "ASH171-A0_S293", NA, NA, NA, NA, NA, "ASH171…
    ## $ Cells_ml                <dbl> 130000, 134000, 128000, 155000, 155000, 20000…
    ## $ Cells_ml_sd             <dbl> 20900, 27600, 22200, 25200, 31900, 49100, 597…

``` r
####

class_metadata <- read_excel("~/GITHUB/acidd_remins/input/144L_2018_BactAbund.xlsx", sheet = "Metadata")

class_data <- read_excel("~/GITHUB/acidd_remins/input/144L_2018_BactAbund.xlsx", sheet = "Data")

class_joined <- left_join(class_metadata, class_data)
```

    ## Joining, by = c("Bottle", "Timepoint")

``` r
###

mud_data <- read_excel("~/GITHUB/acidd_remins/input/Mud191_BactAbund.xlsx", sheet = "Sheet1") %>% 
  group_by(Bottle, Timepoint) %>% 
  mutate(mean = mean(Cells_ml, na.rm = T),
         Cells_ml_sd = sd(Cells_ml, na.rm = T)) %>% 
  select(-Cells_ml) %>% 
  distinct() %>% 
  rename(Cells_ml = mean)


combined_data <- full_join(acidd_joined, class_joined) %>% 
  full_join(., mud_data)
```

    ## Joining, by = c("Experiment", "Location", "Temperature_C", "Depth", "Bottle", "Timepoint", "Treatment", "Target_DOC_Amendment_uM", "Inoculum_L", "Media_L", "Datetime", "TOC_Sample", "Parallel_Sample", "Cell_Sample", "DNA_Sample", "DNA_SampleID", "Cells_ml")

    ## Joining, by = c("Experiment", "Location", "Temperature_C", "Depth", "Bottle", "Timepoint", "Treatment", "Target_DOC_Amendment_uM", "Inoculum_L", "Media_L", "Datetime", "Cells_ml", "Cells_ml_sd")

``` r
glimpse(combined_data)
```

    ## Rows: 298
    ## Columns: 20
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH1…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San D…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, …
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", …
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3,…
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, …
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, …
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, …
    ## $ Datetime                <chr> "2017-12-16T21:30", "2017-12-17T10:00", "2017…
    ## $ TOC_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALS…
    ## $ DOC_Sample              <lgl> TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE,…
    ## $ Parallel_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALS…
    ## $ Cell_Sample             <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FAL…
    ## $ DNA_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE…
    ## $ Nutrient_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE…
    ## $ DNA_SampleID            <chr> "ASH171-A0_S293", NA, NA, NA, NA, NA, "ASH171…
    ## $ Cells_ml                <dbl> 130000, 134000, 128000, 155000, 155000, 20000…
    ## $ Cells_ml_sd             <dbl> 20900, 27600, 22200, 25200, 31900, 49100, 597…

# Prepare data

Convert date and time column values from characters to dates, add
columns with time elapsed for each experiment, and convert cells/ml to
cells/l, subset data to select only VOI & drop na’s

``` r
cells <- combined_data %>% 
  mutate(Datetime = ymd_hm(Datetime),
         cells = Cells_ml * 1000,
         sd_cells = Cells_ml_sd * 1000) %>%
  group_by(Experiment, Treatment, Bottle) %>% 
  mutate(interv = interval(first(Datetime), Datetime),
         hours = as.numeric(interv/3600),
         days = hours/24) %>% 
  ungroup() %>% 
  select(Experiment:DNA_SampleID, hours, days, cells, sd_cells) %>% 
  drop_na(cells)

glimpse(cells)
```

    ## Rows: 258
    ## Columns: 22
    ## $ Experiment              <chr> "ASH171", "ASH171", "ASH171", "ASH171", "ASH1…
    ## $ Location                <chr> "San Diego", "San Diego", "San Diego", "San D…
    ## $ Temperature_C           <dbl> 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 1…
    ## $ Depth                   <dbl> 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, …
    ## $ Bottle                  <chr> "A", "A", "A", "A", "A", "A", "A", "B", "B", …
    ## $ Timepoint               <dbl> 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 0, …
    ## $ Treatment               <chr> "Control", "Control", "Control", "Control", "…
    ## $ Target_DOC_Amendment_uM <dbl> 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 10,…
    ## $ Inoculum_L              <dbl> 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, …
    ## $ Media_L                 <dbl> 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, …
    ## $ Datetime                <dttm> 2017-12-16 21:30:00, 2017-12-17 10:00:00, 20…
    ## $ TOC_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALS…
    ## $ DOC_Sample              <lgl> TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE,…
    ## $ Parallel_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALS…
    ## $ Cell_Sample             <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRU…
    ## $ DNA_Sample              <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE…
    ## $ Nutrient_Sample         <lgl> TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE…
    ## $ DNA_SampleID            <chr> "ASH171-A0_S293", NA, NA, NA, NA, NA, "ASH171…
    ## $ hours                   <dbl> 0.0, 12.5, 25.5, 48.0, 71.5, 95.0, 118.5, 0.0…
    ## $ days                    <dbl> 0.0000000, 0.5208333, 1.0625000, 2.0000000, 2…
    ## $ cells                   <dbl> 1.30e+08, 1.34e+08, 1.28e+08, 1.55e+08, 1.55e…
    ## $ sd_cells                <dbl> 2.09e+07, 2.76e+07, 2.22e+07, 2.52e+07, 3.19e…

# Plot growth curves

``` r
# custom.colors <- c("Control" = "#377EB8", "Ash Leachate" = "#4DAF4A", "Santa Barbara" = "#E41A1C", "San Diego" = "#FF7F00")
levels <- c("Control", "Ash Leachate", "Mud Leachate", "UV-Mud Leachate", "Glucose_Nitrate_Phosphate", "San Diego", "Santa Barbara", "Campus Point", "Goleta Beach")

custom.colors <- c("Control" = "#377EB8", "Ash Leachate" = "#4DAF4A", "Mud Leachate" = "#E41A1C", "UV-Mud Leachate" = "#FF7F00", "Glucose_Nitrate_Phosphate" = "#313695", "Santa Barbara" = "#E41A1C", "San Diego" = "#FF7F00")

ba_curves <- cells %>% 
  filter(!Treatment == "Amendment") %>% 
  # mutate(dna = ifelse(DNA_Sample == T, "*", NA)) %>% 
  ggplot(aes(x = days, y = cells, group = interaction(Experiment, Treatment, Bottle))) +
  geom_errorbar(aes(ymin = cells - sd_cells, ymax = cells + sd_cells, color = factor(Treatment, levels = levels)), width = 0.1) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1, alpha = 0.7) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size = 3, color = "black", shape = 21) +
  # geom_text(aes(label = dna), size = 12,  color = "#E41A1C") +
  labs(x = "Days", y = expression(paste("Bacterioplankton Abundance, Cells L"^-1)), fill = "") +
  guides(color = F) +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  facet_grid(~factor(Location, levels = levels), scales = "free") +
  theme_classic2()  
```

``` r
# custom.colors <- c("Control" = "#377EB8", "Ash Leachate" = "#4DAF4A", "Santa Barbara" = "#E41A1C", "San Diego" = "#FF7F00")

ba_curves_mean <- cells %>% 
  filter(!Treatment == "Amendment") %>% 
  group_by(Experiment, Treatment, Timepoint) %>% 
  mutate(mean = mean(cells, na.rm = T),
         sd = sd(cells, na.rm = T)) %>% 
  ungroup() %>% 
  select(-c(cells, sd_cells)) %>% 
  distinct() %>% 
  ggplot(aes(x = days, y = mean, group = interaction(Experiment, Treatment, Bottle))) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = factor(Treatment, levels = levels)), width = 0.1) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1, alpha = 0.7) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size = 3, color = "black", shape = 21) +
  # geom_text(aes(label = dna), size = 12,  color = "#E41A1C") +
  labs(x = "Days", y = expression(paste("Bacterioplankton Abundance, Cells L"^-1)), fill = "") +
  guides(color = F) +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  facet_grid(~factor(Location, levels = levels), scales = "free") +
  theme_classic2()  
```

``` r
# custom.colors <- c("Control" = "#377EB8", "Ash Leachate" = "#4DAF4A", "Santa Barbara" = "#E41A1C", "San Diego" = "#FF7F00")

sd_ba_curves <- cells %>% 
  filter(!Treatment == "Amendment", Location == "San Diego") %>% 
  group_by(Experiment, Treatment, Timepoint) %>% 
  mutate(mean = mean(cells, na.rm = T),
         sd = sd(cells, na.rm = T)) %>% 
  ungroup() %>% 
  select(-c(cells, sd_cells)) %>% 
  distinct() %>% 
  ggplot(aes(x = days, y = mean, group = interaction(Experiment, Treatment, Bottle))) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = factor(Treatment, levels = levels)), width = 0.1) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1, alpha = 0.7) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size = 3, color = "black", shape = 21) +
  # geom_text(aes(label = dna), size = 12,  color = "#E41A1C") +
  labs(x = "Days", y = expression(paste("Bacterioplankton Abundance, Cells L"^-1)), fill = "") +
  guides(color = F) +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  theme_classic2()  
```

``` r
# custom.colors <- c("Control" = "#377EB8", "Ash Leachate" = "#4DAF4A", "Santa Barbara" = "#E41A1C", "San Diego" = "#FF7F00")

sb_ba_curves <- cells %>% 
  filter(!Treatment == "Amendment", Location == "Santa Barbara") %>% 
  group_by(Experiment, Treatment, Timepoint) %>% 
  mutate(mean = mean(cells, na.rm = T),
         sd = sd(cells, na.rm = T)) %>% 
  ungroup() %>% 
  select(-c(cells, sd_cells)) %>% 
  distinct() %>% 
  ggplot(aes(x = days, y = mean, group = interaction(Experiment, Treatment, Bottle))) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = factor(Treatment, levels = levels)), width = 0.1) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1, alpha = 0.7) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size = 3, color = "black", shape = 21) +
  # geom_text(aes(label = dna), size = 12,  color = "#E41A1C") +
  labs(x = "Days", y = expression(paste("Bacterioplankton Abundance, Cells L"^-1)), fill = "") +
  guides(color = F) +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  theme_classic2()  
```

``` r
# custom.colors <- c("Control" = "#377EB8", "Ash Leachate" = "#4DAF4A", "Santa Barbara" = "#E41A1C", "San Diego" = "#FF7F00")

cp_ba_curves <- cells %>% 
  filter(!Treatment == "Amendment", Location == "Campus Point") %>% 
  group_by(Experiment, Treatment, Timepoint) %>% 
  mutate(mean = mean(cells, na.rm = T),
         sd = sd(cells, na.rm = T)) %>% 
  ungroup() %>% 
  select(-c(cells, sd_cells)) %>% 
  distinct() %>% 
  ggplot(aes(x = days, y = mean, group = interaction(Experiment, Treatment, Bottle))) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = factor(Treatment, levels = levels)), width = 0.1) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1, alpha = 0.7) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size = 3, color = "black", shape = 21) +
  # geom_text(aes(label = dna), size = 12,  color = "#E41A1C") +
  labs(x = "Days", y = expression(paste("Bacterioplankton Abundance, Cells L"^-1)), fill = "") +
  guides(color = F) +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  theme_classic2()  
```

``` r
# custom.colors <- c("Control" = "#377EB8", "Ash Leachate" = "#4DAF4A", "Santa Barbara" = "#E41A1C", "San Diego" = "#FF7F00")

gb_ba_curves <- cells %>% 
  filter(!Treatment == "Amendment", Location == "Goleta Beach") %>% 
  group_by(Experiment, Treatment, Timepoint) %>% 
  mutate(mean = mean(cells, na.rm = T),
         sd = sd(cells, na.rm = T)) %>% 
  ungroup() %>% 
  select(-c(cells, sd_cells)) %>% 
  distinct() %>% 
  ggplot(aes(x = days, y = mean, group = interaction(Experiment, Treatment, Bottle))) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd, color = factor(Treatment, levels = levels)), width = 0.1) +
  geom_line(aes(color = factor(Treatment, levels = levels)), size = 1, alpha = 0.7) +
  geom_point(aes(fill = factor(Treatment, levels = levels)), size = 3, color = "black", shape = 21) +
  # geom_text(aes(label = dna), size = 12,  color = "#E41A1C") +
  labs(x = "Days", y = expression(paste("Bacterioplankton Abundance, Cells L"^-1)), fill = "") +
  guides(color = F) +
  scale_color_manual(values = custom.colors) +
  scale_fill_manual(values = custom.colors) +
  theme_classic2()  
```

``` r
saveRDS(ba_curves, "~/GITHUB/acidd_remins/output/ba_curves.rds")
saveRDS(ba_curves_mean, "~/GITHUB/acidd_remins/output/ba_curves_mean.rds")
saveRDS(sd_ba_curves, "~/GITHUB/acidd_remins/output/sd_ba_curves.rds")
saveRDS(sb_ba_curves, "~/GITHUB/acidd_remins/output/sb_ba_curves.rds")
saveRDS(cp_ba_curves, "~/GITHUB/acidd_remins/output/cp_ba_curves.rds")
saveRDS(gb_ba_curves, "~/GITHUB/acidd_remins/output/gb_ba_curves.rds")

saveRDS(cells, "~/GITHUB/acidd_remins/input/tidy_BactA_Remin_Master.rds")
```
