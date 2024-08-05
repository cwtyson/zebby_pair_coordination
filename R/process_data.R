## Script to process various 'raw' datasets for use in the analysis
library(tidyverse)

## Process tag data ###########

## Tag data from field 
tags_all <- readxl::read_excel("/Users/tyson/Library/CloudStorage/GoogleDrive-cwtyson@gmail.com/My Drive/Zebby_tracking_field_data/tags/zebby_tag_log_20240805.xlsx") %>% 
  janitor::clean_names()

## AKDES 2022
pair_akdes_2022 <- readRDS("./data/akdes.RDS")

## AKDES 2023
pair_akdes_2023 <- readRDS("./data/akdes_2023.RDS")

## Combine
akdes <- c(pair_akdes_2022,pair_akdes_2023)

## Get list of tags with tracking data
tracked_tags <- names(akdes)

## Process to keep groups
tags_p <- tags_all %>% 
  mutate(partner = ifelse(partner == "NA",NA,partner),
         tag = ifelse(tag == "NA",NA,tag)) %>% 
  filter(!is.na(partner)) %>%
  filter(!is.na(tag)) %>% 
  
  ## Only tags that were tracked
  filter(tag %in% tracked_tags) %>% 
  
  filter(sex != "Juvenile") %>% 
  filter(species == "ZF") %>% 
  rowwise() %>% 
  mutate(group_bands = list(sort(c(bird_band,partner)))) %>% 
  group_by(group_bands) %>% 
  mutate(count = n()) %>% 
  filter(count == 2) %>%
  mutate(group = cur_group_id()) %>% 
  mutate(year = format(mdy(date), "%Y")) %>% 
  mutate(section = ifelse(nchar(section) > 1,
                          substring(section,9,9),
                          section)) %>% 
  ungroup() %>% 
  select(year,
         tag,
         bird_band,
         sex,
         section,
         group) %>% 
  arrange(group)

write_csv(tags_p, "./data/pair_tags.csv")

## Process AKDEs #############

## Read in all AKDEs and filter to keep pairs

## AKDES 2022
pair_akdes_2022 <- readRDS("./data/akdes.RDS")

## AKDES 2023
pair_akdes_2023 <- readRDS("./data/akdes_2023.RDS")

## Combine
akdes <- c(pair_akdes_2022,pair_akdes_2023)

## Tags
tags_p <- read_csv("./data/pair_tags.csv")

## Pair tags
pair_tags <- tags_p$tag

## Keep AKDEs from paired birds
pair_akdes <- akdes[names(akdes) %in% pair_tags]

## AKDE df
akde_df <- data.frame(akde = "yes",
                      tag = names(pair_akdes))

## Check 
tags_J <- tags_p %>% 
  left_join(akde_df)

## Save
saveRDS(pair_akdes,"./data/pair_akdes.RDS")
