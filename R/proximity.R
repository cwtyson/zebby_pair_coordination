## Distances between individuals
library(ctmm)
library(tidyverse)

## Read in tag log
tags <- readr::read_csv("./data/tag_data.csv", show_col_types = FALSE)

## Get telemetry objects
tracks <- readRDS("data/detections/pair_detections.RDS")

## Get telemetry objects
fits <- readRDS("data/ctmm_mods.RDS")

## Combine
top_mod_list <- list()
for(i in 1:length(fits)){
  
  top_mod <- fits[[i]][[1]]
  
  top_mod_list[[top_mod@info$identity]] <- top_mod
  
}

names(top_mod_list) <- gsub("tag_","",names(top_mod_list))


## Tag df
tag_df <- tags %>% 
  mutate(ind = 1:n())

## Combinations to calculate
combos <- expand.grid(x = 1:length(tracks),
                      y = 1:length(tracks)) %>% 
  filter(x != y) %>% 
  rowwise() %>% 
  mutate(order = paste(min(x,y),max(x,y))) %>% 
  distinct(order,.keep_all = T) %>% 
  arrange(x,y) %>% 
  select(x,y) %>% 
  left_join(tag_df %>% 
              rename(ind_tag = tag,
                     ind_section = section,
                     ind_group = group,
                     ind_sex = sex), 
            join_by(x == ind))  %>% 
  left_join(tag_df %>% 
              rename(partner_tag = tag,
                     partner_section = section,
                     partner_group = group,
                     partner_sex = sex), 
            join_by(y == ind)) %>% 
  mutate(same_group = ifelse(ind_group == partner_group & ind_group != "solo", "Breeding pair","Neighbor"),
         same_section = ifelse(ind_section == partner_section, "yes", "no"))   %>% 
  
  ## Keep neighboring birds to compare 
  filter(same_section == "yes")


## Set progress bar
pb <- txtProgressBar(min = 0, max = length(combos$x), style = 3)

prox_df <- data.frame()
for(r in 1:nrow(combos)){
  
  ## Progress bar
  Sys.sleep(0.1)
  setTxtProgressBar(pb, which(r == 1:nrow(combos)))
  
  rows <- combos[r,]
  
  tag_f = rows$ind_tag
  partner = rows$partner_tag

  PROXIMITY <- proximity(CTMM = top_mod_list[c(tag_f,partner)],
                         data = tracks[c(tag_f,partner)],
                         GUESS = ctmm(error=TRUE)) 
  
  ## Proximity of pair
  prox_pair_df <- PROXIMITY %>%
    t() %>%
    data.frame() %>%
    dplyr::mutate(ind = tag_f,
                  partner = partner)
  
  prox_df <- bind_rows(prox_df, prox_pair_df)
  
}

close(pb)

write_csv(prox_df, "./outputs/proximity_df.csv")
