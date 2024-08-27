## Distances between individuals
library(ctmm)
library(tidyverse)

## Read in tag log
tags <- readr::read_csv("./data/pair_tags.csv", show_col_types = FALSE) %>% 
  filter(year == "2022")

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
                     ind_sex = sex,
                     year = year),
            join_by(x == ind))  %>%
  left_join(tag_df %>%
              rename(partner_tag = tag,
                     partner_section = section,
                     partner_group = group,
                     partner_sex = sex) %>% 
              select(-year),
            join_by(y == ind)) %>%
  mutate(same_group = ifelse(ind_group == partner_group & ind_group != "solo", "Breeding pair","Neighbor"),
         same_section = ifelse(ind_section == partner_section, "yes", "no"))   %>%

  ## Keep neighboring birds to compare
  filter(same_section == "yes") %>% 
  select(-bird_band.x, -bird_band.y)

write_csv(combos, "./outputs/groups_2022.csv")

# 
# ## Set progress bar
# pb <- txtProgressBar(min = 0, max = length(combos$x), style = 3)
# 
# dists_list <- list()
# sim_dists_list <- list()
# for(r in 1:nrow(combos)){
#   
#   ## Progress bar
#   Sys.sleep(0.1)
#   setTxtProgressBar(pb, which(r == 1:nrow(combos)))
#   cat("\n")
#   rows <- combos[r,]
#   
#   tag_f = rows$ind_tag
#   partner = rows$partner_tag
#   
#   times_tag <- tracks[c(tag_f)][[1]]$t
#   times_partner <- tracks[c(partner)][[1]]$t
#   times <- intersect(times_tag,times_partner)
#   
#   cat("Calculating observed distance","\n")
#   DISTS <- distances(tracks[c(tag_f,partner)],
#                      top_mod_list[c(tag_f,partner)],
#                      t=times)
#   
#   dists_df <- as_tibble(DISTS) %>% 
#     mutate(combo = r)
#   
#   dists_list[[r]] <- dists_df
#   
#   ## Simulate distances
#   tag_time = tracks[tag_f][[1]]$t
#   partner_time = tracks[partner][[1]]$t
# 
#   sim_dists <- list()
#   for(i in 1:10){
# 
#     cat("Simulating distance, rep:", i, "\n")
# 
#     tag_sim <- simulate(top_mod_list[c(tag_f)][[1]],
#                         t=tag_time)
# 
#     partner_sim <- simulate(top_mod_list[c(partner)][[1]],
#                             t=partner_time)
# 
#     dists <- distances(list(tag_sim, partner_sim),
#                        top_mod_list[c(tag_f,partner)])
# 
#     sim_dists[[i]] <- dists
# 
#   }
# 
#   sim_dists_df <- do.call(rbind, sim_dists) %>%
#     mutate(combo = r)
# 
#   sim_dists_list[[r]] <- sim_dists_df
# 
#   saveRDS(dists_list, "./outputs/simultaneous_observed_distances_list.RDS")
#   saveRDS(sim_dists_list, "./outputs/simulated_distances_list.RDS")
# 
# }
# 
# close(pb)

## 2023 ########

## Read in tag log
tags <- readr::read_csv("./data/pair_tags.csv", show_col_types = FALSE) %>% 
  filter(year == "2023")

## Get telemetry objects
tracks <- readRDS("/Users/tyson/Documents/git/zebby_movement_analysis/outputs/ctmm/trajectory/outliers_removed/ml_trajectory_list_2023.RDS")

## Get telemetry objects
fits <- readRDS("/Users/tyson/Documents/git/zebby_movement_analysis/outputs/ctmm/ctmm_fits/ctmm_fits_ml_2023.RDS")

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
                     ind_sex = sex,
                     year = year), 
            join_by(x == ind))  %>% 
  left_join(tag_df %>% 
              rename(partner_tag = tag,
                     partner_section = section,
                     partner_group = group,
                     partner_sex = sex) %>% 
              select(-year), 
            join_by(y == ind)) %>% 
  mutate(same_group = ifelse(ind_group == partner_group & ind_group != "solo", "Breeding pair","Neighbor"),
         same_section = ifelse(ind_section == partner_section, "yes", "no"))   %>% 
  
  ## Keep neighboring birds to compare 
  filter(same_section == "yes") %>% 
  select(-bird_band.x, -bird_band.y)

write_csv(combos, "./outputs/groups_2023.csv")

## Set progress bar
pb <- txtProgressBar(min = 0, max = length(combos$x), style = 3)

dists_list <- list()
sim_dists_list <- list()
for(r in 1:nrow(combos)){
  
  ## Progress bar
  Sys.sleep(0.1)
  setTxtProgressBar(pb, which(r == 1:nrow(combos)))
  cat("\n")
  rows <- combos[r,]
  
  tag_f = rows$ind_tag
  partner = rows$partner_tag
  
  times_tag <- tracks[c(tag_f)][[1]]$t
  times_partner <- tracks[c(partner)][[1]]$t
  times <- intersect(times_tag,times_partner)
  
  cat("Calculating observed distance","\n")
  DISTS <- distances(tracks[c(tag_f,partner)],
                     top_mod_list[c(tag_f,partner)],
                     t=times)
  
  dists_df <- as_tibble(DISTS) %>% 
    mutate(combo = r)
  
  dists_list[[r]] <- dists_df
  
  ## Simulate distances
  tag_time = tracks[tag_f][[1]]$t
  partner_time = tracks[partner][[1]]$t
  
  sim_dists <- list()
  for(i in 1:10){
    
    cat("Simulating distance, rep:", i, "\n")
    
    tag_sim <- simulate(top_mod_list[c(tag_f)][[1]],
                        t=tag_time)
    
    partner_sim <- simulate(top_mod_list[c(partner)][[1]],
                            t=partner_time)
    
    dists <- distances(list(tag_sim, partner_sim),
                       top_mod_list[c(tag_f,partner)])
    
    sim_dists[[i]] <- dists
    
  }
  
  sim_dists_df <- do.call(rbind, sim_dists) %>%
    mutate(combo = r)
  
  sim_dists_list[[r]] <- sim_dists_df
  
  saveRDS(dists_list, "./outputs/simultaneous_observed_distances_list_2023.RDS")
  saveRDS(sim_dists_list, "./outputs/simulated_distances_list_2023.RDS")
  
}

close(pb)

