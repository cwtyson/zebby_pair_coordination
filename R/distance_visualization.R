## Distances between individuals 
library(tidyversxe)

## Map distances between individuals #######

## Get telemetry objects
tracks <- readRDS("data/detections/pair_detections.RDS")

## Groups
combos <- read_csv("./outputs/groups.csv") %>% 
  mutate(combo = 1:n())

## Observed distances
obs_dist <- readRDS("./outputs/observed_distances_list.RDS")

obs_dist_df <- do.call(rbind, obs_dist) %>% 
  rename(obs_est = est)

## Map individual points and distances
for(r in 1:nrow(combos)){
  
  rows <- combos[r,]
  
  tag_f = rows$ind_tag
  partner = rows$partner_tag
  
  
  
  
  tracks_tag <- tracks[c(tag_f)][[1]]
  
  tracks_partner <- tracks[c(partner)][[1]]
  
  
  
  ## Distances
  dist_df <- obs_dist_df %>% 
    filter(combo == r)
  
  ggplot() +
    geom_point(data = tracks_tag,aes(x=x,
                              y=y)) +
    geom_point(data = tracks_partner,aes(x=x,
                                     y=y),
               color = "gold")
  
  
  
}

tracks[c(tag_f,partner)]


## Plot distances between individuals #######

## Groups
combos <- read_csv("./outputs/groups.csv") %>% 
  mutate(combo = 1:n())

## Observed distances
obs_dist <- readRDS("./outputs/simultaneous_observed_distances_list.RDS")

obs_dist_df <- do.call(rbind, obs_dist) %>% 
  rename(obs_est = est)

## Simulated distances
sim_dist <- readRDS("./outputs/simulated_distances_list.RDS")

sim_dist_df <- do.call(rbind, sim_dist) %>% 
  rename(sim_est = est) 

sim_dist_sum_df <- sim_dist_df  %>% 
  ## Get average estimate distance from simulations
  group_by(combo,t) %>%
  summarize(sim_est = mean(sim_est))

dist_df <- obs_dist_df %>% 
  left_join(sim_dist_sum_df,
            by = c("combo","t")) %>% 
  left_join(combos)

## Reformat
dist_df_f <- dist_df %>% 
  mutate(dist_diff = obs_est - sim_est,
         timestamp = as.POSIXct(timestamp,tz = "Australia/Broken_Hill"))

## Get sunrise times
ss_times <- data.frame(date = seq.Date(from =  as.Date(lubridate::floor_date(min(dist_df_f$timestamp,
                                                                                 na.rm = TRUE),
                                                                             "day")),
                                       to = as.Date(lubridate::ceiling_date(max(dist_df_f$timestamp,
                                                                                na.rm = TRUE),
                                                                            "day")),
                                       by = "day")) %>%
  with(., suncalc::getSunlightTimes(date = date,
                                    lat = -31.088747,
                                    lon = 141.684423,
                                    tz = 'Australia/Broken_Hill',
                                    keep = c('sunrise', 'sunset'))) %>%
  dplyr::select(date, sunrise) 

## Add sunrise times
dist_ss <- dist_df_f %>% 
  mutate(date = lubridate::floor_date(timestamp, unit = "day")) %>% 
  left_join(ss_times) %>% 
  
  ## Time since sunrise
  mutate(t_ss = round(as.numeric(difftime(timestamp, sunrise, unit = "mins"))))

set.seed(4)
color_df <- tags %>% 
  distinct(group) %>% 
  slice_sample(prop= 1) %>% 
  mutate(color = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[seq(1, 100, length.out = 12)])

(dist_plot <- ggplot() +
    stat_smooth(aes(x=t_ss/60,
                    y=obs_est,
                    group=ind_tag,
                    color=as.character(ind_group)),
                data = dist_ss %>% 
                  filter(same_group == "Breeding pair")) +
    # stat_smooth(aes(x=timestamp,
    #                 y=obs_est,
    #                 group=ind_tag),
    #             color = grey(0.5),
    #             alpha = 0.7,
    #             data= dist_df_f %>%
    #               filter(same_group == "Neighbor")) +
    stat_smooth(aes(x=t_ss/60,
                    y=obs_est,
                    group = same_group,
                    linetype = same_group),
                color = "black",
                size = 2,
                data=dist_ss) +
    facet_grid(ind_section~.) +
    scale_color_manual(breaks = color_df$group,
                       values = color_df$color) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none") +
    labs(x = "Hours since sunrise", y = "Distance (m) between individuals"))

