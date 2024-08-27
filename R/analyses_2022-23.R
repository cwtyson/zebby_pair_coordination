## Analyses

## Housekeeping ######
library(tidyverse)
library(ggbeeswarm)
library(glmmTMB)
library(DHARMa)
library(lme4)
library(lmerTest)
library(ggmap)
library(patchwork)
library(ggspatial)
library(ctmm)
library(mgcv)
library(sf)
library(spatsoc)

# ## Need to supply a Google map API key to create the maps
# source("google_key.R")
# ggmap::register_google(key = google_key)  

## Summarize detections #########

## Read in tag log
tags_p <- readr::read_csv("./data/pair_tags.csv", show_col_types = FALSE)

##.Color
set.seed(4)
color_df <- tags %>% 
  distinct(section,group) %>% 
  group_by(section) %>%
  slice_sample(prop= 1) %>% 
  mutate(groups = n(),
         color = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[seq(1, 100, length.out = groups)])

## Get tracks - 2022 
tracks_2022 <- readRDS("./data/detections/pair_tracks_2022.RDS")

## Get tracks - 2023
tracks_2023 <- readRDS("./data/detections/pair_tracks_2023.RDS")

## Combine
tracks <- c(tracks_2022,tracks_2023)

## Combine
tracks_df <- lapply(tracks, function(x) x %>%
                      data.frame(tag = gsub("tag_","", x@info$identity))) %>%
  do.call(rbind, .) %>%
  data.frame() %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(32754) %>% 
  dplyr::select(tag,dt=timestamp) %>% 
  dplyr::arrange(dt) %>% 
  mutate(x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
         y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2]) %>% 
  sf::st_drop_geometry() %>% 
  left_join(tags_p)

## Summaryof detections throughout day
tracks_df_sum <- tracks_df %>% 
  mutate(day = as.numeric(format(dt,"%j")),
         hour = as.numeric(format(dt,"%H"))) %>% 
  group_by(group,tag,day,hour) %>% 
  summarise(locs = n(),
            # tot_locs = 240,
            p_locs = locs/240) ## 240 detections possible per hour (4 per minute)


ggplot(tracks_df_sum) +
  stat_smooth(aes(x=hour,
                  y=p_locs,
                  group=tag,
                  color=tag)) +
  facet_wrap(~group)


## Timeline of tracking
tracks_df_timeline_sum <- tracks_df %>% 
  group_by(year) %>% 
  mutate(day = as.numeric(format(dt,"%j")),
         min_day = min(day),
         day_adj = day-min_day) %>% 
  distinct(year,section,group,tag,day_adj) %>% 
  arrange(year, day_adj)

ggplot(tracks_df_timeline_sum) +
  geom_line(aes(x=day_adj,
                 y=group,
                group = tag,
                 color=as.factor(group)))+
  facet_grid(year+section~.,scales="free",drop=T) +
  theme_minimal() +
  scale_color_manual(breaks = color_df$group,
                     values = color_df$color)  +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")
  

## Keep same days for each pair
tracks_df_f <- tracks_df %>% 
  group_by(group,tag) %>%
  mutate(min_dt=min(dt),
         max_dt=max(dt)) %>%
  group_by(group) %>%
  mutate(max_min_dt = max(min_dt),
         min_max_dt = min(max_dt)) %>%
  filter(dt>max_min_dt) %>%
  filter(dt<min_max_dt)

## Roiling join with partners detections to get the nearest time partner was detected
nearest_df <- data.frame()
for(g in unique(tracks_df_f$group)){
  
  # g = unique(tracks_df_f$group)[1]
  
  cat(g,"\n")
  tg <- tracks_df_f %>% 
    filter(group == g)
  
  partners <- unique(tg$tag)
  
  tgc <- tg %>% 
    group_by(tag) %>% 
    summarise(count=n()) %>% 
    arrange(desc(count))
  
  tg1 <- tg %>% 
    filter(tag == tgc$tag[1]) %>%
    transmute(tag,sex,dt,tag_dt=dt)
  
  tg2 <- tg %>% 
    filter(tag == tgc$tag[2])  %>%
    transmute(partner=tag,partner_sex=sex,dt,partner_dt=dt)
  
  
  tg1 <- data.table::data.table(tg1, key = c("dt"))
  tg2 <- data.table::data.table(tg2, key = c("dt"))
  
  ## Rolling join node log to node records
  tgp <- tg2[tg1, roll = "nearest"]
  
  ## Find the difference
  tgpp <- tgp %>% 
    rowwise() %>% 
    mutate(diff = tag_dt-partner_dt,
           diff2 = as.numeric(diff),
           diffc = abs(diff2))
  
  nearest_df <- bind_rows(nearest_df,tgpp)
  
}

## Get proportion that are less than 60
nearest_df_sum <- nearest_df %>%
  group_by(group) %>%  
  summarise(sumless60 = sum(ifelse(diffc<60,T,F)),
            count = n(),
            prop = sumless60/count)

mean(nearest_df_sum$prop)

## Look at how long between partners are detected
partner_det <- tracks_df %>% 
  arrange(group,dt) %>% 
  # slice(1:1000) %>% 
  
  ## Change in tag?
  group_by(group) %>% 
  # rowwise() %>% 
  mutate(diff_tag = ifelse(lag(tag)==tag,"no","yes")) %>% 
  
  ## Filter just different tags
  filter(diff_tag == "yes") %>% 
  
  ## Day
  mutate(day = format(dt,"%j")) %>% 
  
  ## Time between detection
  group_by(group,day) %>%
  mutate(diff_time = as.numeric(lead(dt) - dt))

quantile(partner_det$diff_time, 0.90,na.rm=T)/60
hist(partner_det$diff_time,breaks=10000,xlim=c(0,500))


## Summarize time between detections
partner_det_sum <- partner_det %>% 
  group_by(group) %>% 
  summarise(mean_dt = as.numeric(mean(diff_time,na.rm=T)))

median(partner_det_sum$mean_dt)

## Summarize detections to see how often detected at same time
tracks_df_sum <- tracks_df %>% 
  
  ## Keep same days for each pair
  group_by(group,tag) %>%
  mutate(min_dt=min(dt),
         max_dt=max(dt)) %>%
  group_by(group) %>%
  mutate(max_min_dt = max(min_dt),
         min_max_dt = min(max_dt)) %>%
  filter(dt>max_min_dt) %>%
  filter(dt<min_max_dt) %>%
  group_by(group,dt) %>% 
  summarise(int_count = n()) %>% 
  group_by(group,int_count) %>% 
  summarise(count = n()) %>% 
  group_by(group) %>% 
  mutate(tot_ints = sum(count),
         prop = count/tot_ints)

median(tracks_df_sum[tracks_df_sum$int_count == 2,]$prop)

ggplot(tracks_df_sum) +
  geom_boxplot(aes(x=as.factor(int_count),y=prop))

## Summarize detections
tracks_int <- tracks_df %>% 
  mutate(dt_r = floor_date(dt,
                           unit = "day")) %>% 
  group_by(tag,
           dt_r) %>%
  arrange(tag,dt) %>% 
  mutate(int = as.numeric(lead(dt) - dt),
         dets = n(),
         days = n_distinct(dt_r,na.rm=T))

median(tracks_int$int,na.rm=T)
mean(tracks_int$int,na.rm=T)
sd(tracks_int$int,na.rm=T)


## Total detections per individual
tracks_int <- tracks_df %>% 
  group_by(tag) %>% 
  summarise(tot_dets = n())

median(tracks_int$tot_dets,na.rm=T)
mean(tracks_int$tot_dets,na.rm=T)

rm(list = ls())

##  Distance between pairs and ridge plots of distance between pairs and neighbors #####

## Read in tag log
tags_p <- readr::read_csv("./data/pair_tags.csv", show_col_types = FALSE)

## Get tracks - 2022 
tracks_2022 <- readRDS("./data/detections/pair_tracks_2022.RDS")

## Get tracks - 2023
tracks_2023 <- readRDS("./data/detections/pair_tracks_2023.RDS")

## Combine
tracks <- c(tracks_2022,tracks_2023)

## Combine
tracks_df <- lapply(tracks, function(x) x %>%
                      data.frame(tag = gsub("tag_","", x@info$identity))) %>%
  do.call(rbind, .) %>%
  data.frame() %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(32754) %>% 
  dplyr::select(tag,dt=timestamp) %>% 
  dplyr::arrange(dt) %>% 
  mutate(x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
         y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2]) %>% 
  sf::st_drop_geometry()

## Add group information
pair_tracks <- tracks_df %>%
  left_join(tags_p) %>%
  na.omit() %>%
  ungroup() %>%
  select(year,
         x,
         y,
         dt,
         tag,
         group,
         section,
         sex) %>%
  arrange(tag,
          dt) %>%
  data.frame()

## Convert to dt
tracks_dt <- data.table::setDT(pair_tracks)

## Group times - simultaneous fixes
group_times(tracks_dt, 
            datetime = 'dt')

## Get distance between simultaneous fixes for pairs
pair_distances_dt <- edge_dist(tracks_dt,
                               threshold = 2000,
                               id = "tag",
                               timegroup = "timegroup",
                               coords = c("x","y"),
                               fillNA = FALSE,
                               returnDist = TRUE,
                               splitBy = c("year","group"))

## Rearrange IDs
dyad_id(pair_distances_dt, id1 = 'ID1', id2 = 'ID2')

## Distinct
pair_distances_dt_p <- pair_distances_dt %>% 
  distinct(year,group, timegroup, distance) %>% 
  left_join(tags_p %>% 
              select(group,
                     section) %>% 
              distinct(group, section, .keep_all = T)) %>% 
  left_join(tracks_dt %>% 
              select(group,
                     timegroup,
                     dt) %>% 
              distinct())

## For each group, what is the proportion of fixes that are within X meters
pair_dist_sum <- pair_distances_dt_p %>% 
  group_by(group, section) %>% 
  summarise(tot = n(),
            together = sum(distance < 120),
            prop = round(together/tot,2),
            median = round(median(distance))) %>% 
  arrange(desc(prop))

## Change group order
pair_distances_dt_p$section <- factor(pair_distances_dt_p$section,
                                      levels = rev(unique(pair_dist_sum$section)))


## Get distance between simultaneous fixes for birds in the same section
section_distances_dt <- edge_dist(tracks_dt,
                                  threshold = 2000,
                                  id = "tag",
                                  timegroup = "timegroup",
                                  coords = c("x","y"),
                                  fillNA = FALSE,
                                  returnDist = TRUE,
                                  splitBy = c("year","section"))

## Reorder IDs
dyad_id(section_distances_dt, id1 = 'ID1', id2 = 'ID2')

## Add group information and exclude pairs
section_distances_dt_p <- section_distances_dt %>% 
  left_join(tags_p %>% 
              select(ind = tag,
                     ind_group = group,
                     ind_section = section),
            by = c("ID1" = "ind")) %>% 
  left_join(tags_p %>% 
              select(partner = tag,
                     partner_group = group,
                     partner_section = section), 
            by = c("ID2" = "partner")) %>% 
  mutate(same_group = ifelse(ind_group == partner_group & ind_group != "solo", "yes","no"),
         same_section = ifelse(ind_section == partner_section, "yes", "no")) %>% 
  filter(same_group == "no") %>% 
  filter(partner_group != "solo" ) %>%  
  distinct(year,section,ind_group, timegroup, dyadID, distance) %>% 
  select(year,
         group = ind_group,
         section,
         distance) %>% 
  arrange(group) %>% 
  group_by(section) %>% 
  tidyr::expand(year,group, distance)

## Median distances
median(pair_distances_dt_p$distance)
median(section_distances_dt_p$distance)

## By group
pair_distances_sum <- pair_distances_dt_p %>% 
  group_by(group) %>% 
  summarise(median = median(distance))

## By section
section_distances_dt_p_sum_22_23 <- section_distances_dt_p %>% 
  group_by(section) %>% 
  summarise(median = median(distance))

## Set group levels
pair_distances_dt_p_index <- pair_distances_dt_p %>% 
  group_by(year,group) %>% 
  mutate(median = median(distance)) %>% 
  arrange(median) %>% 
  group_by(median) %>% 
  mutate(num = cur_group_id()) %>% 
  ungroup() %>% 
  mutate(num_f = factor(as.character(num),levels = unique(num), ordered =T))

## Expand distances between neighbors
section_distances_dt_p_expand <- pair_distances_dt_p_index %>% 
  distinct(year, group, section,num_f) %>% 
  select(-group) %>% 
  left_join(section_distances_dt_p,by=c("year","section")) %>% 
  left_join(pair_distances_dt_p_index %>% 
              select(year,section,group) %>% 
              distinct() %>% 
              mutate(present = "yes") ) %>% 
  na.omit()

ggplot() + 
  ggridges::geom_density_ridges(alpha = .6,
                                point_alpha = 0.3,
                                height = 1,
                                scale = 0.9,
                                data = pair_distances_dt_p_index,
                                quantile_lines = TRUE, quantiles = 2,
                                aes(x = distance, 
                                    y = num_f, 
                                    fill = as.character(group))) +
  
  ggridges::geom_density_ridges(alpha = 0.1,
                                point_alpha = 0.3,
                                linetype = 2,
                                scale = 0.9,
                                color = grey(0.5),
                                height = 1,
                                data = section_distances_dt_p_expand,
                                quantile_lines = TRUE, 
                                quantiles = 2,
                                aes(x = distance, y = num_f)) +
  scale_x_continuous(limits = c(0, 1000),
                     expand = expansion(mult = c(0.01,0.03))) +
  scale_fill_manual(breaks = as.character(color_df$group),
                    values = color_df$color,
                    guide = "none") +
  labs(x = "Distance between individuals (m)",
       y = NULL) +
  theme_classic(base_size = 14) +
  facet_grid(year+paste("Site ", section)~.,scales="free",switch = "y") +
  theme(legend.position = "none",text = element_text(size = 28),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_y_discrete(expand = c(0, 0))

ggsave("./plots/ridge_plot_all.jpg",
       width = 7,
       height = 7,
       dpi=500,
       scale = 2)

rm(list = ls())

## Home range overlap #######

## Read in tag log
tags <- readr::read_csv("./data/pair_tags.csv", show_col_types = FALSE)

## Read in overlap 2022
overlap_df_2022 <-read_csv("./outputs/hr_overlap.csv") 

## Read in overlap 2023
overlap_df_2023 <-read_csv("./outputs/hr_overlap_2023.csv")

## Combine
overlap_df <- bind_rows(overlap_df_2022,overlap_df_2023)

## Join tag data
overlap_df_p <- overlap_df %>% 
  left_join(tags %>% 
              select(tag,
                     tag_year = year,
                     ind_sex = sex,
                     ind_group = group,
                     ind_section  = section), by = c("ind" = "tag")) %>% 
  left_join(tags %>% 
              select(tag,
                     partner_year = year,
                     partner_sex = sex,
                     partner_group = group,
                     partner_section  = section), by = c("partner" = "tag")) %>% 
  na.omit() %>% 
  mutate(group = case_when(ind_sex == "Male" & partner_sex == "Male" ~ "Male",
                           ind_sex == "Male" & partner_sex == "Female" ~ "Female/Male",
                           ind_sex == "Female" & partner_sex == "Male" ~ "Female/Male",
                           ind_sex == "Female" & partner_sex == "Female" ~ "Female",
                           ind_sex == "Female" & partner_sex == "Juvenile" ~ "Female/Juvenile",
                           ind_sex == "Juvenile" & partner_sex == "Female" ~ "Female/Juvenile",
                           ind_sex == "Juvenile" & partner_sex == "Male" ~ "Male/Juvenile",
                           ind_sex == "Male" & partner_sex == "Juvenile" ~ "Male/Juvenile",
                           ind_sex == "Juvenile" & partner_sex == "Juvenile" ~ "Juvenile"),
         same_group = ifelse(ind_group == partner_group & ind_group != "solo", "Breeding pair","Neighbor"),
         same_section = ifelse(ind_section == partner_section, "yes", "no"))   %>% 
  
  ## Keep neighboring birds to compare 
  filter(same_section == "yes") %>% 
  
  ## Replace 1 with large value
  mutate(est = ifelse(est == 1,0.99999,est))

median(overlap_df_p[overlap_df_p$same_group == "Breeding pair",]$est)
median(overlap_df_p[overlap_df_p$same_group == "Neighbor",]$est)

## beta regression of BC values for breeding pairs
mod <- glmmTMB::glmmTMB(est ~ same_group  + (1|ind) + (1|partner), 
                        family=beta_family(link="logit"),
                        data = overlap_df_p)
res <- simulateResiduals(mod, plot = T)

## Model with pair type
mod_null <- glmmTMB::glmmTMB(est ~ (1|ind) + (1|partner), 
                             family=beta_family(link="logit"),
                             data = overlap_df_p)

## LRT
anova(mod,mod_null)

## Get model response to plot
rg_mod <- emmeans::ref_grid(mod, 
                            type = "response",
                            at = list(same_group = c("Breeding pair","Neighbor"))) %>% 
  data.frame()

##.Color
set.seed(4)
color_df <- tags %>% 
  distinct(section,group) %>% 
  group_by(section) %>%
  slice_sample(prop= 1) %>% 
  mutate(groups = n(),
         color = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[seq(1, 100, length.out = groups)])

my_comparisons <- list( c("Breeding pair", "Neighbor"))

hr_data <- overlap_df_p %>%
  filter(ind_group != "solo" & partner_group != "solo")

(hr_overlap <- 
    ggplot() +
    geom_violin(data= hr_data,
                aes(x = same_group, y = est, 
                    group = same_group)) +
    geom_quasirandom(alpha = 0.7, 
                     aes(size = same_group,
                         color = as.factor(ind_group),
                         shape= as.factor(ind_section),
                         x = same_group, y = est, 
                         group = same_group),
                     data = hr_data %>% 
                       filter(same_group == "Breeding pair")) +
    
    geom_quasirandom(alpha = 0.5, 
                     aes(size = same_group,
                         x = same_group, y = est, 
                         shape= as.factor(ind_section),
                         group = same_group),
                     color = "black",
                     data = hr_data %>% 
                       filter(same_group == "Neighbor")) +
    
    scale_size_manual(values = c(4, 8),
                      breaks = c("Neighbor","Breeding pair"),
                      guide = "none") +
    
    geom_point(aes(x = same_group,
                   y = response),
               size = 8,
               data = rg_mod) +
    geom_segment(aes(x = same_group,
                     xend = same_group,
                     y = response-SE,
                     yend = response+SE),
                 linewidth = 1,
                 data = rg_mod) +
    scale_color_manual(breaks = color_df$group,
                       values = color_df$color,
                       guide = "none") +
    scale_shape_manual(name = "Site",values = c(19,17,15)) +
    theme_classic(base_size = 16) +
    theme(legend.position = "bottom") +
    guides(shape = guide_legend(override.aes = list(size = 4))) +
    labs(x = NULL, y = "Home range overlap"))


## Plot home ranges #########

## Read in tag log
tags <- readr::read_csv("./data/pair_tags.csv", show_col_types = FALSE)

## Get the akdes from birds in pairs  
pair_akdes <- readRDS("./data/pair_akdes.RDS")

## Summarize
hr_sum <- lapply(pair_akdes, function(x) summary(x, units = FALSE)$CI %>% 
                   data.frame()) %>% 
  do.call(rbind,.) %>% 
  mutate(tag = gsub("tag_", "",rownames(.))) %>% 
  left_join(tags) %>% 
  ## Put in sq. km
  mutate(low = low/1000000,
         est = est/1000000,
         high = high/1000000)

## Mean values
mean(hr_sum$est)
mean(hr_sum$low)
mean(hr_sum$high)

## Get HRs
hr_percs <- rep(c(0.25,0.5,0.95), length(pair_akdes))
kdes_tag <- rep(pair_akdes, each = 3)  %>%
  map2(.x = .,
       .y = hr_percs,
       .f = ~ ctmm::SpatialPolygonsDataFrame.UD(.x,
                                                level.UD = .y,
                                                level = 0.9) %>%
         sf::st_as_sf() %>%
         sf::st_transform(4326) %>%
         transmute(tag = gsub("tag_","",.x@info$identity),
                   hr = factor(.y,
                               ordered = TRUE,
                               levels = rev(c(0.25,0.5,0.95)))))

## Combine
hr_polys <- mapedit:::combine_list_of_sf(kdes_tag)
keep <- seq(2,nrow(hr_polys),by=3)
hr_polys <- hr_polys[keep,]

## Add tag info
hr_polys_j <- hr_polys %>% 
  left_join(tags) %>% 
  filter(hr == 0.95)

## Encounter areas for each section ######
sec_b <- readRDS("outputs/pkde_B.RDS")
sec_d <- readRDS("outputs/pkde_D.RDS")
sec_e <- readRDS("outputs/pkde_E.RDS")

secs <- list(sec_b,sec_d,sec_e)

## Change to sf, getting different contours for encounters
hr_percs <- rep(c(0.95),3)
enc_secs <- map2(.x = hr_percs,
                 .y = secs,
                 .f = ~ ctmm::SpatialPolygonsDataFrame.UD(.y,
                                                          level.UD = .x,
                                                          level = 0.9) %>%
                   sf::st_as_sf() %>%
                   sf::st_transform(4326) %>%
                   transmute(level = factor(.x,
                                            ordered = TRUE,
                                            levels = rev(c(0.25,0.5,0.95)))))

## Combine
enc_secs_polys <- mapedit:::combine_list_of_sf(enc_secs) %>% 
  mutate(section = rep(c("B","D","E"), each = 3))

## Retain estimated mean values
keep <- seq(2,nrow(enc_secs_polys),by=3)
enc_secs_polys <- enc_secs_polys[keep,]

## Get separate maps of FG
centers_list <- list(c(141.765, -30.94914),
                     c(141.77125, -30.94914),
                     c(141.770, -30.948))

# centers2 <- enc_secs_polys %>% 
#   aggregate(.,
#             by = list(.$section),
#             function(x) x = x[1]) %>%
#   st_centroid() %>% 
#   mutate(x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
#          y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2]) %>% 
#   select(x,y) %>% 
#   st_drop_geometry() 
# centers_list <- split(centers2, seq(nrow(centers2)))

ggmap::register_google(key = "AIzaSyBDdg9Ll4Efy8jYujrMClHZRWq2H45C-G4")
fg_maps <- map2(.x = list(15,16,16),
                .y = centers_list,
                .f = function(x,y)
                  ggmap::get_map(location = y,
                                 maptype = "satellite",
                                 zoom = x))

# saveRDS(fg_maps, "./data/fg_maps.RDS")

# fg_maps <- readRDS("./data/fg_maps.RDS")

## Read in nest coords
nest_coords <- sf::read_sf("./data/nest_boxes_and_polygon.GPX") %>% 
  select(nest = name) %>% 
  filter(grepl("GH",nest)) %>% 
  mutate(area = substring(nest,3,3))

##.Color
set.seed(4)
color_df <- tags %>% 
  distinct(section,group) %>% 
  group_by(section) %>%
  slice_sample(prop= 1) %>% 
  mutate(groups = n(),
         color = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[seq(1, 100, length.out = groups)])

## Plot home ranges
hr_overlap_maps <- map2(.x = as.list(overlap_df_p %>% 
                                       filter(same_group == "Breeding pair") %>% 
                                       arrange(desc(est)) %>%
                                       distinct(ind_section) %>% 
                                       pull(ind_section)),
                        .y = fg_maps,
                        .f = function(x,y)
                          ggmap(y) +
                          geom_sf(data = enc_secs_polys %>% 
                                    filter(section == x),
                                  color = "white",
                                  linewidth = 1,
                                  fill = "white",
                                  alpha = 0.3,
                                  inherit.aes = FALSE) +
                          geom_sf(data = hr_polys_j %>% 
                                    filter(section == x),
                                  aes(color = as.factor(group),
                                      linetype = sex),
                                  linewidth = 1,
                                  fill = NA,
                                  inherit.aes = FALSE) +
                          scale_color_manual(breaks = color_df$group,
                                             values = color_df$color) +
                          geom_sf(size = 2,
                                  color = grey(0.1),
                                  alpha = 0.9,
                                  data = nest_coords %>%
                                    dplyr::filter(area %in% c("B","D","E")),
                                  aes(shape = area),
                                  inherit.aes = FALSE) +
                          theme_minimal(base_size = 16) +
                          scale_linetype_manual(breaks = c("Male","Female"),
                                                values = c(1,4)) +
                          theme(axis.text = element_blank(),
                                axis.ticks = element_blank(),
                                legend.position = "none",
                                strip.text = element_blank(),
                                plot.margin = margin(2,2,2,2,"pt"),
                                plot.title = element_text(hjust = 0.01, vjust = -7, color="white")) +
                          labs(x = NULL, y= NULL,title = paste0("Site ", x)) +
                          annotation_scale(location = "br", width_hint = 0.4,text_col="white",pad_y=unit(0.4, "cm"),text_cex=0.9))

# hr_overlap_maps[[4]] <- hr_overlap
plot_list <- list()
plot_list[[1]] <- hr_overlap_maps[[3]]
plot_list[[2]] <- hr_overlap_maps[[1]]
plot_list[[3]] <- hr_overlap_maps[[2]]
plot_list[[4]] <- hr_overlap
patchwork::wrap_plots(plot_list, 2,2,
                      guides = "keep") +
  plot_annotation(tag_levels = 'A')

ggsave("./plots/groups_hrs_all.jpg",
       width = 7,
       height = 7,
       dpi=500,
       scale = 2)


## Proximity analysis #########

overlap_df <- read_csv("./outputs/tracking_time_overlap.csv")

## Get tracks - 2022 
tracks_2022 <- readRDS("./data/detections/pair_tracks_2022.RDS")

## Get tracks - 2023
tracks_2023 <- readRDS("./data/detections/pair_tracks_2023.RDS")

## Combine
tracks <- c(tracks_2022,tracks_2023)

## Change to data frame
tracks_df <- lapply(tracks, function(x) x %>%
                      data.frame(tag = gsub("tag_","", x@info$identity))) %>%
  do.call(rbind, .) %>%
  data.frame() %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), 
               crs = 4326) %>%
  sf::st_transform(3308) %>% 
  dplyr::select(tag,dt=timestamp) %>% 
  dplyr::arrange(tag, dt) %>% 
  mutate(x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
         y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2]) %>% 
  sf::st_drop_geometry() %>% 
  left_join(tags %>% 
              select(tag,
                     group,
                     sex)) 

## Read in proximity estimates between pairs and neighbors - 2022
prox_df_22 <- read_csv("./outputs/proximity_df.csv")

prox_df_23 <- vroom::vroom(list.files("/Users/tyson/Documents/git/zebby_movement_analysis/outputs/ctmm/proximity/2023/",
                                      full.names = T,
                                      recursive = T)) 

## Read in tag log
tags <- readr::read_csv("./data/pair_tags.csv", show_col_types = FALSE)

## Proximity
prox_df <- bind_rows(prox_df_22,prox_df_23) %>% 
  left_join(tags %>% 
              select(tag,
                     year,
                     ind_sex = sex,
                     ind_group = group,
                     ind_section  = section), by = c("ind" = "tag")) %>% 
  left_join(tags %>% 
              select(tag,
                     partner_sex = sex,
                     partner_group = group,
                     partner_section  = section), by = c("partner" = "tag")) %>% 
  # na.omit() %>% 
  mutate(group = case_when(ind_sex == "Male" & partner_sex == "Male" ~ "Male",
                           ind_sex == "Male" & partner_sex == "Female" ~ "Female/Male",
                           ind_sex == "Female" & partner_sex == "Male" ~ "Female/Male",
                           ind_sex == "Female" & partner_sex == "Female" ~ "Female",
                           ind_sex == "Female" & partner_sex == "Juvenile" ~ "Female/Juvenile",
                           ind_sex == "Juvenile" & partner_sex == "Female" ~ "Female/Juvenile",
                           ind_sex == "Juvenile" & partner_sex == "Male" ~ "Male/Juvenile",
                           ind_sex == "Male" & partner_sex == "Juvenile" ~ "Male/Juvenile",
                           ind_sex == "Juvenile" & partner_sex == "Juvenile" ~ "Juvenile"),
         same_group = ifelse(ind_group == partner_group & ind_group != "solo", "Breeding pair","Neighbor"),
         same_section = ifelse(ind_section == partner_section, "yes", "no")) %>% 
  na.omit() %>% 
  left_join(overlap_df)

# ## Tracking time and overlap
# ggplot(prox_df) +
#   geom_point(aes(x=track_time,y=est, color=same_group))


## Check number of groups
n_distinct(prox_df$ind_group)

## Distribution 
hist(prox_df$est,breaks=50)

## LMER of proximity
prox_mod <- lme4::lmer(est ~ same_group + ind_section + (1|ind) + (1|partner),
                       
                       data = prox_df)

influential <- cooksD[(cooksD > (3 * mean(cooksD, na.rm = TRUE)))]
influential <- as.numeric(names(influential))

## Remove outliers if any and rerun
prox_mod_no <- lme4::lmer(est ~ same_group + ind_section + (1|ind) + (1|partner), data = prox_df[-influential,],weights=track_time,)
plot(prox_mod_no)

## Compare - give very similar results
summary(prox_mod)
summary(prox_mod_no)

coef(summary(as(prox_mod_no,"merModLmerTest")))
coef(summary(as(prox_mod,"merModLmerTest")))

prox_mod_null1 <- lme4::lmer(est ~  ind_section + (1|ind) + (1|partner),
                             data = prox_df)
prox_mod_null2 <- lme4::lmer(est ~  same_group + (1|ind) + (1|partner),
                             data = prox_df)

## LRT
anova(prox_mod, prox_mod_null1)
anova(prox_mod,prox_mod_null2)

## Estimates
prox_mod_preds <- emmeans::ref_grid(prox_mod, 
                                    type = "response",
                                    at = list(same_group = c("Breeding pair","Neighbor"))) %>% 
  data.frame()

## Plot 
(prox_estimates <- ggplot() +
    
    # ## Neighbors
    # geom_violin(aes(x= same_group, y = est),
    #             data=prox_df,
    #             fill = NA) +
    # 
    geom_point(aes(x= same_group, y = est,shape=as.factor(year),
                   size=track_time),
               color = grey(0.6),
               alpha = 0.7,
               position = position_jitter(width = 0.2, height = 0, seed = 124),
               # size=2,
               data=prox_df %>% 
                 filter(same_group == "Neighbor"))+
    geom_linerange(aes(x=same_group,ymin = low, ymax = high,y=est),
                   color = grey(0.6),
                   alpha = 0.7,
                   position = position_jitter(width = 0.2, height = 0, seed = 124),
                   data=prox_df %>% 
                     filter(same_group == "Neighbor")) +
    
    ## Breeding pairs
    geom_point(aes(x= same_group, y = est,color = as.character(ind_group),shape=as.factor(year),
                   size=track_time),
               position = position_jitter(width = 0.2, height = 0, seed = 124),
               # size = 3,
               data=prox_df %>% 
                 filter(same_group == "Breeding pair"))+
    geom_linerange(aes(x=same_group,ymin = low, ymax = high,,y=est),
                   
                   position = position_jitter(width = 0.2, height = 0, seed = 124),
                   data=prox_df %>% 
                     filter(same_group == "Breeding pair")) +
    
    
    
    geom_point(aes(x= same_group,
                   y = prediction),
               size = 5,
               data = prox_mod_preds)+
    geom_linerange(aes(x = same_group,
                       ymin = prediction-SE,
                       ymax = prediction+SE),
                   linewidth = 1,
                   data = prox_mod_preds) +
    
    geom_hline(yintercept = 1,
               linetype = 2) +
    scale_color_manual(breaks = color_df$group,
                       values = color_df$color) +
    scale_y_continuous(limits = c(0,1.25)) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") +
    labs(x = NULL,
         y = "Proximity estimate") +
    facet_grid(~paste("Site", ind_section)))

## Distance visulalization ######

## Tags
tags_p <- read_csv("./data/pair_tags.csv")

## Groups
combos_2022 <- read_csv("./outputs/groups_2022.csv") %>% 
  mutate(combo = 1:n()) %>% 
  mutate(year = 2022)

combos_2023 <- read_csv("./outputs/groups_2023.csv") %>% 
  mutate(combo = 1:n())

combos <- bind_rows(combos_2022,combos_2023)

## Observed distances
obs_dist_2022 <- readRDS("./outputs/simultaneous_observed_distances_list_2022.RDS")
obs_dist_2023 <- readRDS("./outputs/simultaneous_observed_distances_list_2023.RDS")

obs_dist_df_2022 <- do.call(rbind, obs_dist_2022) %>% 
  rename(obs_est = est) %>% 
  mutate(year = 2022)

obs_dist_df_2023 <- do.call(rbind, obs_dist_2023) %>% 
  rename(obs_est = est) %>% 
  mutate(year = 2023)

obs_dist_df <- bind_rows(obs_dist_df_2022,obs_dist_df_2023)


obs_dist_sum <- obs_dist_df %>% 
  left_join(combos) %>% 
  group_by(ind_tag,partner_group, ind_section, same_group) %>% 
  summarise(mean_dist = mean(obs_est),
            median_dist = median(obs_est))

## Simulated distances
sim_dist_2022 <- readRDS("./outputs/simulated_distances_list_2022.RDS")
sim_dist_2023 <- readRDS("./outputs/simulated_distances_list_2023.RDS")

## Process
sim_dist_df_2022 <- do.call(rbind, sim_dist_2022) %>% 
  rename(sim_est = est) %>% 
  mutate(year = 2022)

sim_dist_df_2023 <- do.call(rbind, sim_dist_2023) %>% 
  rename(sim_est = est) %>% 
  mutate(year = 2023)

sim_dist_df <- bind_rows(sim_dist_df_2022,sim_dist_df_2023)

sim_dist_sum_df <- sim_dist_df  %>% 
  ## Get average estimate distance from simulations
  group_by(year,combo,t) %>%
  summarize(sim_est = mean(sim_est))

dist_df <- obs_dist_df %>% 
  left_join(sim_dist_sum_df,
            by = c("year", "combo","t")) %>% 
  left_join(combos) %>% 
  
  ## Update group numbers from tags
  select(-ind_group,
         -partner_group) %>% 
  
  left_join(tags_p %>% 
              transmute(ind_tag = tag,
                        ind_group = group),
            by = "ind_tag") %>% 
  
  left_join(tags_p %>% 
              transmute(partner_tag = tag,
                        partner_group = group),
            by = "partner_tag") %>% 
  na.omit()

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

## Data to model
dist_mod_dat <- dist_ss %>% 
  filter(same_group == "Breeding pair") %>% 
  transmute(obs_est,
            obs_est_t = log10(obs_est),
            t_ss,
            group = as.factor(ind_group))

## Hierarchical GAM with smoother for each group
hgam <- gam(obs_est_t ~ s(t_ss, k=20, m=2) + 
              s(t_ss, group, k=20, bs="fs", m=2), 
            data = dist_mod_dat, 
            family = "gaussian",
            method = "REML")

## Model checks
gam.check(hgam)
k.check(hgam)

## Summary
summary(hgam)

## Plot
(dist_plot <- ggplot() +
    stat_smooth(aes(x=t_ss/60,
                    y=obs_est,
                    group=ind_tag,
                    # linetype = as.factor(year),
                    color=as.character(ind_group)),
                method = 'gam',
                data = dist_ss %>% 
                  filter(same_group == "Breeding pair")) +
    scale_color_manual(breaks = color_df$group,
                       values = color_df$color) +
    
    stat_smooth(aes(x=t_ss/60,
                    y=obs_est,
                    group = same_group,
                    linetype = same_group),
                color = "black",
                method = 'gam',
                linewidth = 2,
                data=dist_ss) +
    facet_grid(paste("Site", ind_section)~.) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = seq(0,12,by=2)) +
    labs(x = "Hours since sunrise", y = "Separation distance (m)"))

prox_plots <- list()
prox_plots[[1]] <- prox_estimates
prox_plots[[2]] <- dist_plot

patchwork::wrap_plots(prox_plots, nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = "A")

ggsave("./plots/proximity_plot_all.jpg",
       width = 7,
       height = 7, 
       dpi=500,
       scale = 2)
