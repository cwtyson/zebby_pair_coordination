## Reunion spots - When pairs separate where is the next place they reunite?

# Housekeeping ------------------------------------------------------------

library(tidyverse)
library(cluster)
library(sf)
library(ctmm)
library(spatsoc)

## Get reunion points ######

## Read in tag log
tags_p <- readr::read_csv("./data/tag_data.csv", show_col_types = FALSE)

## Get tracks
tracks <- readRDS("./data/detections/pair_detections.RDS")

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
  select(x,
         y,
         dt,
         tag,
         group,
         section,
         sex) %>%
  arrange(tag,
          dt) %>%
  data.frame()

## Use same time periods for each pair
pair_time_filter <- pair_tracks %>% 
  group_by(group, tag) %>% 
  summarise(min_date = min(dt),
            max_date = max(dt)) %>% 
  group_by(group) %>% 
  summarise(min_date = max(min_date),
            max_date = min(max_date),
            duration = max_date-min_date)

## Join and filter by dates
pair_tracks <- pair_tracks %>% 
  left_join(pair_time_filter) %>% 
  filter(dt < max_date & dt > min_date)

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
                               splitBy = c("group"))

## Reorder IDs
dyad_id(pair_distances_dt, id1 = 'ID1', id2 = 'ID2')

## Distinct
pair_distances_dt_p <- pair_distances_dt %>% 
  distinct(group, timegroup, distance) %>% 
  left_join(tags_p %>% 
              select(group,
                     section) %>% 
              distinct(group, section, .keep_all = T)) %>% 
  left_join(tracks_dt %>% 
              select(group,
                     timegroup,
                     dt,
                     x,
                     y) %>% 
              distinct(group,
                       timegroup,.keep_all = T))

## Summarise proportion of detections that are within 40 meters
dist_sum <- pair_distances_dt_p %>% 
  mutate(close = distance < 60) %>% 
  group_by(group) %>% 
  summarise(dets = n(),
            close = sum(close),
            prop = close/dets)
  


## Get the akdes from birds in pairs  
pair_akdes <- readRDS("./data/akdes.RDS")

## Get HRs
hr_percs <- rep(c(0.25), length(pair_akdes))
kdes_tag <- pair_akdes  %>%
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
core_hr_polys <- mapedit:::combine_list_of_sf(kdes_tag) %>% 
  left_join(tags_p)

## Empty data frames
hangouts_reunion_sum <- data.frame()
hangouts_reunions_all <- data.frame()

## Values to try
sep_val <- seq(0.5,0.75,by=0.05)
con_val <- seq(0.05,0.25,by=0.05)

## All values
cutoffs <- expand_grid(sep_val,con_val)

## For each set of values
for(r in 1:nrow(cutoffs)){
  
  cat(r,"\n")
  
  sep_val <- cutoffs[r,]$sep_val
  con_val <- cutoffs[r,]$con_val  
  
  ## Identify distant points
  pair_distances_dt <- pair_distances_dt_p  %>% 
    ## Add time information
    
    mutate(day = format(dt, "%j")) %>% 
    mutate(top_90 = quantile(distance, sep_val),
           top_10 = quantile(distance, con_val),
           distant = ifelse(distance > top_90, "yes", "no"),
           close = ifelse(distance < top_10, "yes", "no"),
           dc = paste(distant, close, sep = "-")) 
  
  ## Get status
  status <- pair_distances_dt %>% 
    select(group, day, distance, distant, close, dt) %>%
    group_by(group, day) %>% 
    pivot_longer(cols = c("distant", "close"),
                 names_to = "status") %>% 
    
    ## Keep 'contact' points
    filter(value == "yes") %>% 
    select(-value) %>% 
    rename(status_dt = dt) %>% 
    group_by(group, day) %>% 
    mutate(switch = ifelse(lag(status) != status, T, F)) %>% 
    filter(switch == T  | is.na(switch)) %>% 
    select(-switch) %>% 
    arrange(group, day, status_dt)
  
  ## Rolling join 
  pd_status <- pair_distances_dt %>% 
    left_join(status, join_by(group, day, closest(dt >= status_dt)))
  
  ## Identify reunion points
  pd_status <- pd_status %>% 
    group_by(group, day) %>% 
    
    ## If any points distant and changes from distant to close, then that is a reunion
    mutate(reunion = case_when(any(status == "distant") & (status == "close" & lag(status) == "distant") ~ "yes",
                               TRUE ~ NA)) %>% 
    ungroup() %>% 
    select(section,
           group,
           timegroup,
           dt,
           distance= distance.x,
           x,y,
           status,
           close,
           reunion)
  
  ## Get close points
  reunion_pts <- pd_status %>% 
    sf::st_as_sf(coords = c("x", "y"), 
                 crs = 32754) %>% 
    st_transform(4326) %>% 
    
    # Keep only close points (these are also reunion points)
    filter(close == "yes") %>%

    select(group,
           distance,
           close,
           reunion)
  
  ## Polygons based on HRs
  
  ## Get reunion points that are in hangouts
  hangouts_reunions_list <- list()
  for(gf in unique(reunion_pts$group)){
    
    # gf = unique(reunion_pts$group)[1]
    
    reunions_gf <- reunion_pts %>% 
      filter(group == gf) %>% 
      ungroup() %>% 
      mutate(pt = 1:n())
    
    recursion_polys_gf <- core_hr_polys %>% 
      filter(group == gf) 
    
    ## Join recursion polygons and identify if they overlap with reunion points
    hangouts_reunions_gf <- st_join(reunions_gf,
                                    recursion_polys_gf) %>% 
      group_by(pt) %>% 
      # mutate(cluster_group = ifelse(n_distinct(sex) == 1 & sex == "Male", "Male",
      #                               ifelse(n_distinct(sex) == 1 & sex == "Female", "Female", 
      #                                      ifelse(n_distinct(sex) == 2, "Both", NA)))) %>% 
      distinct(pt, .keep_all = T) %>% 
      mutate(hangout = ifelse(!is.na(group.y),"yes","no"),
             reunion = ifelse(is.na(reunion),"no",reunion)) %>% 
      select(group =group.x,
             hangout,
             reunion)
    
    print(ggplot(hangouts_reunions_gf) +

      geom_sf(data = recursion_polys_gf) +
      geom_sf(data = hangouts_reunions_gf,
              aes(color = hangout))) +
      facet_grid(~reunion)

    # print(table(hangouts_reunions_gf$hangout, hangouts_reunions_gf$reunion))
    # 
    hangouts_reunions_list[[gf]] <- hangouts_reunions_gf
    
  }
  
  
  hangouts_reunions <- do.call(rbind,
                               hangouts_reunions_list)
  
  ## Reformat points
  hangouts_reunions_ref <- hangouts_reunions %>% 
    mutate(x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
           y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2]) %>% 
    st_drop_geometry() %>% 
    mutate(r = r)
  
  hangouts_reunions_all <- bind_rows(hangouts_reunions_ref,
                                     hangouts_reunions_all)
  
  ## Summarise number of reunion points that were in hangout zones
  hangouts_reunion_sum_i <- hangouts_reunions  %>% 
    st_drop_geometry() %>%
  
    # mutate(cluster_group = ifelse(is.na(cluster_group),"non-hangout", cluster_group)) %>%
    group_by(group, hangout, reunion) %>%
    summarise(count = n())%>% 
    mutate(dets = sum(count),
           prop = count/dets) %>% 
    distinct() %>%  
    mutate(r = r,
           sep_quan = sep_val,
           con_quan = con_val)
  
  hangouts_reunion_sum <- rbind(hangouts_reunion_sum_i,
                                hangouts_reunion_sum)
  
}

saveRDS(hangouts_reunion_sum, "./outputs/hangout_revisits_all.RDS")
saveRDS(hangouts_reunions_all, "./outputs/hangouts_reunions_all_raw.RDS")


## Reunions and hangouts ######

## Reunion spots relative to 'hangout' zones
rh <- readRDS("./outputs/hangout_revisits_all.RDS")

## Look at proportion of reunions in the hangout zones compared to contact points (all points are contact points) 
rh_sum <- rh %>% 
  filter(reunion == "yes")


all_reunions <- readRDS("./outputs/hangouts_reunions_all_raw.RDS") %>% 
  sf::st_as_sf(coords = c("x", "y"), 
               crs = 4326) 


all_reunions_sum <- all_reunions %>% 
  st_drop_geometry() %>% 
  ## By group, by replicate: How many contacts occurred in/out of the core area vs reunions
  group_by(group,r) %>% 
  summarise(contacts = n(),
            contacts_h = sum(hangout == "yes"),
            reunions = sum(reunion == "yes"),
            reunions_h = sum(reunion == "yes" & hangout == "yes"),
            prop_contacts = contacts_h/contacts,
            prop_reunion = reunions_h/reunions) %>% 
  group_by(group) %>% 
  summarise(mean_prop_contacts = mean(prop_contacts, na.rm=T),
            mean_prop_reunion = mean(prop_reunion, na.rm=T),
            mean_reunion = mean(reunions_h, na.rm = T)) %>% 
  pivot_longer(cols = c(mean_prop_contacts,mean_prop_reunion),
               names_to = "location",
               values_to = "value")

## Anova comparing values
mod <- aov(value ~ location, data = all_reunions_sum)
anova(mod,mod1)
res <- simulateResiduals(mod, plot = T)

all_reunions_sum_j <- all_reunions_sum %>% 
  mutate(location = case_when(location  == "mean_prop_contacts" ~ "All contacts",
                              location  == "mean_prop_reunion" ~ "Reunion")) %>% 
  left_join(tags %>% 
              distinct(group, section))


## Get model response to plot
rg_mod <- emmeans::ref_grid(mod, 
                            type = "response") %>% 
  data.frame() %>% 
  mutate(location = case_when(location  == "mean_prop_contacts" ~ "All contacts",
                              location  == "mean_prop_reunion" ~ "Reunion"))

set.seed(4)
color_df <- tags %>% 
  distinct(group) %>% 
  slice_sample(prop= 1) %>% 
  mutate(color = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[seq(1, 100, length.out = 12)])

my_comparisons <- list( c("Breeding pair", "Neighbor"))



## Violin plot
ggplot(all_reunions_sum_j ) +
  geom_violin(aes(x = location,
                  y=value)) +
  geom_quasirandom(aes(x = location,
                       y=value,
                       color = as.factor(group),
                       shape= as.factor(section)),
                   size = 4,
                   alpha = 0.8,
                   data = all_reunions_sum_j) +
  
  geom_point(aes(x = location,
                 y = prediction),
             size = 5,
             data = rg_mod) +
  geom_segment(aes(x = location,
                   xend = location,
                   y = prediction-SE,
                   yend = prediction+SE),
               size = 1,
               data = rg_mod) +
  
  
  scale_color_manual(breaks = color_df$group,
                     values = color_df$color,
                     guide = "none") +
  scale_shape_manual(name = "Site",values = c(15,17,19)) +
  theme_classic(base_size = 32) +
  theme(legend.position = "bottom") +
  guides(shape = guide_legend(override.aes = list(size = 8))) +
  labs(x = NULL, y = "Proportion of points in core area")

ggsave("/Users/tyson/Documents/academia/conferences/Ethology/2024/reunion_compare.jpg",
       width = 7,
       height = 7, 
       dpi=500,
       scale = 1.2)


summary(mod)

## Read in overlap df
overlap_df <-read_csv("./outputs/hr_overlap.csv")

## Read in tag log
tags <- readr::read_csv("./data/tag_data.csv", show_col_types = FALSE)

## Join tag data
overlap_df <- overlap_df %>% 
  left_join(tags %>% 
              select(tag,
                     ind_sex = sex,
                     ind_group = group,
                     ind_section  = section), by = c("ind" = "tag")) %>% 
  left_join(tags %>% 
              select(tag,
                     partner_sex = sex,
                     partner_group = group,
                     partner_section  = section), by = c("partner" = "tag")) %>% 
  na.omit()

## Plot amount of overlap between sexes
overlap_df_p <- overlap_df %>% 
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
  filter(same_group == "Breeding pair")

## Get the akdes from birds in pairs  
pair_akdes <- readRDS("./data/akdes.RDS")

## Get HRs
hr_percs <- rep(c(0.25), length(pair_akdes))
kdes_tag <- pair_akdes  %>%
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
core_hr_polys <- mapedit:::combine_list_of_sf(kdes_tag) %>% 
  left_join(tags)

##.Color
set.seed(4)
color_df <- tags %>% 
  distinct(group) %>% 
  slice_sample(prop= 1) %>% 
  mutate(color = wesanderson::wes_palette("Zissou1", 100, type = "continuous")[seq(1, 100, length.out = 12)])


ggmap(fg_map) +
  geom_sf(data = core_hr_polys,
          aes(linetype = sex),
          linewidth = 2,
          fill = NA,
          inherit.aes = FALSE) +
  geom_sf(data= all_reunions,
          aes(color = as.factor(hangout)),
          alpha = 0.3,
          inherit.aes = FALSE) +
  theme_minimal() +
  scale_linetype_manual(breaks = c("Male","Female"),
                        values = c(1,4)) +
  theme(axis.text.x = element_text(angle = 45,vjust=-0.001),
        axis.text.y = element_text(angle = 45),
        # axis.ticks = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        plot.margin = margin(2,2,2,2,"pt"),
        plot.title = element_text(hjust = 0.01, vjust = -7, color="white")) +
  labs(x = NULL, y= NULL) +
  annotation_scale(location = "br", width_hint = 0.4,text_col="white",pad_y=unit(0.4, "cm")) +
  facet_wrap(~group)

## Show example of one pair
group_f = 5

## Combine
core_hr_polys_f <- mapedit:::combine_list_of_sf(kdes_tag) %>% 
  left_join(tags) %>% 
  filter(group == group_f)

all_reunions_f <- all_reunions %>% 
  filter(group == group_f) %>% 
  filter(r == 30) %>% 
  mutate(reunion = case_when(reunion  == "yes" ~ "Reunion",
                             reunion  == "no" ~ "Contact"))

mean_coordinates <- all_reunions_f %>% 
  transmute(x = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,1],
            y = as.matrix((sf::st_coordinates(.data$geometry)), ncol = 2)[,2]) %>% 
  st_drop_geometry() %>% 
  summarise(x = mean(x),
            y = mean(y))

reunion_map <- get_map(location = c(lon = mean_coordinates$x,
                                    lat = mean_coordinates$y),
                       zoom = 16,
                       maptype = "satellite")

ggmap(reunion_map) +
  
  geom_sf(data= all_reunions_f,
          aes(color = as.factor(hangout)),
          alpha = 0.8,
          size = 2,
          inherit.aes = FALSE) +
  geom_sf(data = core_hr_polys_f,
          aes(linetype = sex),
          color = "white",
          linewidth = 2,
          fill = NA,
          inherit.aes = FALSE) +
  theme_minimal(base_size = 32) +
  scale_color_manual(breaks = c("yes","no"),
                     values = c(wesanderson::wes_palette("Zissou1", 100, type = "continuous")[50],
                                wesanderson::wes_palette("Zissou1", 100, type = "continuous")[90])) +
  scale_linetype_manual(breaks = c("Male","Female"),
                        values = c(1,4)) +
  theme(legend.position = "none",
        # strip.text = element_blank(),
        plot.margin = margin(2,2,2,2,"pt"),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.01, vjust = -7, color="white")) +
  labs(x = NULL, y= NULL) +
  annotation_scale(location = "br", width_hint = 0.4,text_col="white",pad_y=unit(0.4, "cm")) +
  facet_wrap(~reunion)

ggsave("/Users/tyson/Documents/academia/conferences/Ethology/2024/proximity_points_plot.jpg",
       width = 7,
       height = 7, 
       dpi=500,
       scale = 2)


## All reunion points

all_reunions_plot <- all_reunions %>% 
  left_join(tags %>% 
              distinct(group,section))

## Get map of FG
fg_maps <- map2(.x = hr_polys_j %>%
                  distinct(section) %>%
                  pull(section),
                .y = list(16,16,15),
                .f = function(x,y)
                  ggmap::get_map(location = as.numeric(sf::st_centroid(hr_polys_j) %>%
                                                         filter(section == x) %>%
                                                         sf::st_coordinates() %>%
                                                         data.frame() %>%
                                                         summarise(mean_x = mean(X),
                                                                   mean_y = mean(Y))),
                                 maptype = "satellite",
                                 zoom = y))
centers <- list(c(141.771, -30.949135),
                c(141.77, -30.94914),
                c(141.768, -30.950))

fg_maps <- map2(.x = centers,
                .y = list(16,16,15),
                .f = function(x,y)
                  ggmap::get_map(location = x,
                                 maptype = "satellite",
                                 zoom = y))

ggmap(fg_maps[[1]]) +
  # geom_sf(data = enc_secs_polys %>% 
  #           filter(section == x),
  #         color = "white",
  #         linewidth = 1,
  #         fill = "white",
  #         alpha = 0.3,
  #         inherit.aes = FALSE) +
  geom_sf(data = hr_polys_j %>% 
            filter(section == "E") %>% 
            filter(group != 11),
          aes(linetype = sex),
          linewidth = 2,
          fill = NA,
          color = "white",
          inherit.aes = FALSE) +
  geom_sf(data = all_reunions_plot %>% 
            filter(section == "E") %>% 
            filter(group != 11) %>% 
            filter(reunion == "yes") %>% 
            filter(r == 30),
          size = 4,
          inherit.aes = FALSE,
          aes(color = as.factor(group))) +
  scale_color_manual(breaks = color_df$group,
                     values = color_df$color) +
  theme_minimal(base_size = 16) +
  scale_linetype_manual(breaks = c("Male","Female"),
                        values = c(1,4)) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        strip.text = element_blank(),
        plot.margin = margin(2,2,2,2,"pt"),
        plot.title = element_text(hjust = 0.01, vjust = -7, color="white")) +
  labs(x = NULL, y= NULL) +
  facet_wrap(~group) +
  annotation_scale(location = "br", width_hint = 0.4,text_col="white",pad_y=unit(0.4, "cm"))


ggsave("/Users/tyson/Documents/academia/conferences/Ethology/2024/reunion_points_map.jpg",
       width = 7,
       height = 7, 
       dpi=500,
       scale = 2)


