## Recursions

# Housekeeping ------------------------------------------------------------

library(recurse)
library(tidyverse)
library(ggmap)
library(factoextra)
library(cluster)
library(sf)

## Get grid points
gps <- st_read("./data/grid_points.GPX") %>% 
  select(geometry)

fg_map <- readRDS("./data/fg_map.RDS")

ss_colors <- paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging",
                                    n = 100)


# Revisits per pair -------------------------------------------------------

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
         id = tag,
         group) %>%
  arrange(id,
          dt) %>%
  data.frame()

## Use same time periods for each pair
pair_time_filter <- pair_tracks %>% 
  group_by(group, id) %>% 
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

radius <- seq(10,50,by=10)
quant <- seq(0.7,0.9,by=0.05)
combos <- expand_grid(radius,quant)

for(r in 1:nrow(combos)){
  
  cat(r,"\n")
  
  radius = combos[r,]$radius
  quant = combos[r,]$quant
  
  ## Nest and identify revisits
  pair_revisits <- pair_tracks %>%
    select(-min_date,
           -max_date,
           -duration) %>% 
    nest(.by = id) %>%
    mutate(data =  map(data, ~.x %>%
                         mutate(revisits = getRecursions(data.frame(.),
                                                         radius = radius,
                                                         threshold = 1,
                                                         timeunits = "mins")[[1]]))) %>%
    unnest(cols = c(data))
  
  ## Identify clusters
  pair_clusters <- pair_revisits %>%
    
    ## Scale revisits
    group_by(id) %>%
    mutate(revisits_scl = as.numeric(scale(revisits)),
           top = quantile(revisits_scl, quant)) %>%
    
    ## Keep most revisited
    filter(revisits_scl > top) %>%
    
    ## For each group, identify clusters
    group_by(id) %>%
    mutate(cluster = as.character(kmeans(x = as.matrix(x,y,ncol =2), centers = 2)$cluster)) %>%
    arrange(group,
            id,
            revisits) %>%
    left_join(tags_p %>%
                select(id = tag,
                       sex,
                       group,
                       section),
              by = c("group", "id")) %>% 
    sf::st_as_sf(coords = c("x", "y"), crs = 32754) %>%
    sf::st_transform(4326)
  
  ## Get data frame of sunrise and sunset times
  ss_times <- data.frame(date = seq.Date(from =  as.Date(lubridate::round_date(min(pair_clusters$dt,
                                                                                   na.rm = TRUE),
                                                                               "day")),
                                         to = as.Date(lubridate::round_date(max(pair_clusters$dt,
                                                                                na.rm = TRUE),
                                                                            "day")),
                                         by = "day")) %>%
    with(., suncalc::getSunlightTimes(date = date,
                                      lat = -31.088747,
                                      lon = 141.684423,
                                      tz = 'Australia/Broken_Hill',
                                      keep = c('sunrise'))) %>%
    mutate(yd = format(date,"%Y-%j")) %>%
    select(yd, sunrise)
  
  ## Add sunrise
  pair_clusters_f <- pair_clusters %>%
    mutate(yd = format(dt,"%Y-%j")) %>%
    left_join(ss_times,
              by = c("yd")) %>%
    mutate(post_ss = round(as.numeric(difftime(dt, sunrise, unit = "mins"))))
  
  ## Get polygons around recursion points
  recursion_polys <- pair_clusters_f %>% 
    dplyr::group_by(group, id, sex, cluster) %>%
    dplyr::summarize(geometry = sf::st_union(geometry)) %>%
    sf::st_convex_hull() 
  
  ## Save
  saveRDS(recursion_polys, paste0("./outputs/pair_revisits_",radius,"_",quant,".RDS"))
  
  
}
