## Presentation plots

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

fg_map <- readRDS("/Users/tyson/Documents/git/zebby_movement_analysis/plots/fg_map.RDS")

## Register map API key
ggmap::register_google(key = "AIzaSyCjJjtwV1nnlkRyq_dAUACftUqGdApvchg")

## Plot home ranges #########

## Get grid points
gps <- sf::st_read("/Users/tyson/Documents/git/zebby_movement_analysis/data/grid_points.GPX") %>% 
  select(geometry)

## Read in tag log
tags <- readr::read_csv("./data/tag_data.csv", show_col_types = FALSE)

## Get the akdes from birds in pairs  
pair_akdes <- readRDS("./data/akdes.RDS")

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

# as.numeric(sf::st_centroid( hr_polys %>% 
#                               left_join(tags) %>% 
#                               filter(hr == 0.25) %>%
#                               filter(group == x)) %>%
#              sf::st_coordinates() %>%
#              data.frame() %>%
#              summarise(mean_x = mean(X),
#                        mean_y = mean(Y)))

center <- c(141.767222,-30.949444)

## Get map of FG
fg_maps <- map2(.x = hr_polys_j %>%
                  distinct(group) %>%
                  pull(group),
                .y = list(15),
                .f = function(x,y)
                  ggmap::get_map(location = center,
                                 maptype = "satellite",
                                 zoom = y))


## Reunion spots relative to 'hangout' zones
rh <- readRDS("./outputs/hangouts_reunions_all_raw.RDS") %>% 
  sf::st_as_sf(coords = c("x", "y"), 
               crs = 4326)


## Look at proportion of reunions in the hangout zones compared to contact points (all points are contact points) 
rh_sum <- rh %>% 
  # filter(reunion == "yes") %>% 
  left_join(tags %>% 
              distinct(group,section))

## Plot overlaps by sex
hr_overlap_maps <- map2(.x = as.list(overlap_df_p %>% 
                                       filter(same_group == "Breeding pair") %>% 
                                       arrange(desc(est)) %>%
                                       distinct(ind_group) %>% 
                                       pull(ind_group)),
                        .y = fg_maps,
                        .f = function(x,y)
                          ggmap(y) +
                          
                          geom_sf(data = hr_polys_j %>% 
                                    filter(group == x),
                                  aes(color = as.factor(group),
                                      size = 4,
                                      linetype = sex),
                                  linewidth = 2,
                                  fill = NA,
                                  inherit.aes = FALSE) +
                          
                          geom_sf(data = rh_sum %>%
                                    filter(group == x),
                                  color = "white",
                                  fill = NA,
                                  inherit.aes = FALSE) +
                          
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
                          annotation_scale(location = "br", width_hint = 0.4,text_col="white",pad_y=unit(0.4, "cm")))

patchwork::wrap_plots(hr_overlap_maps)

ggsave("/Users/tyson/Documents/academia/conferences/Ethology/2024/contact_points_map.jpg",
       width = 7,
       height = 7, 
       dpi=500,
       scale = 2)

