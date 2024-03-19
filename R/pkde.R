## Encounters in sub-colonies

## Housekeeping ######
library(ctmm)

## Read in data ######

## Get telemetry objects
tracks <- readRDS("data/detections/pair_detections.RDS")

## Get telemetry objects
akdes <- readRDS("data/akdes.RDS")

## Read in tag log
tags <- readr::read_csv("./data/tag_data.csv", show_col_types = FALSE)

## For each sub-colony, get encounter area:
for(sc in unique(tags$section)){
  
  ## Tags
  tags_sc = tags[tags$section == sc,]$tag
  
  tracks_sc <- tracks[names(tracks) %in% tags_sc]
  akdes_sc <- akdes[names(akdes) %in% tags_sc]
  
  
  ## Get encounters
  cat("Get pkde", sc)
  pkde <- ctmm::pkde(data = tracks_sc, 
                    UD = akdes_sc)
  
  ## Save
  saveRDS(pkde, paste0("./outputs/pkde_",
                      sc,
                      ".RDS"))
}
