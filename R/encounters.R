## Encounters in sub-colonies

## Housekeeping ######
library(ctmm)

## Read in data ######

## Read in AKDE list
akde_list <- readRDS("./data/akdes.RDS")

## Read in tag log
tags <- readr::read_csv("./data/tag_data.csv", show_col_types = FALSE)

## For each sub-colony, get encounter area:
for(sc in unique(tags$section)){
  
  ## Tags
  tags_sc = tags[tags$section == sc,]$tag
  
  ## Keep AKDEs from that section
  akde_sc <- akde_list[names(akde_list) %in% tags_sc]
  
  ## Get encounters
  cat("Get encounters", sc)
  enc <- ctmm::encounter(akde_sc, normalize = TRUE)
  
  ## Save
  saveRDS(enc, paste0("./outputs/encounters_",
                      sc,
                      ".RDS"))
}
