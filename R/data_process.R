library(tidyverse)


tags <- read_csv("./data/tag_data.csv") 

tags <- tags %>% 
  select(-start_time)

write_csv(tags,"./data/tag_data.csv")

traj <- readRDS("./data/detections/ml_trajectory_list.RDS")

keep <- which(names(traj) %in% tags$tag)

##
traj_f <- traj[keep]

saveRDS(traj_f,
        "./data/detections/pair_detections.RDS")

## Remove

## AKDEs
ctmm <- readRDS("/Users/tyson/Documents/git/zebby_movement_analysis/outputs/ctmm/ctmm_fits/ctmm_fits_ml.RDS")
names <- unlist(lapply(ctmm, function(x) gsub("tag_","", x[[1]]@info$identity)))
names(ctmm) <- names
keep <- which(names(ctmm) %in% tags$tag)

ctmm_f <- ctmm[keep]
saveRDS(ctmm_f,"./data/ctmm_mods.RDS")


## AKDEs
akdes <- readRDS("/Users/tyson/Documents/git/zebby_movement_analysis/outputs/ctmm/akdes/ml_ctmm_akdes.RDS")
keep <- which(names(akdes) %in% tags$tag)

akdes_f <- akdes[keep]
saveRDS(akdes_f,"./data/akdes.RDS")

