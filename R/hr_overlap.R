## Get the HR overlap between each pair of individuals

## Housekeeping ######
library(tidyverse)
library(ctmm)
library(sf)

## Read in AKDE list
akde_list <- readRDS("./data/akdes.RDS")

## Get overlap for each 
combos <- expand.grid(x = 1:length(akde_list),
                      y = 1:length(akde_list)) %>% 
  filter(x != y) %>% 
  rowwise() %>% 
  mutate(order = paste(min(x,y),max(x,y))) %>% 
  distinct(order,.keep_all = T) %>% 
  arrange(x,y) %>% 
  select(x,y)

## Set progress bar
pb <- txtProgressBar(min = 0, max = length(combos$x), style = 3)

overlap_df <- data.frame()
for(row in 1:nrow(combos)){
  
  ## Progress bar
  Sys.sleep(0.1)
  setTxtProgressBar(pb, which(row == 1:nrow(combos)))
  
  rows <- combos[row,]
  
  m <- ctmm::overlap(akde_list[as.numeric(rows)])
  df <- data.frame(ind = dimnames(m$CI)[[1]][1],
                   partner = dimnames(m$CI)[[1]][2],
                   low = round(m$CI[3],2),
                   est = round(m$CI[6],2),
                   high = round(m$CI[10],2))
  
  overlap_df <- rbind(overlap_df,df)
  
}

close(pb)

write_csv(overlap_df, "./outputs/hr_overlap.csv")
