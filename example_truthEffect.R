
# Load Required Packages --------------------------------------------------

# library(tidyverse)
# library(BayesFactor)


# Create Help Functions ---------------------------------------------------

load_data <- function(dat_url){
  read_csv(url(dat_url)) %>% 
    dplyr::filter(phase == 2 & filter == "selected") %>% 
    dplyr::mutate(`T` = ifelse(status, "true", "false"),
                  S = as.factor(dense_rank(subject)),
                  R = ifelse(repetition == "No", "novel", "repeated"),
                  Y = trating,
                  W = groupdescription) %>% 
    dplyr::select(S, `T`, R, W, Y) %>% 
    data.frame()
}

get_bfs <- function(data, rse){
  res1 <- lmBF(Y ~ S + `T`, data = data,
               whichRandom = "S", rscaleEffects = rse[1:2])
  res2 <- lmBF(Y ~ S + `T` + R + S:R, data = data,
               whichRandom = "S", rscaleEffects = rse[1:4])
  res3 <- lmBF(Y ~ S + `T` + W, data = data,
               whichRandom = "S", rscaleEffects = rse[c(1,2,5)])
  res4 <- lmBF(Y ~ S + `T` + R + S:R + W, 
               data = data, whichRandom = "S", rscaleEffects = rse[1:5])
  res5 <- lmBF(Y ~ S + `T` + R + S:R + W +
                 R:W + S:R:W, data = data,
               whichRandom = "S", rscaleEffects = rse)
  bf1 <- extractBF(res1, logbf = T, onlybf = T)
  bf2 <- extractBF(res2, logbf = T, onlybf = T)
  bf3 <- extractBF(res3, logbf = T, onlybf = T)
  bf4 <- extractBF(res4, logbf = T, onlybf = T)
  bf5 <- extractBF(res5, logbf = T, onlybf = T)
  
  bfs <- c(
    "Y~S+T" = bf1 - bf1,
    "Y~S+T+R+SR" = bf2-bf1,
    "Y~S+T+W" = bf3-bf1,
    "Y~S+T+R+SR+W" = bf4-bf1,
    "Y~S+T+R+SR+W+RW+SRW" = bf5-bf1
  )
  return(bfs)
}


# Analyze Data ------------------------------------------------------------

# load data from OSF
data <- load_data("https://osf.io/qz2rg/download") 

# choose prior settings (see paper for details)
rse <- c(S = 2, "T" = 2, R = .5, SR = .3,
         W = .5, RW = .3, SRW = .15)

# compute log Bayes factors against null model
b <- get_bfs(data, rse)