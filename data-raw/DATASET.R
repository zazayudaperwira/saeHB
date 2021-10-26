## code to prepare `DATASET` dataset goes here

setwd("E:/1. ZAZA YUDA/BPS/Project Package/saeHBp")


dataBeta       <- read.csv("dataBeta.csv",sep=";")
dataBetaNs     <- read.csv("dataBetaNs.csv",sep=";")
dataNormal     <- read.csv("dataNormal.csv",sep=";")
dataNormalNs   <- read.csv("dataNormalNs.csv",sep=";")

usethis::use_data(dataBeta, overwrite = TRUE)
usethis::use_data(dataBetaNs, overwrite = TRUE)
usethis::use_data(dataNormal, overwrite = TRUE)
usethis::use_data(dataNormalNs, overwrite = TRUE)
