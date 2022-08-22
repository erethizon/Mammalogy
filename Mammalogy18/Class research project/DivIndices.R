##Using package Vegan to calculate diversity indices
rm(list = ls())
library(vegan)
library(dplyr)
library(ggplot2)

#import some species data
library(readr)
Seedlings <- read.csv("~/FE2017/FinalDataAnalysis/seedling_sp.csv")
#matrix of species by site for seedling trees; numbers are the # of seedlings 
#of the species at the site
Trees <- read.csv("~/FE2017/FinalDataAnalysis/tree_sp.csv")
#matrix of species by site for mature trees; numbers are the # of trees of
#each species at each site
TS <- read.csv("~/FE2017/FinalDataAnalysis/TS_speciesrichness.csv")
#matrix of species by site, trees and seedlings combined

glimpse(Seedlings)
glimpse(Trees)
glimpse(TS)
#need to delete first row of TS
TS<-TS[-1,]
#shannonindex
SeedShannon<-diversity(Seedlings[,2:22], index = "shannon")
SeedShannon
#build some columns to put together into output
Names<-TS$Site[2:12]
Names
levels(Names)
Sites<-c("Degrasse1", "Degrasse2", "Degrasse3", "Donnerville1", "Donnerville2", 
         "Donnerville3", "Peavine1", "Peavine2", "Peavine3",
         "SouthHammond1", "SouthHammond2", "SouthHammond3")
SeedSimpson<-diversity(Seedlings[,2:22], index = "simpson")
TreesShannon<-diversity(Trees[,2:22], index = "shannon")
TreesSimpson<-diversity(Trees[,2:22], index = "simpson")
AllShannon<-diversity(TS[,2:22], index = "shannon")
AllSimpson<-diversity(TS[,2:22], index = "simpson")

DivIndices<-as.data.frame(cbind(Sites,SeedShannon,SeedSimpson,
         TreesShannon,TreesSimpson,AllShannon,AllSimpson))
write.csv(DivIndices, file="div_indices.csv")

#Now beta diversity
Beta<-betadiver(Seedlings[,2:22],method = "w" )
Beta
