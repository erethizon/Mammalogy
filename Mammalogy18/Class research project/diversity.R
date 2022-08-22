#testing vegan

library(vegan)
library(ggplot2)
library(dplyr)
rm(list =ls())


#Diversity

#import data
library(readr)
SR <- read.csv("~/rstudioshared/Biol4018_F17_ForestEcology/Projects/Forest Analysis/In class/TS_speciesrichness.csv")#delete first column

SR<-SR[,-1] #dumps first column of site names

Shannon<-diversity(SR, index = "shannon")
Shannon

Simpson<-diversity(SR, index = "simpson")
Simpson

diversity(SR, index = "invsimpson")

#vector of forests
Forests<-c("Degrasse", "Degrasse", "Degrasse",
           "Donnerville", "Donnerville", "Donnerville",
           "Peavine", "Peavine", "Peavine",
           "S. Hammond", "S. Hammond", "S. Hammond")
DivIndices<-as.data.frame(cbind(Forests, Shannon, Simpson))
