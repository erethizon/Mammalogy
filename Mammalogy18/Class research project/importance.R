#Testing some ideas for working up importance values for forest data
library(dplyr)
library(ggplot2)
rm(list =ls())

#import data
library(readr)
DF <- read.csv("~/FE2017/FinalDataAnalysis/cleanSeedlings.csv")#for seedlings

DF <- read.csv("~/FE2017/FinalDataAnalysis/cleanTrees.csv")#for trees
glimpse(DF)

levels(DF$Species)
#now fix DBH (so that any trees that are deleted don't count against #'s above)
levels(DF$DBH) #DBH values <10 should not be there
Dump<-which(DF$DBH =="< 10" | DF$DBH =="<10")

DF<-DF[-Dump,]
DF$DBH<-as.character(DF$DBH)
DF$DBH<-as.numeric(DF$DBH)

#DENSITY AND RELATIVE DENSITY
#1) Calculate total number of trees in plot by quadrant:
     #Step 1
          #ask user for area sampled (area per quadrant, assign as APQ)
          APQ<-100 #100 for trees in a plot, 30 for seedlings in a plot
     #Step 2: determine number of trees or seedlings per quadrant
          PerP<-summarize(group_by(DF,Forest, PlotNum),
                         NumPerP = length(Species))

     #Step 3: get total density per quadrant =  #trees/area sampled;
          PerP<-mutate(PerP, TotDens = NumPerP/APQ)

     #Step 4: Calculate density for each species where density = #trees or seeds/area
          #Step 4A: Need number of trees or seedlings by species for each quadrant
          ImportancePerP<-summarize(group_by(DF, Forest, PlotNum, Species),
                            NumTorS = length(TSID))
          #Step 4B: now get density
          ImportancePerP<-mutate(ImportancePerP, Density = NumTorS/APQ)

          #Step 4C: now add total density per Q to SpPerQ

          #use loop
          i<-1
          stop <-nrow(ImportancePerP)  #sets counting var as equal to number of rows that need to be processed
          for (i in 1:stop){

               if (ImportancePerP$PlotNum[i] == 1){
                    ImportancePerP$TotDen[i]= PerP$TotDens [which (PerP$PlotNum == 1)]
               }
               else
                    if(ImportancePerP$PlotNum[i] == 2){
                         ImportancePerP$TotDen[i]= PerP$TotDens[which (PerP$PlotNum == 2)]
                    }

               else ImportancePerP$TotDen[i]= PerP$TotDens[which (PerP$PlotNum == 3)]
                    
                i = i+1
               }

    #Step 5: calculate relative density for each species =
     #(density for a species/total density) *100
         ImportancePerP$RelDen<-(ImportancePerP$Density/ImportancePerP$TotDen)*100

#DOMINANCE AND RELATIVE DOMINANCE
#1) calculate dominance for each species = (total of Basal Area)/area sampled
     #Step 1: Calculate basal area for each tree in DF
         DF$RAD<-DF$DBH/2 #determine radius from DBH
         DF$BA<-pi*(DF$RAD*DF$RAD) #calculates basal area per tree as pi*radius squared

     #Step 2: Get total basal area and dominance (BA/area sampled) per quadrant
         nextTab<-DF %>% group_by(Forest, PlotNum) %>% summarize(
               TotBA = sum(BA), 
               Dom = sum(BA)/APQ)
      
         PerP$TotBA<-nextTab$TotBA
         PerP$Dom<-nextTab$Dom
         
     #Step 3: get basal area and dominance per species
         BA_Dom_perSp<-summarize(group_by(DF, Forest, PlotNum, Species),
                                 BA = sum(BA), Dom = sum(BA)/APQ)

         ImportancePerP$BA<-BA_Dom_perSp$BA
         ImportancePerP$Dom<-BA_Dom_perSp$Dom

     #Step4: Add total dominance to ImportancePerP
         #use loop
         i<-1
         for (i in 1:stop){
              if (ImportancePerP$PlotNum[i] == 1){
                   ImportancePerP$TotDom[i]= PerP$Dom [which (PerP$PlotNum == 1)]
              }
              else
                   if(ImportancePerP$PlotNum[i] == 2){
                        ImportancePerP$TotDom[i]= PerP$Dom[which (PerP$PlotNum == 2)]
                   }

              else ImportancePerP$TotDom[i]= PerP$Dom[which (PerP$PlotNum == 3)]
                   
           i = i+1
           }

     #Step 5: calculate relative domianance = (dom for a species/total dom)*100
         ImportancePerP$RelDom<-(ImportancePerP$Dom/ImportancePerP$TotDom)*100

#FREQUENCY AND RELATIVE FREQUENCY
#Step 1: calculate frequency for each species = # plots in which species occurs/total number plots sampled
     #Step1A: determine total number of quadrants per plot
         TotNumPlots <-nrow(PerP)
     #Step1B: determine, for each species, the number of plots in which the species occurs
          SpeciesByPlot<-summarize(group_by(DF, Forest, Species, PlotNum), Plots = length(PlotNum) )
          SpeciesByPlot<-mutate(SpeciesByPlot, count = 1)
          SpByP<-summarize(group_by(SpeciesByPlot, Species), 
                           NumPlots = sum(count))
          SpByP<-mutate(SpByP, TotPlots = nrow(PerP))
     #Step 1C: Determine relative frequency per species = (freq value for a single species/totl of freq values for all spp)*100
          SpByP<-mutate(SpByP, RelFreq = (NumPlots/TotPlots)*100)
     #Step 1D: Assign relative frequency to each species in each plot in ImportancePerP         
          #use a loop
          #use if...then inside of for loop
          choices<-levels(DF$Species)
          runs<-nrow(ImportancePerP)
          numSp<-nrow(SpByP)
          i<-1
          k<-1
          ImportancePerP$Species<-as.character(ImportancePerP$Species)
          #begin main loop - step through once for each row in ImportancePerP
             for(i in 1:runs){#iterates as many times as there are rows in ImportancePerP
                for(k in 1:numSp){
                     #is it the species? Test condition:
                   Value<-ImportancePerP$Species[i] == choices[k]
            
                   if (Value){
                     ImportancePerP$RelFreq[i]= SpByP$RelFreq [which (SpByP$Species == choices[k])]
                   }
                   k <-k+1 }
                i<-i+1
                
             }
               
#CALCULATE IMPORTANCE VALUE FOR EACH SPECIES IN EACH QUADRANT
#Step 1: calculate importance value for each species = relative density + relative dominance + relative frequency
ImportancePerP<-mutate(ImportancePerP, IV = RelDen+RelDom+RelFreq)

##NOW GET AVERAGE IMPORATANCE VALUE PER SPECIES
AvgImpBySpecies<-summarize(group_by(ImportancePerP, Species),
                    AvgIV=mean(IV,na.rm=T),
                    sdIV=sd(IV,na.rm=T),
                    seIV=sd(IV,na.rm=T)/sqrt(sum(!is.na(IV)))
                    )

write.csv(ImportancePerP, file = "Trees/TreeIV/TREE_IV.csv")
ggplot(AvgImpBySpecies, aes(x = Species, y = AvgIV))+
     geom_col()+
   coord_flip()


