#prep files for analysis of importance values

rm(list = ls())
#1) Open a tree or seedling file
#2) split out into separate csv files per forest plot 

#list files in directory; save names
SeedFiles<-list.files(path = "Seedlings/", pattern = "\\.csv$")


#create a loop that opens each seedling file, determines the number of plots,
#and then splits out the data for each plot and saves as a .csv
stop<-length(SeedFiles)
i<-1
for (i in 1:stop){
   #open file
   tempDF<-read.csv(paste("Seedlings/",SeedFiles[i], sep = '')) #read in file
   tempDF$PlotNum<-as.factor(tempDF$PlotNum)
   runs<-levels(tempDF$PlotNum)
   Forests<-levels(tempDF$Forest)
   stop2<-length(runs)
   x<-1
      for (x in 1:stop2){
         tempDF2<-filter(tempDF, PlotNum == runs[x])
         name<-paste("Seedlings/ByPlot/",Forests,"_P",runs[x],".csv", sep = "")
         write.csv(tempDF2, file = name)
         x<-x+1
         rm(tempDF2)
      }
   i<-i+1
}

#create a loop that opens each tree file, determines the number of plots,
#and then splits out the data for each plot and saves as a .csv
TreeFiles<-list.files(path="Trees/", pattern = "\\.csv$")

stop<-length(TreeFiles)
i<-1
for (i in 1:stop){
   #open file
   tempDF<-read.csv(paste("Trees/",TreeFiles[i], sep = '')) #read in file
   tempDF$PlotNum<-as.factor(tempDF$PlotNum)
   runs<-levels(tempDF$PlotNum)
   Forests<-levels(tempDF$Forest)
   stop2<-length(runs)
   x<-1
   for (x in 1:stop2){
      tempDF2<-filter(tempDF, PlotNum == runs[x])
      name<-paste("Trees/ByPlot/",Forests,"_P",runs[x],".csv", sep = "")
      write.csv(tempDF2, file = name)
      x<-x+1
      rm(tempDF2)
   }
   i<-i+1
}


