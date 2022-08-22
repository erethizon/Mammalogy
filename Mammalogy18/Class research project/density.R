#start function here
      my_density<-function(DF, APQ, OutputName) {
         
      #User provides a DF with columns Quadrant, Species and a value for 
      #APQ (area sampled)
         
     #1) determine number of trees or seedlings per quadrant
          PerQ<-summarize(group_by(DF,Quadrant),
                         NumPerQ = length(Species))

     #2) get total density per quadrant =  #trees or seedling/area sampled;
          PerQ<-mutate(PerQ, TotDens = NumPerQ/APQ)

     #3) Calculate density for each species where density = #trees or seedlings
          #/area
          #Step A: Need number of trees/seedlings by species for each quadrant
          ImportancePerQ<-summarize(group_by(DF, Quadrant, Species),
                            NumTorS = length(TSID))
          #Step B: now get density
          ImportancePerQ<-mutate(ImportancePerQ, Density = NumTorS/APQ)

          #Step C: now add total density per Q to SpPerQ

             #use loop
             stop <-nrow(ImportancePerQ)  #sets counting var as equal to number of rows that need to be processed
             for (i in 1:stop){
   
                  if (ImportancePerQ$Quadrant[i] == "Q1"){
                       ImportancePerQ$TotDen[i]= PerQ$TotDens [which (PerQ$Quadrant == "Q1")]
                  }
                  else
                       if(ImportancePerQ$Quadrant[i] == "Q2"){
                            ImportancePerQ$TotDen[i]= PerQ$TotDens[which (PerQ$Quadrant == "Q2")]
                       }
   
                  else
                       if(ImportancePerQ$Quadrant[i] == "Q3"){
                            ImportancePerQ$TotDen[i]= PerQ$TotDens[which (PerQ$Quadrant == "Q3")]
                       }
                  else ImportancePerQ$TotDen[i] = PerQ$TotDens[which (PerQ$Quadrant == "Q4")]
                  i = i+1}

    #4) calculate relative density for each species =
     #(density for a species/total density) *100
         ImportancePerQ$RelDen<-(ImportancePerQ$Density/ImportancePerQ$TotDen)*100
    #5) Write output
         write.csv(ImportancePerQ, file = (paste(OutputName,"_density_rel_density.csv", sep = "")))
         #print output to console
         return(ImportancePerQ)
      }
      