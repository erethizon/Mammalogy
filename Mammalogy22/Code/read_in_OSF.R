##write a function called read_in_OSF to make a temporary space to hold a file downloaded from the Open Science Framework (OSF) and then read the file into R.

##Dependency is httr package. It must be installed

read_in_OSF <- function(fileURL,filename, UserName, Password){
  require("httr")
  file_location<-tempdir()
  getOSF <- GET(
    paste0(fileURL, "/?action=download"),
    write_disk(paste0(file_location,"/",filename), overwrite = TRUE),
    authenticate(UserName, PW))
    OSF_data<-read.csv(paste0(file_location,"/", filename))
    return(OSF_data)
}