# Search, download and plot GBIF data using rgbif library
# Dr Andrew Cowley, ESI, Unversity of Exeter, UK
# (c) 2017

# More information on the rgbif package can be found at the following URLs:
# - https://www.gbif.org/tool/81747/rgbif-an-interface-to-the-gbif-api-for-the-r-statistical-programming-environment
# - https://ropensci.org/tutorials/rgbif_tutorial.html
# - https://cran.r-project.org/web/packages/rgbif/vignettes/rgbif_vignette.html 

# To Do : Handle more errors

rm(list=ls())

if (! requireNamespace('rgbif')) {
  stop('You need to install the rgbif package to use this function')
} else {
  library("rgbif")
}

#---------- Configuration options ----------
my_loops=2000
my_limit=1000
my_start=8
my_finish=9 # set to 0 for all records
root="~/Desktop/GBIF_er403" #setwd("~/UoE_U_Drive/GBIF_auto")
species_file="speciestest.txt"

out_folder="species_lists/plant_dist_raw/"
gbif_unknown="gbif_unknown.txt"
gbif_200k_plus="gbif_200k_plus.txt"
gbif_data_null="gbif_data_null.txt"
#-------------------------------------------

setwd(root) #Set working directory

sp_list <- read.delim(species_file, sep="\t", strip.white=TRUE) #Import list of species
cat("species_file has", nrow(sp_list), "record(s)\n")
if(my_finish == 0){my_finish=nrow(sp_list)}

#----- Loop : species -----
for (i in my_start:my_finish){
  sp_name<-as.character(sp_list[i,1]) # assume name is in first column
  cat("-----",sp_name," ( i =",i,") -----\n")
  
  #key <- name_suggest(q=sp_name, rank='species')$key[1] # Throws a warning if no matching species i.e. key is unitialised

  ns <- name_suggest(q=sp_name, rank='species')
  if(length(ns) == 0){ # no records
    cat("Oooops!, '", sp_name, "' cannot be found! ... writing to '", gbif_unknown, "'\n", sep='')
    write(sp_name,file=gbif_unknown,append=TRUE)
  } else { # records found
    my_key=ns$key[1] #Take the int key of the first obsveration
    dat <- occ_search(taxonKey=my_key, fields=c('name', 'key', 'speciesKey', 'decimalLatitude','decimalLongitude', 'country'), limit=1, start=1) # This potential wastes a query to gbif, but it's small ;-)
    cat(dat$meta$count,"record(s) in",ceiling(dat$meta$count/my_limit),"x",my_limit,"record search(es) with (max",my_loops,"searches)\n")
    
    if(dat$meta$count > 200000) {
      print(cat('Panic!!! Too many records ... sending to',gbif_200k_plus))
      write(paste(sp_name,",",dat$meta$count,sep=''),file=gbif_200k_plus,append=TRUE)
    } else if(is.null(dat$data)) {
      print(cat('Panic!!! Data is NULL ... sending to',gbif_data_null))
      write(paste(sp_name,",",dat$meta$count,sep=''),file=gbif_data_null,append=TRUE)
    } else {
      my_start=0
      for (j in 1:my_loops){ #Loop : search records
        dat <- occ_search(taxonKey=my_key, fields=c('name', 'scientificName', 'key', 'publishingCountry', 'speciesKey', 'decimalLatitude','decimalLongitude', 'country'), limit=my_limit, start=my_start)
        
        if(exists("mydata")) {
          mydata <- rbind(mydata, dat$data)
        } else {
          mydata <- dat$data
        }
        
        cat(".")
        
        if(isTRUE(dat$meta$endOfRecords)) break
        
        my_start=j*my_limit
      }
      
      cat("\n")
      mydata %>% print(n = 5)
      
      #Create filename
      file_plot=paste(out_folder,sp_name,'.png', sep='')
      my_plot <- gbifmap(mydata)
      png(filename=file_plot)
      plot(my_plot)
      dev.off()
      
      #save data to directory
      file_data=paste(out_folder,sp_name,'.csv', sep='')
      write.table(mydata, file_data, row.names = FALSE, sep = ",")
      
      rm(mydata)
      cat("---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|\n") # end of species marker (divided into 10s)

    }
  }
}
