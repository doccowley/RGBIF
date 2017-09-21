#Looping species and Paging within loop

rm(list=ls())

if (! requireNamespace('rgbif')) {
  stop('You need to install the rgbif package to use this function')
} else {
  library("rgbif")
}

#Configuration options
my_loops=2000
my_limit=1000
root="~/Desktop/GBIF_er403"
species_file="speciestest.txt"

#Outputs
out_folder="species_lists/plant_dist_raw/"
gbif_unknown="gbif_unknown.txt"
gbif_200k_plus="gbif_200k_plus.txt"
gbif_data_null="gbif_data_null.txt"

#Set working directory
#setwd("~/UoE_U_Drive/GBIF_auto")
setwd(root)

#Import list of species
sp_list <- read.delim(species_file, sep="\t", strip.white=TRUE)
print(sp_list)

#i=9

#Loop : species
for (i in 1:nrow(sp_list)){
  sp_name<-as.character(sp_list[i,1]) # assume name is in first column
  print(sp_name)
  
  key <- name_suggest(q=sp_name, rank='species')$key[1]

  if(is.null(key)){
    cat("Oooops!, '", sp_name, "' cannot be found!\n", sep='')
    write(sp_name,file=gbif_unknown,append=TRUE)
    #break
  } else { # Key is found
    dat <- occ_search(taxonKey=key, fields=c('name', 'key', 'speciesKey', 'decimalLatitude','decimalLongitude', 'country'), limit=1, start=1) # This potential wastes a query to gbif, but it's small ;-)
    cat("'",sp_name,"' has",dat$meta$count,"record(s) in",ceiling(dat$meta$count/my_limit),"x",my_limit,"record search(es) with (max",my_loops,"searches)\n")
    
    if(dat$meta$count > 200000) {
      print(cat('Panic!!! Too many records ... sending to',gbif_200k_plus))
      write(paste(sp_name,",",dat$meta$count,sep=''),file=gbif_200k_plus,append=TRUE)
    } else if(is.null(dat$data)) {
      print(cat('Panic!!! Data is NULL ... sending to',gbif_data_null))
      write(paste(sp_name,",",dat$meta$count,sep=''),file=gbif_data_null,append=TRUE)
    } else {
      my_start=0
      for (j in 1:my_loops){ #Loop : search records
        dat <- occ_search(taxonKey=key, fields=c('name', 'scientificName', 'key', 'publishingCountry', 'speciesKey', 'decimalLatitude','decimalLongitude', 'country'), limit=my_limit, start=my_start)
        
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
