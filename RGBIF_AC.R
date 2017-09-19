# Search, download and plot GBIF data using rgbif library
# Dr Andrew Cowley, ESI, Unversity of Exeter, UK
# (c) 2017

# More information on the rgbif package can be found at the following URLs:
# - https://www.gbif.org/tool/81747/rgbif-an-interface-to-the-gbif-api-for-the-r-statistical-programming-environment
# - https://ropensci.org/tutorials/rgbif_tutorial.html
# - https://cran.r-project.org/web/packages/rgbif/vignettes/rgbif_vignette.html 

# To Do : Handle errors

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
out_folder="species_lists/plant_dist_raw/"
gbif_unknown="gbif_unknown.txt"

#Set working directory
#setwd("~/UoE_U_Drive/GBIF_auto")
setwd(root)

#Import list of species
sp_list <- read.delim(species_file, sep="\t", strip.white=TRUE)
print(sp_list)

#Loop : species
for (i in 1:nrow(sp_list)){
  sp_name<-as.character(sp_list[i,1]) # assume name is in first column
  print(sp_name)
  
  key <- name_suggest(q=sp_name, rank='species')$key[1]
  if(is.null(key)){
    cat("Oooops!, '", sp_name, "' cannot be found!\n", sep='')
    write(sp_name,file=gbif_unknown,append=TRUE)
    #break
  } else {
    #Loop : search records
    my_start=0
    for (j in 1:my_loops){
      dat <- occ_search(taxonKey=key, fields=c('name', 'scientificName', 'key', 'publishingCountry', 'speciesKey', 'decimalLatitude','decimalLongitude', 'country'), limit=my_limit, start=my_start)
      
      if(exists("mydata")) {
        mydata <- rbind(mydata, dat$data)
      } else {
        cat(dat$meta$count, "record(s) in", ceiling(dat$meta$count/my_limit), "search(es) for this species (max", my_loops, "searches)\n")
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
