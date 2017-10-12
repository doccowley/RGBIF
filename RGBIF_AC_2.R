# Search, download and plot GBIF data using rgbif library
# Dr Andrew Cowley, ESI, Unversity of Exeter, UK
# (c) 2017

# More information on the rgbif package can be found at the following URLs:
# - https://www.gbif.org/tool/81747/rgbif-an-interface-to-the-gbif-api-for-the-r-statistical-programming-environment
# - https://ropensci.org/tutorials/rgbif_tutorial.html
# - https://cran.r-project.org/web/packages/rgbif/vignettes/rgbif_vignette.html 

# To Do : Handle more errors
{
  rm(list=ls())
  
  if (! requireNamespace('rgbif')) {
    stop('You need to install the rgbif package to use this function')
  } else {
    library("rgbif")
  }
  
  #---------- Configuration options ----------
  my_loops=2000
  my_limit=1000
  my_skip=0 # number of species to skip
  my_nrows=5 # set to -1 for all species
  my_hasCoord=TRUE
  my_hasGeoIssue=FALSE
  occ_search_fields=c('name', 'key', 'speciesKey', 'decimalLatitude','decimalLongitude', 'country')
  
  root="~/Desktop/GBIF_er403" #setwd("~/UoE_U_Drive/GBIF_auto")
  species_file="species_zhang2.txt"
  async_download=20000 # lower limit : 0 = attempt occ_download for all with at least 1 record, 200000 take over from occ_search (limited to 200k)
  search_download=20000 # upper limit : 0 = don't use occ_search, API max = 200k. After ~ 17500 records asyn_download becomes faster
  
  out_folder="data/" #"species_lists/plant_dist_raw/"
  gbif_unknown="gbif_unknown.txt"
  gbif_data_null="gbif_data_null.txt"
  
  gbif_acct=c('user', 'pass', 'email@domain') # username, passwd, email
  #-------------------------------------------
  
  async_list <- list()
  setwd(root) #Set working directory
  
  sp_list <- read.delim(species_file, sep="\t", header=FALSE,col.names=c("species","count","occ_download_key","occ_download","occ_search"),skip=my_skip,nrows=my_nrows,strip.white=TRUE) #Import list of species
  cat("Loaded", nrow(sp_list), "record(s) from species_file\n")
  
  #----- Loop : Order Species by record count -----
  cat("Checking + Sorting ...\n")
  for (i in 1:nrow(sp_list)){
    sp_name<-as.character(sp_list[i,1]) # assume name is in first column
  
    ns <- name_suggest(q=sp_name, rank='species')
    if(length(ns) == 0){ # no records
      sp_list[i,2]=0
    } else { # records found
      my_key=ns$key[1] #Take the int key of the first obsveration
      
      #search method
      dat <- occ_search(taxonKey=my_key, fields=occ_search_fields, limit=1, start=1, hasCoordinate = my_hasCoord, hasGeospatialIssue = my_hasGeoIssue) # This potential wastes a query to gbif, but it's small ;-)
      cat(i,',',sp_name,',',dat$meta$count,',',length(ns$key),',',"\n",sep='')
      sp_list[i,2]=dat$meta$count
      
      ##count method (doesn't have info on geospatial issues)
      #dat <- occ_count(taxonKey=my_key, georeferenced = TRUE)#, protocol = 'DWC_ARCHIVE')
      #cat(i,',',sp_name,',',dat,',',length(ns$key),',',"\n",sep='')
      #sp_list[i,2]=dat
    }
  }
  #sp_list <- sp_list[rev(order(sp_list$count)),c(1,2,3)] # Sor tlist by record count descending
  sp_list <- sp_list[order(sp_list$count),c(1,2,3,4,5)] # Sort list by record count ascending
  
  cat("Writing to file 'counted",species_file,"'\n",sep='')
  write.csv(sp_list,file=paste('counted_',species_file,sep=''))
  
  #stop('Stop')

  #----- Loop : species -----
  cat("sp_name,i,dat$meta$count,length(ns$key),j\n")
  for (i in 1:length(sp_list)){
    sp_name<-as.character(sp_list[i,1]) # assume name is in first column
    cat('---------------------------------------- ',sp_name,' ----------------------------------------\n')
    
    ns <- name_suggest(q=sp_name, rank='species')
    if(length(ns) == 0){ # no records
      cat("Oooops!, '", sp_name, "' cannot be found! ... writing to '", gbif_unknown, "'\n", sep='')
      write(sp_name,file=gbif_unknown,append=TRUE)
      cat(sp_name,',',i,',','-',',',length(ns$key),"\n",sep='')
    } else { # records found
      sp_name = paste(gsub(" ", "-", sp_name))
      my_key=ns$key[1] #Take the int key of the first obsveration
      dat <- occ_search(taxonKey=my_key, fields=occ_search_fields, limit=1, start=1, hasCoordinate = my_hasCoord, hasGeospatialIssue = my_hasGeoIssue) # This potential wastes a query to gbif, but it's small ;-)
      cat(i,',',sp_name,',',dat$meta$count,',',length(ns$key),',',"\n",sep='')
      
      if(is.null(dat$data)) {
        print(cat('Panic!!! Data is NULL ... sending to',gbif_data_null))
        write(paste(sp_name,",",dat$meta$count,sep=''),file=gbif_data_null,append=TRUE)
      } else {
        if(dat$meta$count > async_download){ #Code to run asynchronous download ... it works :-)
          cat('Occ_Download\n')
          async_file_check=Sys.glob(paste(root,'/',out_folder,gsub(" ", "-", sp_name),"*.zip",sep=''))
          if(length(async_file_check)>0){
            print("Zip file already exists ... skipping")
            sp_list[i,4] = "x"
          } else {
            gbif_dl <- occ_download(paste('taxonKey =',my_key), paste('hasCoordinate =', my_hasCoord), paste('hasGeospatialIssue =', my_hasGeoIssue), user=gbif_acct[1], pwd=gbif_acct[2], email=gbif_acct[3])
            sp_list[i,3] = gbif_dl[1]
            gbif_meta = occ_download_meta(sp_list[i,3])
            cat(sp_name, gbif_dl[1], '\n')
            
            my_time1=proc.time() # Start timer
            while(gbif_meta$status != 'SUCCEEDED') {
              Sys.sleep(5)
              gbif_meta = occ_download_meta(sp_list[i,3])
              cat('.')
            }
            cat('\n')
            my_time2 = proc.time() - my_time1 # Get elapsed time
            sp_list[i,4] = round(unname(my_time2[3]),3)
            
            gbif_dl_get <- occ_download_get(key=sp_list[i,3],path=root,overwrite=TRUE)
            file.rename(gbif_dl_get[1],paste(root,'/',out_folder,sp_name,'_',gbif_meta$key,".zip",sep=''))
            
            # Other functions        
            # odl = occ_download_list(user=gbif_user,pwd=gbif_pass)
            # occ_download_cancel(key='0002134-171002173027117',user=gbif_user,pwd=gbif_pass)
            # subset.data.frame(odl$results,select=c('key','status'))
            
            Sys.sleep(10)
          }
        } else {
          sp_list[i,4] = "n"
        }
        if(dat$meta$count <= 200000) { 
          cat('Occ_Search\n')
          
          file_data_check=Sys.glob(paste(out_folder,gsub(" ", "_", sp_name),"*.csv",sep=''))
          if(length(file_data_check)==0){
            #Create filenames
            file_data=paste(root,'/',out_folder,sp_name,'.csv', sep='')
            file_plot=paste(root,'/',out_folder,sp_name,'.png', sep='')

            cat(dat$meta$count,"record(s) in",ceiling(dat$meta$count/my_limit),"x",my_limit,"record search(es) with (max",my_loops,"searches)\n")
            
            my_time1=proc.time()
            cat("---------|---------|---------|---------|---------|---------|---------|---------|---------|---------|\n") # end of species marker (divided into 10s)
            for (j in 0:(my_loops-1)){ #Loop : search records
              my_start=j*my_limit
              dat <- occ_search(taxonKey=my_key, fields=c('name', 'scientificName', 'key', 'publishingCountry', 'speciesKey', 'decimalLatitude','decimalLongitude', 'country'), limit=my_limit, start=my_start, hasCoordinate = my_hasCoord, hasGeospatialIssue = my_hasGeoIssue)
              
              if(exists("mydata")) {
                mydata <- rbind(mydata, dat$data)
              } else {
                mydata <- dat$data
              }
              cat(".")
              
              if(isTRUE(dat$meta$endOfRecords)) break
            }
            cat("\n")
            my_time2 = proc.time() - my_time1 # Get elapsed time
            sp_list[i,5] = round(unname(my_time2[3]),3)
            
            #mydata %>% print(n = 5) # Print the first 5 records to screen
      
            #save data to directory
            write.table(mydata, file_data, row.names = FALSE, sep = ",")
    
            #Create and save plot
            my_plot <- gbifmap(mydata)
            png(filename=file_plot)
            plot(my_plot)
            dev.off()
    
            rm(mydata)
          } else {
            cat("Have already downloaded data for this species :-)\n")
            sp_list[i,5] = "n"
          }
        } else {
          cat("More than 200000 records for", sp_name, "\n")
        }
      }
    }
    #if(readline(prompt = "Continue? (y/n) ... ")=="n") { break }
  }
  cat("Updating file 'counted",species_file,"'\n",sep='')
  write.csv(sp_list,file=paste('counted_',species_file,sep=''))  
}
