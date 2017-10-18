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
  
  if (! requireNamespace('rgbif')) { stop('You need to install the rgbif package to use this function') } else { library('rgbif') }
  if (! requireNamespace('tibble')) { stop('You need to install the tibble package to use this function') } else { library('tibble') }  
  
  #---------- Configuration options ----------
  my_limit=1000
  my_loops=2000
  my_nrows=5 # set to -1 for all species
  my_skip=60 # number of species to skip
  
  out_folder="data/" #"species_lists/plant_dist_raw/"
  root="~/Desktop/GBIF_er403" #setwd("~/UoE_U_Drive/GBIF_auto")
  species_file="species_zhang2.txt" # Single column file of species to retrieve data from GBIF
  
  gbif_acct=c('user', 'pass', 'mail') # username, passwd, email
  my_hasCoord=TRUE
  my_hasGeoIssue=FALSE
  occ_async_dl_lim=20000 # lower limit : 0 = attempt occ_download for all with at least 1 record, 200000 take over from occ_search (limited to 200k)
  occ_search_dl_lim=20000 # upper limit : 0 = don't use occ_search, API max = 200k
  occ_search_fields=c('name', 'key', 'speciesKey', 'decimalLatitude','decimalLongitude', 'country')
  #-------------------------------------------
  
  setwd(root) #Set working directory
  
  sp_list <- read.delim(species_file, sep="\t", header=FALSE,col.names=c('species','i','j','ns_key','ns_can_name','ns_rank','count','occ_download_key','occ_download','occ_search'),skip=my_skip,nrows=my_nrows,strip.white=TRUE) #Import list of species
  cat("Loaded", nrow(sp_list), "record(s) from species_file\n")
  
  #----- Loop : Order Species by record count -----
  cat("Checking species ...\n")
  i=1
  while(i>0){ # Loop through a list that is growing (break if i > length of list) ... Allows for multiple taxonomy keys
    sp_name<-as.character(sp_list[i,'species']) # assume name is in first column
    ns <- name_suggest(q=sp_name)#, rank='species') #returns name usages by prefix matching against scientific name. Results ordered by relevance.
    
    cat(i, nrow(sp_list), as.character(sp_list[i,'species']), '\n')
    
    if(nrow(ns) == 0){ # no records
      sp_list[i,'count']=0
    } else { # records found
      for(j in 1:length(ns$key)){
        dat <- occ_search(taxonKey=ns$key[j], fields=occ_search_fields, limit=1, start=1, hasCoordinate = my_hasCoord, hasGeospatialIssue = my_hasGeoIssue) # This potential wastes a query to gbif, but it's small ;-)
        cat(i,j,nrow(sp_list),sp_name,ns$key[j],ns$canonicalName[j],ns$rank[j],dat$meta$count,'\n')

        if(j>1){ # Create new row
          sp_list <- add_row(sp_list, species=sp_name, i=i, j=j, ns_key=ns$key[j], ns_can_name=ns$canonicalName[j], ns_rank=ns$rank[j], count=dat$meta$count, .after = i+j-2)
        } else { # Is there a more efficient way of doing this ???
          sp_list[i,'i']=i
          sp_list[i,'j']=j
          sp_list[i,'ns_key']=ns$key[j]
          sp_list[i,'ns_can_name']=ns$canonicalName[j]
          sp_list[i,'ns_rank']=ns$rank[j]
          sp_list[i,'count']=dat$meta$count
        }
      }
    }
    i=i+j
    if(i>nrow(sp_list)) break
  }

  rm(ns,dat,my_nrows,my_skip,sp_name)
  
  cat("Writing to file 'count_",species_file,"'\n",sep='')
  write.csv(sp_list,file=paste('count_',species_file,sep=''))
  
  cat("Sorting species results ...\n")
  #sp_list <- sp_list[rev(order(sp_list$count)),c(1,2,3)] # Sort list by record count descending
  sp_list <- sp_list[order(sp_list$count),1:length(sp_list)] # Sort list by record count ascending
  sp_list <- sp_list[rev(order(sp_list$count)),1:length(sp_list)] # Sort list by record count descending
  
  cat("Writing to file 'count_sort_",species_file,"'\n",sep='')
  write.csv(sp_list,file=paste('count_sort_',species_file,sep=''))
  
  #stop('Stop')

  #----- Loop : species -----
  for (i in 1:nrow(sp_list)){
    sp_name <- gsub(" ", "-", as.character(sp_list[i,'species'])) # assume name is in first column
    my_key <- as.character(sp_list[i,"ns_key"])
    count <- as.integer(sp_list[i,"count"])
    cat('------------------------------ ',sp_name,'_',my_key,' (',count,') ------------------------------\n',sep='')
    
    if(count>0){ # Check that records are available to download
      if(count > occ_async_dl_lim){ #Code to run asynchronous download ... it works :-)
        cat('Occ_Download\n')
        
        async_file_check=Sys.glob(paste(root,'/',out_folder,sp_name,'_',my_key,'_',"*.zip",sep=''))
        if(length(async_file_check)>0){
          print("Zip file already exists ... skipping")
          sp_list[i,'occ_download'] = "skip"
        } else {
          gbif_dl <- occ_download(paste('taxonKey =',my_key), paste('hasCoordinate =', my_hasCoord), paste('hasGeospatialIssue =', my_hasGeoIssue), user=gbif_acct[1], pwd=gbif_acct[2], email=gbif_acct[3])
          sp_list[i,'occ_download_key'] = gbif_dl[1]
          gbif_meta = occ_download_meta(sp_list[i,'occ_download_key'])
          cat(sp_name, gbif_dl[1], '\n')
          
          my_time1=proc.time() # Start timer
          while(gbif_meta$status != 'SUCCEEDED') {
            Sys.sleep(5)
            gbif_meta = occ_download_meta(sp_list[i,'occ_download_key'])
            cat('.')
          }
          cat('\n')
          my_time2 = proc.time() - my_time1 # Get elapsed time
          sp_list[i,'occ_download'] = round(unname(my_time2[3]),3)
          
          gbif_dl_get <- occ_download_get(key=sp_list[i,'occ_download_key'],path=root,overwrite=TRUE)
          file.rename(gbif_dl_get[1],paste(root,'/',out_folder,sp_name,'_',my_key,'_',gbif_meta$key,".zip",sep=''))
          
          # Other functions        
          # odl = occ_download_list(user=gbif_user,pwd=gbif_pass)
          # occ_download_cancel(key='0002134-171002173027117',user=gbif_user,pwd=gbif_pass)
          # subset.data.frame(odl$results,select=c('key','status'))
          
          Sys.sleep(5)
        }
      } else {
        sp_list[i,'occ_download'] = "n"
      }      
      
      if(count <= occ_search_dl_lim) { 
        cat('Occ_Search\n')
        
        #Create filenames
        file_data=paste(root,'/',out_folder,sp_name,'_',my_key,'.csv', sep='')
        file_plot=paste(root,'/',out_folder,sp_name,'_',my_key,'.png', sep='')
        
        file_data_check=Sys.glob(file_data)
        if(length(file_data_check)==0){
          
          cat(count,"record(s) in",ceiling(count/my_limit),"x",my_limit,"record search(es) with (max",my_loops,"searches)\n")
          
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
          sp_list[i,'occ_search'] = round(unname(my_time2[3]),3)
          
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
          sp_list[i,'occ_search'] = "skip"
        }
      } else {
        cat("More than 200000 records for", sp_name, "\n")
        sp_list[i,'occ_search'] = "n"
      }
    } else {
      sp_list[i,'occ_download']='-'
      sp_list[i,'occ_search']='-'
    }
  }
  cat("Updating file 'count_sort_",species_file,"'\n",sep='')
  write.csv(sp_list,file=paste('count_sort_',species_file,sep=''))  
}
