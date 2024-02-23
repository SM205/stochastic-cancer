path = folder of the study containing markercount files in Predictions folder
setwd(path)
temp = list.files(pattern="*.csv")
list2env(
  lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
         function(i){read.csv(i,row.names=1,header=TRUE)}), envir = .GlobalEnv)
setwd(path+"/Predictions")
temp2 = list.files(pattern="*.csv")
list2env(
  lapply(setNames(temp2, make.names(gsub("*.csv$", "", temp2))), 
         function(i){read.csv(i, row.names =1)}), envir = .GlobalEnv)
setwd(path+"/final_data")

#common selection list
selected_cell_types <- c(List of selected cell types)

filter_sample1 <- which(sample1_table$cell_type_pred %in% selected_cell_types)
sample1_filtered <- sample1[,filter_sample1]
write.csv(sample1_filtered,"filtered_sample1.csv")

filter_sample2 <- which(sample2_table$cell_type_pred %in% selected_cell_types)
sample2_filtered <- sample2[,filter_sample2]
write.csv(sample2_filtered,"filtered_sample2.csv")

# This pattern can be continued for all samples in a study