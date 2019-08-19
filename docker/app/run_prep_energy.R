library(dplyr)
library(tidyr)

run_prep_energy <- function(datapath, dataset_name, pipeline, large_peak){

  # Output folder
  tmp_path <- getwd()
  output_path <- "output_prep_energy"
  if (!file.exists(paste0(tmp_path,"/", output_path))) dir.create(paste0(tmp_path,"/", output_path))
  if (!file.exists(paste0(output_path, "/", dataset_name))) dir.create(paste0(output_path, "/", dataset_name))
  
  # Import Data
  dir_data <- paste0(datapath)
  datasets <- dir(dir_data)
  file_format <- "vertical"
  
  data <- data.frame()
  
  if (file_format == "vertical"){
    for (d in datasets){
      temp <- read.csv(paste0(dir_data, "/", d), header = F, stringsAsFactors = F)[,1:2]
      colnames(temp)[1] <- "Date"
      temp$Date <- as.character(temp$Date)
      temp$Time <- as.numeric(substr(temp$Date,9,12))
      temp$Date <- as.Date(substr(temp$Date,1,8), format = '%Y%m%d')
      data <- rbind(data, cbind(id = d, temp %>% spread(Time,V2)))
    }
  }else{
    for (d in datasets){
      data <- rbind(data, cbind(id = d, read.csv(paste0(dir_data, "/", d), header = F, stringsAsFactors = F)))
    }
    colnames(data)[2] <- "Date"
    data$Date <- as.Date(data$Date)
  }
  
  mes_per_day <- ncol(data) - 2
  
  write.table(data, file = paste0(output_path, "/", dataset_name,"/all_data_", dataset_name, ".txt"), sep = "\t", row.names = F)
  ## Preprocessing 
  
  
  # Negative values
  if (1 %in% pipeline){
    for (i in 3:ncol(data)){
      data[which(!is.na(data[,i] < 0) & (data[,i] < 0)),i] <- data[which(!is.na(data[,i] < 0) & (data[,i] < 0)),i] * (-1)
    }
  }
  
  # Delete the days for which we have no measurements
  if (2 %in% pipeline){
    data <- data %>% filter(rowSums(is.na(data[ , 3:ncol(data)])) != mes_per_day)
  }
  
  # Delete the days that have only one load value
  if (3 %in% pipeline){
    data <- data %>% filter(apply(data[ , 3:ncol(data)], 1, function(x)length(unique(x))) != mes_per_day)
  }
  
  # Replace Nas
  na_indexes <- which(rowSums(is.na(data[ , 3:ncol(data)])) > 0)
  
  # Upward or downward shift for the na_indexes measurements
  for (i in na_indexes){
    correct <- (data %>% filter(id == data$id[i]) %>% filter(Date > data$Date[i]))
    
    if (nrow(correct) != 0){
      correct <- (correct %>% filter(!(row.names(correct) %in% as.character(na_indexes)) ))[1,which(is.na(data[i,]))]
    }else{
      correct <- (data %>% filter(id == data$id[i]) %>% filter(Date < data$Date[i]))
      correct <- (correct %>% filter(!(row.names(correct) %in% as.character(na_indexes)) ))[,which(is.na(data[i,]))]
      correct <- correct[nrow(correct),]
    }
    data[i, which(is.na(data[i,]))] <- correct
  }
  
  # Find which days are not close to the min,max,average
  if (4 %in% pipeline){
    data <- data %>% mutate(min = apply(data[ , 3:ncol(data)], 1, function(x)min(x, na.rm = T)))
    data <- data %>% mutate(max = apply(data[ , 3:ncol(data)], 1, function(x)max(x, na.rm = T)))
    data <- data %>% mutate(mean = apply(data[ , 3:ncol(data)], 1, function(x)mean(x, na.rm = T)))
    
    wrong_index_min <- which(abs(data$min - mean(data$min)) > mean(data$min))
    wrong_index_max <- which(abs(data$max - mean(data$max)) > mean(data$max))
    wrong_index_mean <- which(abs(data$mean - mean(data$mean)) > mean(data$mean) * 0.7)
    
    # Upward or downward shift for the wrong_index_min measurements
    for (i in wrong_index_min){
      diff <- (data %>% filter(id == data$id[i]) %>% filter(Date > data$Date[i]))
      diff <- (diff %>% filter(!(row.names(diff) %in% as.character(wrong_index_min)) ))$min[1] - data$min[i]
      data[i, 3:(mes_per_day + 2)] <- data[i, 3:(mes_per_day + 2)] - diff
    }
    
    # Correct the wrong_index_max measurements replacing them with the correct
    # measurements that occur after them or before them.
    for (i in wrong_index_max){
      if (i %in% wrong_index_mean){
        to_replace <- data %>% filter(id == data$id[i]) %>% filter(Date > data$Date[i])
        to_replace <- to_replace %>% filter(!(row.names(to_replace) %in% as.character(wrong_index_max)))
        if (nrow(to_replace) == 0){
          to_replace <- data %>% filter(id == data$id[i]) %>% filter(Date < data$Date[i])
          to_replace <- to_replace %>% filter(!(row.names(to_replace) %in% as.character(wrong_index_max)))
        }
        data[i, 3:(mes_per_day + 2)] <- to_replace[1, 3:(mes_per_day + 2)]
      }
    }
  }
  
  
  # Continues same values
  if (5 %in% pipeline){
    index_continues_values <- which(apply(data[ , 3:ncol(data)], 1, function(x)length(unique(x))) < (1+0.001*mes_per_day))
    
    for (i in index_continues_values){
      count <- 0
      cont_values_id <- c()
      for (j in 2:mes_per_day){
        if (!is.na(data[i, j+2]) & !is.na(data[i, j+2-1]) & (data[i, j+2] == data[i, j+2-1])){
          count <- count + 1
          cont_values_id <- c(cont_values_id, j+2)
        }else{
          count <- 0
          cont_values_id <- c()
        }
      }
      if (count > 0.75*mes_per_day){
        to_replace_prev <- (data %>% filter(id == data$id[i]) %>% filter(Date < data$Date[i])) 
        to_replace_prev <- to_replace_prev %>% filter(!(row.names(to_replace_prev) %in% as.character(index_continues_values)))
        if (nrow(to_replace_prev) > 0) to_replace_prev <- to_replace_prev[nrow(to_replace_prev),]
        to_replace_after <- (data %>% filter(id == data$id[i]) %>% filter(Date > data$Date[i]))
        to_replace_after <- to_replace_after %>% filter(!(row.names(to_replace_after) %in% as.character(index_continues_values)))
        if (nrow(to_replace_after) > 0) to_replace_after <- to_replace_after[1,]
        if (nrow(to_replace_prev) > 0 & nrow(to_replace_after) > 0){
          data[i, cont_values_id] <- (to_replace_prev[1, cont_values_id] + to_replace_after[1, cont_values_id])/2
        }else if (nrow(to_replace_after) > 0){
          data[i, cont_values_id] <- to_replace_after[1, cont_values_id]
        }else if (nrow(to_replace_prev) > 0){
          data[i, cont_values_id] <- to_replace_prev[1, cont_values_id]
        }
      }
    }
  }
  
  
  
  # Find big peaks
  if (6 %in% pipeline){
    big_peak_id <- which(data[,3:(mes_per_day+2)] > 100000, arr.ind = T)
    if (nrow(big_peak_id) > 0){
      for (i in 1:nrow(big_peak_id)){
        if ((big_peak_id[i,2] < mes_per_day) & big_peak_id[i,2] > 1) {
          data[,3:(mes_per_day+2)][big_peak_id[i,1],big_peak_id[i,2]] <- (data[,3:(mes_per_day+2)][big_peak_id[i,1],(big_peak_id[i,2] + 1)] +
                                                                            data[,3:(mes_per_day+2)][big_peak_id[i,1],(big_peak_id[i,2] + 1)])/2 
        }else if (big_peak_id[i,2] < mes_per_day) {
          data[,3:(mes_per_day+2)][big_peak_id[i,1],big_peak_id[i,2]] <- data[,3:(mes_per_day+2)][big_peak_id[i,1],(big_peak_id[i,2] + 1)]
        }else if (big_peak_id[i,2] > 1) {
          data[,3:(mes_per_day+2)][big_peak_id[i,1],big_peak_id[i,2]] <- data[,3:(mes_per_day+2)][big_peak_id[i,1],(big_peak_id[i,2] - 1)]
        }
      }
    }
  }
  
  
  # Replace Nas
  na_indexes <- which(rowSums(is.na(data[ , 3:ncol(data)])) > 0)
  
  # Upward or downward shift for the na_indexes measurements
  for (i in na_indexes){
    correct <- (data %>% filter(id == data$id[i]) %>% filter(Date > data$Date[i]))
    
    if (nrow(correct) != 0){
      correct <- (correct %>% filter(!(row.names(correct) %in% as.character(na_indexes)) ))[1,which(is.na(data[i,]))]
    }else{
      correct <- (data %>% filter(id == data$id[i]) %>% filter(Date < data$Date[i]))
      correct <- (correct %>% filter(!(row.names(correct) %in% as.character(na_indexes)) ))[,which(is.na(data[i,]))]
      correct <- correct[nrow(correct),]
    }
    data[i, which(is.na(data[i,]))] <- correct
  }
  
  write.table(data, file = paste0(output_path, "/", dataset_name,"/data_after_preprocessing", ".txt"), sep = "\t", row.names = F)
  
  # Create average plots per household 
  average_plots <- data %>%
    group_by(id) %>% 
    summarise_at(vars(colnames(data)[3:(mes_per_day + 2)]), mean)
  
  if (!file.exists(paste0(output_path, "/", dataset_name, "/", mes_per_day,"topics"))){
    dir.create(paste0(output_path, "/", dataset_name, "/", mes_per_day,"topics"))
  }
  
  write.table(average_plots[,2:ncol(average_plots)]/max(average_plots[,2:ncol(average_plots)]), file = paste0(output_path, "/", dataset_name, "/", mes_per_day,"topics","/model-final.theta"), sep = " ", row.names = F, col.names = F)
  write.table(average_plots[,1], file = paste0(output_path, "/", dataset_name, "/", mes_per_day,"topics","/docIds", ".dat"), sep = "\t", row.names = F, col.names = F)
  
}