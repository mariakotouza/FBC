#load_rData=F
#save_rData=F

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

load_rData <- F

if (load_rData){
  x <- loadRData(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/x","_thrBranch",thr_branch_break,".rData"))
  df <- loadRData(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/df","_thrBranch",thr_branch_break,".rData"))
  Clus <- loadRData(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/Clus","_thrBranch",thr_branch_break,".rData"))
  clep <- loadRData(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/clep","_thrBranch",thr_branch_break,".rData"))
  xN <- loadRData(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/xN","_thrBranch",thr_branch_break,".rData"))
  listb <- loadRData(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/listb",".rData"))
  num_of_levels=ncol(x[colSums(!is.na(x)) > 0])-length(which(!str_detect(names(x),"level")))-1
  flagtic = F
  
  lastlist = listb
  # A list with the permat matrix for every branch
  perlist = lastlist$list
  persimlist = lastlist$listn
}

######################### Distances #########################
#use_only_leafs=T
print("Start with distances!!!!")
lev=num_of_levels
if (use_only_leafs){
  t1 = which(clep == lev)
  # epipleon
  xm = as.numeric(as.data.frame(xN$leaves))
  orio = min(which(clep == lev))
  xm2 = sort(xm[xm < orio])
  useful_nodes = sort(append(xm2,t1,after = length(xm2)))
  bo=length(useful_nodes)
}else{
  bo = length(which(!is.na(clep)))
  useful_nodes=which(!is.na(clep))
}

if(file.exists(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/distances","_thrBranch",thr_branch_break,".rData"))){
  distances <- loadRData(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/distances","_thrBranch",thr_branch_break,".rData"))
}else{
  distances=computeDistances(xN,num_of_levels,df,Clus,use_only_leafs,sim,altsim,clep,num_of_topics)
  save(distances,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/distances","_thrBranch",thr_branch_break,".rData"))
  print("End distances!!!!!")
}

num_of_classes_is_set <- T

# save rData to disk
#if (save_rData){
#  out <- paste0(output_folder,slash_for_topics,data_topic,"/RDATA") 
#  if(!file.exists(out)){ 
#    dir.create(out)
#  }
#  out <- paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo) 
#  if(!file.exists(out)){ 
#    dir.create(out)
#  }
#  save(distances,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/distances.rData"))
#}

################### graph
thrt = "finalID"
#thr = 40
print(thr_graph_connections)
thr = thr_graph_connections  #percentage - remove all those connections that have a dissimilarity percentage > thr %
netyp = "Whole"
net_sil = FALSE

print("Create the network!!!!")
nt=Netw(xN,clep,num_of_levels,thr,thrt,netyp,df,Clus,use_only_leafs, distances$ffg, distances$ffg2, distances$ffg3, distances$ffg2sim, net_sil)

if (use_only_leafs) {u="use_only_leafs"} else u="use_all_nodes"

pal1 <- rainbow(6, alpha=1) 
colors=c("red","blue","green","pink","grey") #red:strong relationship

print("Plot the network!!!!")

#whole
png(paste0(output_folder,slash_for_topics,data_topic,"/","whole_graph_thr ",thr," ",u,".png"), height = 2000, width = 2000)
plot(nt$net0.copy, edge.curved=.1, vertex.label.color = "black",edge.width=1)  #plot the network graph
legend("topleft", inset=c(0.1,0.2), colnames(df[str_which(names(df), "level.")])[1:(lev+1)], pch=21,
       col="#777777", pt.bg=unique(V(nt$net0.copy)$color), pt.cex=2, cex=.8, bty="n", ncol=1)
legend("topright", inset=c(0.1,0.2), c("4-5","3-4","2-3","1-3","0-1"), pch=21,
       col="#777777", pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1,title = "Relationship strength")
dev.off()

#partial
png(paste0(output_folder,slash_for_topics,data_topic,"/","partial_graph_thr ",thr," ",u,".png"), height = 2000, width = 2000)
plot(nt$net1.copy, edge.curved=.1, vertex.label.color = "black",edge.width=1)  #plot the network graph
legend("topleft", inset=c(0.1,0.2), paste(unique(c(0,clep[useful_nodes]))), pch=21,
       col="#777777", pt.bg=unique(V(nt$net0.copy)$color), pt.cex=2, cex=.8, bty="n", ncol=1)
legend("topright", inset=c(0.1,0.2), c("4-5","3-4","2-3","1-3","0-1"), pch=21,
       col="#777777", pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1,title = "Relationship strength")
dev.off()


#matrix.heatmap(nt$thrtyp)

thr_distance <- c(5,10,15,20,25,30,35,40,45,50)

# merging without using standard thresholds
exp_thresholds_merged <- data.frame(matrix(0, nrow = 40, ncol = 10)) 
colnames(exp_thresholds_merged) <- c("thr", "num_of_cl", "Identity","Similarity","Entropy_Id", "Entropy_Sim", "BS", "TopicSim", "TS_notNA","FScore")

thrt = "finalID"
netyp = "Whole"
net_sil = FALSE
iter <- 0
thr <- 5
num_of_cl_merged <- 0

using_standard_thr <- F

if (domain != 'AA'){
  print(paste0("actual num of classes = ",length(unique(actual_classes))))
  if (!using_standard_thr){
    while (num_of_cl_merged != length(unique(actual_classes))){
      iter <- iter + 1
      if (iter > 40) break
      nt=Netw(xN,clep,num_of_levels,thr,thrt,netyp,df,Clus,use_only_leafs, distances$ffg, distances$ffg2, distances$ffg3, distances$ffg2sim,net_sil)
      edge_table <- nt$tempor3 #the matrix with only those edges that have high connectivity 
      if (sum(edge_table,na.rm = T)==0){
        thr <- thr + 5
        next
      }  
      
      # for each row find the nearest neighbor ie the colum with the smallest percentage!!!
      edge_table2 <- edge_table
      
      for (i in 1:nrow(edge_table)){
        if (length(which(edge_table[i,]>0))){
          min_value <- min(edge_table[i,][which(edge_table[i,]>0)], na.rm = T)
          edge_table2[i, which(edge_table[i,] != min_value)] <- 0
        }
        #if you have equal distances then... take into account something more important?
      }
      
      edges <- as.data.frame(which(edge_table2>0,arr.ind=T))
      
      # find all those nodes that are only in pairs!!!
      edges2 <- edges
      edges2$distance <- NA
      edges2$seqnum_row <- NA
      edges2$seqnum_col <- NA
      
      #number of elements per node
      n <- Clus$seqnum[useful_nodes+1]
      
      n_clusters <- data.frame(node_id=useful_nodes,n=n)
      n_clusters <- n_clusters[order(-n),]
      
      # edges2: row=the bigger, column=the smallest cluster of the two
      for (i in 1:nrow(edges)){
        if (n[edges[i,1]]<n[edges[i,2]]){
          edges2[i,1]<-useful_nodes[edges[i,1]]
          edges2[i,2]<-useful_nodes[edges[i,2]]
        }else{
          edges2[i,1]<-useful_nodes[edges[i,2]]
          edges2[i,2]<-useful_nodes[edges[i,1]]
        }
        edges2$distance[i] <- edge_table[min(edges[i,1],edges[i,2]), max(edges[i,1],edges[i,2])]
        edges2$seqnum_row[i] <- Clus$seqnum[(edges2[i,1]+1)]
        edges2$seqnum_col[i] <- Clus$seqnum[(edges2[i,2]+1)]
      }
      
      edges2 <- edges2[order(edges2$distance*edges2$seqnum_row),]
      
      if (num_of_classes_is_set){
        set <- 1
      }else{
        set <- nrow(edges2)
      }
      
      for (set_n in set:nrow(edges2)){
        ################### merge ################### 
        merge_list <- list()
        for (i in 1:set_n){
          if (paste0(edges2[i,2]) %in% names(merge_list)){
            merge_list[[paste0(edges2[i,2])]] <- c(merge_list[[paste0(edges2[i,2])]], edges2[i,1]) 
          }else{
            merge_list[[paste0(edges2[i,2])]] <- edges2[i,1]
          }
        }
        
        #if there are nodes with more than 1 connections
        merge_list2 <- merge_list
        for (i in 1:length(merge_list)){
          if (i<length(merge_list)){
            for (j in (i+1):length(merge_list)){
              merge_list2[[j]][which(merge_list2[[j]] %in% merge_list2[[i]])] <- as.numeric(names(merge_list2)[i])
            }
          }
        }
        
        merge_list3 <- merge_list2
        for (i in 1:length(merge_list)){
          index <- which(merge_list2[[i]] %in% as.numeric(names(merge_list2)))
          if (length(index)>0){
            for (j in 1:length(index)){
              merge_list3[[i]] <- unique(c(merge_list3[[i]],merge_list3[[paste0(merge_list3[[i]][index[j]])]]))
            }
          }
        }
        
        #delete elements from the list  
        to_be_deleted <- c()
        for (i in 1:length(merge_list)){
          if (i<length(merge_list)){
            for (j in (i+1):length(merge_list)){
              if (as.numeric(names(merge_list2)[i]) %in% merge_list2[[j]]){
                to_be_deleted <- c(to_be_deleted,names(merge_list2)[i])
                break()
              }
            }
          }
        }
        
        merge_list_final <- list()
        for (i in 1:length(merge_list)){
          if (!(names(merge_list3)[i] %in% to_be_deleted)){
            merge_list_final[[names(merge_list3)[i]]] <- merge_list3[[i]]
          }
        }
        
        #Merged matrices 
        df$clusters <- apply(df[,4:(ncol(df)-1)],1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
        df_merged <- df[,1:3]
        for (i in 1:length(merge_list_final)){
          df_merged[which(df_merged$clusters %in% merge_list_final[[i]]),3] <- as.numeric(names(merge_list_final)[i])
        }
        
        
        useful_nodes_merged <- unique(df_merged$clusters)
        num_of_cl_merged <- length(useful_nodes_merged)
        
        if ((num_of_classes_is_set) & (num_of_cl_merged == length(unique(actual_classes)))) break 
      }
      
      print(paste0("thr = ", thr, " --- num of merged clusters = ", num_of_cl_merged))
      
      exp_thresholds_merged$thr[iter] <- thr
      exp_thresholds_merged$num_of_cl[iter] <- num_of_cl_merged
      
      if (num_of_cl_merged != length(unique(actual_classes))){
        if ((num_of_cl_merged - length(unique(actual_classes))) > 5){
          thr <- thr +5
          next
        }else if ((num_of_cl_merged - length(unique(actual_classes))) > 1){
          thr <- thr + 1 
        }else if ((num_of_cl_merged - length(unique(actual_classes))) == 1 ){
          thr <- thr + 0.15
        }else if ((length(unique(actual_classes)) - num_of_cl_merged) > 3){
          thr <- thr - 1
        }else if ((length(unique(actual_classes)) - num_of_cl_merged) >= 1){
          thr <- thr - 0.15
        }
        
        if (iter > 1){
          if (length(which(exp_thresholds_merged$thr[iter] == exp_thresholds_merged$thr[1:(iter-1)]))>2){
            break
          }
          if (exp_thresholds_merged$num_of_cl[iter] %in% exp_thresholds_merged$num_of_cl[1:(iter-1)]){
            next
          }
        }
      }
      
      
      ################### metrics ################### 
      merged_mean_metrics=data.frame(matrix(0, nrow = 4, ncol = 15))
      colnames(merged_mean_metrics)=c("ClusterId", "Identity","Similarity","Entropy_Id", "Entropy_Sim", "BS", "seqnum", "level", "Topic1", 
                                      "Topic1_notA","Topic2","Topic2_notA", "Topic3", "Topic3_notA", "TopicSim" )
      
      Clus_merged = as.data.frame(matrix(100, ncol = 6, nrow = num_of_cl_merged))
      names(Clus_merged) = c("ClusterId","Identity","Similarity","Entropy_Id","Entropy_Sim","BS")
      Clus_merged$ClusterId = useful_nodes_merged
      Clus_merged$seqnum[1] = NA
      Clus_merged$level = NA
      Clus_merged$Entropy_Id=0
      Clus_merged$Entropy_Sim=0
      
      Clus_merged$Topic1=NA
      Clus_merged$Topic1_notA=NA
      Clus_merged$Topic2=NA
      Clus_merged$Topic2_notA=NA
      Clus_merged$Topic3=NA
      Clus_merged$Topic3_notA=NA
      Clus_merged$TopicSim=NA
      
      
      #Plot 
      # center=c()
      # d=c()
      for (cl in 1:length(useful_nodes_merged)){
        a <- df_merged %>% filter(clusters==useful_nodes_merged[cl])
        
        #metrics for evaluation 
        matrix_merged=create_freq_matrix_for_merged_clusters(a,groups,let,sim,max_group_length)
        Clus_merged$seqnum[cl] = nrow(a)
        Clus_merged$Identity[cl]=100*length(which(matrix_merged$permat==100))/str_length(udata$AA.JUNCTION[1])
        Clus_merged$Similarity[cl]=100*length(which(matrix_merged$persim==100))/str_length(udata$AA.JUNCTION[1])
        Clus_merged$Entropy_Id[cl]=mean(matrix_merged$permat[nrow(matrix_merged$permat),])
        Clus_merged$Entropy_Sim[cl]=mean(matrix_merged$persim[nrow(matrix_merged$persim),])
        
        Clus_merged$BS[cl]=mean(colSums(computeSimilarityCl(matrix_merged$permat,letter_sim,let)[1:(nrow(matrix_merged$permat)-1),]))
        
        Clus_merged$Topic1[cl]=as.numeric(names(findTopicsPerCl(matrix_merged$udata,3,let))[1])
        Clus_merged$Topic2[cl]=as.numeric(names(findTopicsPerCl(matrix_merged$udata,3,let))[2])
        Clus_merged$Topic3[cl]=as.numeric(names(findTopicsPerCl(matrix_merged$udata,3,let))[3])
        Clus_merged$Topic1_notA[cl]=as.numeric(findTopicsPerCl(matrix_merged$udata,3,let)[1])
        Clus_merged$Topic2_notA[cl]=as.numeric(findTopicsPerCl(matrix_merged$udata,3,let)[2])
        Clus_merged$Topic3_notA[cl]=as.numeric(findTopicsPerCl(matrix_merged$udata,3,let)[3])
        if (!is.na(Clus_merged$Topic2[cl])){
          Clus_merged$TopicSim[cl]=topic_similarity[Clus_merged$Topic1[cl],Clus_merged$Topic2[cl]]
        }else{
          Clus_merged$TopicSim[cl]=1
        }
        
      }
      
      temp=Clus_merged$Topic1_notA
      Clus_merged$TopicSim[which(temp<0.3)]=NA
      Clus_merged$Topic2_notA[which(temp<0.3)]=NA
      
      temp=Clus_merged$Topic2_notA
      #Clus_merged$Topic2_notA[which(temp<0.3)]=NA
      #Clus_merged$Topic3_notA[which(temp<0.3)]=NA
      #Clus_merged$Topic2[which(temp<0.3)]=NA
      #Clus_merged$TopicSim[which(Clus_merged$seqnum<10)]=NA
      #Clus_merged$Topic3[which(temp<0.3)]=NA
      Clus_merged$TopicSim[which(temp<0.3)]=1
      Clus_merged$TopicSim[which(Clus_merged$seqnum<5)]=NA
      
      #mean values
      mean(Clus_merged$Identity)
      mean(Clus_merged$Similarity)
      mean(Clus_merged$Entropy_Id)
      mean(Clus_merged$Entropy_Sim)
      mean(Clus_merged$BS)
      mean(na.omit(Clus_merged$TopicSim))
      
      #Bar plot
      g_merged <- df_merged[order(df_merged$Sequence.ID),]$clusters
      #Bar_plots_general(g_merged, paste0("_merged_thr_",thr,"_thrBranch",thr_branch_break),output_folder,"","",algo,num_of_topics)
      
      ################ FScore ################
      Fscore <- 0
      F2 <- c()
      pres_all <- c()
      for (r in 1:length(labels)){
        # find the documents that belong to class r and save them into dataframe content 
        nr <- length(which(actual_classes==labels[r]))
        FLr <- c()
        pres_all_temp <- c()
        recall_all_temp <- c()
        for (i in unique(g_merged)){
          #find the documents that belong to cluster i
          content <- which(g_merged==i)
          ni <- length(content)
          nri <- length(which(actual_classes[content]==labels[r]))
          precision <- nri/ni
          recall <- nri/nr
          if ((precision+recall)>0){
            F1 <- 2*precision*recall/(precision+recall)
          }else{
            F1 <- 0
          }
          FLr <- c(FLr, F1)
          pres_all_temp <- c(pres_all_temp,precision)
          recall_all_temp <- c(recall_all_temp,recall)
        }
        
        #for each class find the maximum Fscore from the tree levels 
        F2 <- c(F2, max(FLr, na.rm = T))
        pres_all <- c(pres_all, max(pres_all_temp))
        Fscore <- Fscore + max(FLr, na.rm = T) * nr/nrow(udata)
      }
      
      #or find the index of the node (merged_groups) that each node of the first column is related to
      #merged_groups<- unique(edges2$col)
      #merged_groups[edges2 %>% group_indices(col)]
      
      ################## save metrics ###############
      exp_thresholds_merged$Identity[iter] <- mean(Clus_merged$Identity)
      exp_thresholds_merged$Similarity[iter] <- mean(Clus_merged$Similarity)
      exp_thresholds_merged$Entropy_Id[iter] <- mean(Clus_merged$Entropy_Id)
      exp_thresholds_merged$Entropy_Sim[iter] <- mean(Clus_merged$Entropy_Sim)
      exp_thresholds_merged$BS[iter] <- mean(Clus_merged$BS)
      exp_thresholds_merged$TopicSim[iter] <- mean(na.omit(Clus_merged$TopicSim))
      exp_thresholds_merged$Entropy_Id[iter] <-mean(Clus_merged$Entropy_Id)
      exp_thresholds_merged$FScore[iter] <- Fscore
      exp_thresholds_merged$TS_notNA[iter] <- length(which(!is.na(Clus_merged$TopicSim)))
      
      write.table(Clus_merged,paste0(output_folder,slash_for_topics,data_topic,"/Clus_merged_thr_",thr,"_thrBranch",thr_branch_break,".txt"),sep="\t",row.names = F)
      
    }
    # Write table
    #exp_thresholds[,2:ncol(exp_thresholds)]=apply(exp_thresholds[,2:ncol(exp_thresholds)], 2, function(x) formatC(x, format = "f", digits = 3))
    write.table(exp_thresholds_merged,paste0(output_folder,slash_for_topics,data_topic,"/exp_thresholds_merged_thrBranch",thr_branch_break,".txt"),sep="\t",row.names = F)
    
  }
  
  
  # merging using standard thresholds
  thr_distance <- c(5,10,15,20,25,30,35,40,45)
  exp_thresholds_merged <- data.frame(matrix(0, nrow = 40, ncol = 10)) 
  colnames(exp_thresholds_merged) <- c("thr", "num_of_cl", "Identity","Similarity","Entropy_Id", "Entropy_Sim", "BS", "TopicSim", "TS_notNA","FScore")
  
  thrt = "finalID"
  netyp = "Whole"
  net_sil = FALSE
  iter <- 0
  thr <- 5
  num_of_cl_merged <- 0
  
  if (using_standard_thr){
    print(paste0("actual num of classes = ",length(unique(actual_classes))))
    for (thr in thr_distance ){
      iter <- iter + 1
      if (iter > 40) break
      nt=Netw(xN,clep,num_of_levels,thr,thrt,netyp,df,Clus,use_only_leafs, distances$ffg, distances$ffg2, distances$ffg3, distances$ffg2sim,net_sil)
      edge_table <- nt$tempor3 #the matrix with only those edges that have high connectivity 
      if (sum(edge_table,na.rm = T)==0){
        next
      }  
      
      # for each row find the nearest neighbor ie the colum with the smallest percentage!!!
      edge_table2 <- edge_table
      
      #for (i in 1:nrow(edge_table)){
      #if (length(which(edge_table[i,]>0))){
      #min_value <- min(edge_table[i,][which(edge_table[i,]>0)], na.rm = T)
      #edge_table2[i, which(edge_table[i,] != min_value)] <- 0
      #}
      #if you have equal distances then... take into account something more important?
      #}
      
      edges <- as.data.frame(which(edge_table2>0,arr.ind=T))
      
      # find all those nodes that are only in pairs!!!
      edges2 <- edges
      edges2$distance <- NA
      edges2$seqnum_row <- NA
      edges2$seqnum_col <- NA
      edges2$non_unique_row <- NA
      
      #number of elements per node
      n <- Clus$seqnum[useful_nodes+1]
      
      n_clusters <- data.frame(node_id=useful_nodes,n=n)
      n_clusters <- n_clusters[order(-n),]
      
      # edges2: row=the bigger, column=the smallest cluster of the two
      for (i in 1:nrow(edges)){
        if (n[edges[i,1]]<n[edges[i,2]]){
          edges2[i,1]<-useful_nodes[edges[i,1]]
          edges2[i,2]<-useful_nodes[edges[i,2]]
        }else{
          edges2[i,1]<-useful_nodes[edges[i,2]]
          edges2[i,2]<-useful_nodes[edges[i,1]]
        }
        edges2$distance[i] <- edge_table[min(edges[i,1],edges[i,2]), max(edges[i,1],edges[i,2])]
        edges2$seqnum_row[i] <- Clus$seqnum[(edges2[i,1]+1)]
        edges2$seqnum_col[i] <- Clus$seqnum[(edges2[i,2]+1)]
      }
      
      edges2 <- edges2[order(edges2$distance*edges2$seqnum_row),]
      
      for (i in unique(edges2$row)){
        temp_edge <- edges2 %>% filter(row==i)
        if (nrow(temp_edge)>1){
          edges2$non_unique_row[which(edges2$row==temp_edge$row[1])[2:nrow(temp_edge)]] <- 1
        }
      }
      
      edges2 <- rbind(edges2 %>% filter(is.na(edges2$non_unique_row)), edges2 %>% filter(!is.na(edges2$non_unique_row)))
      
      if (num_of_classes_is_set){
        set <- 1
      }else{
        set <- nrow(edges2)
      }
      
      for (set_n in set:nrow(edges2)){
        ################### merge ################### 
        merge_list <- list()
        for (i in 1:set_n){
          if (paste0(edges2[i,2]) %in% names(merge_list)){
            merge_list[[paste0(edges2[i,2])]] <- c(merge_list[[paste0(edges2[i,2])]], edges2[i,1]) 
          }else{
            merge_list[[paste0(edges2[i,2])]] <- edges2[i,1]
          }
        }
        
        #if there are nodes with more than 1 connections
        merge_list2 <- merge_list
        for (i in 1:length(merge_list)){
          if (i<length(merge_list)){
            for (j in (i+1):length(merge_list)){
              merge_list2[[j]][which(merge_list2[[j]] %in% merge_list2[[i]])] <- as.numeric(names(merge_list2)[i])
            }
          }
        }
        
        merge_list3 <- merge_list2
        for (i in 1:length(merge_list)){
          index <- which(merge_list2[[i]] %in% as.numeric(names(merge_list2)))
          if (length(index)>0){
            for (j in 1:length(index)){
              merge_list3[[i]] <- unique(c(merge_list3[[i]],merge_list3[[paste0(merge_list3[[i]][index[j]])]]))
            }
          }
        }
        
        #delete elements from the list  
        to_be_deleted <- c()
        for (i in 1:length(merge_list)){
          if (i<length(merge_list)){
            for (j in (i+1):length(merge_list)){
              if (as.numeric(names(merge_list2)[i]) %in% merge_list2[[j]]){
                to_be_deleted <- c(to_be_deleted,names(merge_list2)[i])
                break()
              }
            }
          }
        }
        
        merge_list_final <- list()
        for (i in 1:length(merge_list)){
          if (!(names(merge_list3)[i] %in% to_be_deleted)){
            merge_list_final[[names(merge_list3)[i]]] <- merge_list3[[i]]
          }
        }
        
        #Merged matrices 
        df$clusters <- apply(df[,4:(ncol(df)-1)],1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
        df_merged <- df[,1:3]
        for (i in 1:length(merge_list_final)){
          df_merged[which(df_merged$clusters %in% merge_list_final[[i]]),3] <- as.numeric(names(merge_list_final)[i])
        }
        
        
        useful_nodes_merged <- unique(df_merged$clusters)
        num_of_cl_merged <- length(useful_nodes_merged)
        
        if ((num_of_classes_is_set) & (num_of_cl_merged == length(unique(actual_classes)))) break 
      }
      
      print(paste0("thr = ", thr, " --- num of merged clusters = ", num_of_cl_merged))
      
      exp_thresholds_merged$thr[iter] <- thr
      exp_thresholds_merged$num_of_cl[iter] <- num_of_cl_merged
      
      if (num_of_cl_merged != length(unique(actual_classes))){
        if ((num_of_cl_merged - length(unique(actual_classes))) > 5){
          next 
        }
        if (exp_thresholds_merged$num_of_cl[iter] %in% exp_thresholds_merged$num_of_cl[1:(iter-1)]){
          next 
        }
      }
      
      
      ################### metrics ################### 
      merged_mean_metrics=data.frame(matrix(0, nrow = 4, ncol = 15))
      colnames(merged_mean_metrics)=c("ClusterId", "Identity","Similarity","Entropy_Id", "Entropy_Sim", "BS", "seqnum", "level", "Topic1", 
                                      "Topic1_notA","Topic2","Topic2_notA", "Topic3", "Topic3_notA", "TopicSim" )
      
      Clus_merged = as.data.frame(matrix(100, ncol = 6, nrow = num_of_cl_merged))
      names(Clus_merged) = c("ClusterId","Identity","Similarity","Entropy_Id","Entropy_Sim","BS")
      Clus_merged$ClusterId = useful_nodes_merged
      Clus_merged$seqnum[1] = NA
      Clus_merged$level = NA
      Clus_merged$Entropy_Id=0
      Clus_merged$Entropy_Sim=0
      
      Clus_merged$Topic1=NA
      Clus_merged$Topic1_notA=NA
      Clus_merged$Topic2=NA
      Clus_merged$Topic2_notA=NA
      Clus_merged$Topic3=NA
      Clus_merged$Topic3_notA=NA
      Clus_merged$TopicSim=NA
      
      
      #Plot 
      # center=c()
      # d=c()
      for (cl in 1:length(useful_nodes_merged)){
        a <- df_merged %>% filter(clusters==useful_nodes_merged[cl])
        
        #metrics for evaluation 
        matrix_merged=create_freq_matrix_for_merged_clusters(a,groups,let,sim,max_group_length)
        Clus_merged$seqnum[cl] = nrow(a)
        Clus_merged$Identity[cl]=100*length(which(matrix_merged$permat==100))/str_length(udata$AA.JUNCTION[1])
        Clus_merged$Similarity[cl]=100*length(which(matrix_merged$persim==100))/str_length(udata$AA.JUNCTION[1])
        Clus_merged$Entropy_Id[cl]=mean(matrix_merged$permat[nrow(matrix_merged$permat),])
        Clus_merged$Entropy_Sim[cl]=mean(matrix_merged$persim[nrow(matrix_merged$persim),])
        
        Clus_merged$BS[cl]=mean(colSums(computeSimilarityCl(matrix_merged$permat,letter_sim,let)[1:(nrow(matrix_merged$permat)-1),]))
        
        Clus_merged$Topic1[cl]=as.numeric(names(findTopicsPerCl(matrix_merged$udata,3,let))[1])
        Clus_merged$Topic2[cl]=as.numeric(names(findTopicsPerCl(matrix_merged$udata,3,let))[2])
        Clus_merged$Topic3[cl]=as.numeric(names(findTopicsPerCl(matrix_merged$udata,3,let))[3])
        Clus_merged$Topic1_notA[cl]=as.numeric(findTopicsPerCl(matrix_merged$udata,3,let)[1])
        Clus_merged$Topic2_notA[cl]=as.numeric(findTopicsPerCl(matrix_merged$udata,3,let)[2])
        Clus_merged$Topic3_notA[cl]=as.numeric(findTopicsPerCl(matrix_merged$udata,3,let)[3])
        if (!is.na(Clus_merged$Topic2[cl])){
          Clus_merged$TopicSim[cl]=topic_similarity[Clus_merged$Topic1[cl],Clus_merged$Topic2[cl]]
        }else{
          Clus_merged$TopicSim[cl]=1
        }
        
      }
      
      temp=Clus_merged$Topic1_notA
      Clus_merged$TopicSim[which(temp<0.3)]=NA
      Clus_merged$Topic2_notA[which(temp<0.3)]=NA
      
      temp=Clus_merged$Topic2_notA
      #Clus_merged$Topic2_notA[which(temp<0.3)]=NA
      #Clus_merged$Topic3_notA[which(temp<0.3)]=NA
      #Clus_merged$Topic2[which(temp<0.3)]=NA
      #Clus_merged$TopicSim[which(Clus_merged$seqnum<10)]=NA
      #Clus_merged$Topic3[which(temp<0.3)]=NA
      Clus_merged$TopicSim[which(temp<0.3)]=1
      Clus_merged$TopicSim[which(Clus_merged$seqnum<5)]=NA
      
      #mean values
      mean(Clus_merged$Identity)
      mean(Clus_merged$Similarity)
      mean(Clus_merged$Entropy_Id)
      mean(Clus_merged$Entropy_Sim)
      mean(Clus_merged$BS)
      mean(na.omit(Clus_merged$TopicSim))
      
      #Bar plot
      g_merged <- df_merged[order(df_merged$Sequence.ID),]$clusters
      #Bar_plots_general(g_merged, paste0("_merged_thr_",thr,"_thrBranch",thr_branch_break),output_folder,"","",algo,num_of_topics)
      
      ################ FScore ################
      Fscore <- 0
      F2 <- c()
      pres_all <- c()
      for (r in 1:length(labels)){
        # find the documents that belong to class r and save them into dataframe content 
        nr <- length(which(actual_classes==labels[r]))
        FLr <- c()
        pres_all_temp <- c()
        recall_all_temp <- c()
        for (i in unique(g_merged)){
          #find the documents that belong to cluster i
          content <- which(g_merged==i)
          ni <- length(content)
          nri <- length(which(actual_classes[content]==labels[r]))
          precision <- nri/ni
          recall <- nri/nr
          if ((precision+recall)>0){
            F1 <- 2*precision*recall/(precision+recall)
          }else{
            F1 <- 0
          }
          FLr <- c(FLr, F1)
          pres_all_temp <- c(pres_all_temp,precision)
          recall_all_temp <- c(recall_all_temp,recall)
        }
        
        #for each class find the maximum Fscore from the tree levels 
        F2 <- c(F2, max(FLr, na.rm = T))
        pres_all <- c(pres_all, max(pres_all_temp))
        Fscore <- Fscore + max(FLr, na.rm = T) * nr/nrow(udata)
      }
      
      #or find the index of the node (merged_groups) that each node of the first column is related to
      #merged_groups<- unique(edges2$col)
      #merged_groups[edges2 %>% group_indices(col)]
      
      ################## save metrics ###############
      exp_thresholds_merged$Identity[iter] <- mean(Clus_merged$Identity)
      exp_thresholds_merged$Similarity[iter] <- mean(Clus_merged$Similarity)
      exp_thresholds_merged$Entropy_Id[iter] <- mean(Clus_merged$Entropy_Id)
      exp_thresholds_merged$Entropy_Sim[iter] <- mean(Clus_merged$Entropy_Sim)
      exp_thresholds_merged$BS[iter] <- mean(Clus_merged$BS)
      exp_thresholds_merged$TopicSim[iter] <- mean(na.omit(Clus_merged$TopicSim))
      exp_thresholds_merged$Entropy_Id[iter] <-mean(Clus_merged$Entropy_Id)
      exp_thresholds_merged$FScore[iter] <- Fscore
      exp_thresholds_merged$TS_notNA[iter] <- length(which(!is.na(Clus_merged$TopicSim)))
      
      write.table(Clus_merged,paste0(output_folder,slash_for_topics,data_topic,"/Clus_merged_thr_",thr,"_thrBranch",thr_branch_break,".txt"),sep="\t",row.names = F)
      
    }
    
    # Write table
    #exp_thresholds[,2:ncol(exp_thresholds)]=apply(exp_thresholds[,2:ncol(exp_thresholds)], 2, function(x) formatC(x, format = "f", digits = 3))
    write.table(exp_thresholds_merged,paste0(output_folder,slash_for_topics,data_topic,"/exp_thresholds_merged_with_standard_thr_thrBranch_",thr_branch_break,".txt"),sep="\t",row.names = F)
    
  }
}

######################### Silhouette #########################
#if(net_sil == TRUE){
#  use_only_leafs=F
#  distances=computeDistances(xN,num_of_levels,df,lastlist,Clus,use_only_leafs)
#  Clus$Silhouette=NA
#  newffg = normalize(distances$ffg)
#  newffg2 = normalize(distances$ffg2)
#  newffg2sim = normalize(distances$ffg2sim)
#  newffg3 = newffg * 0.5 + newffg2 * 0.25 + newffg2sim * 0.25  
#  nnnn = normalize(newffg3)
  
  #this indecates how similar the clusters that belong to the same level are. when you go down the hierarchy of the tree the clusters 
  #that belong to each level become more similar, so they have a smaller silhouette value
  
#  if (use_only_leafs){
#    jhj = silhouette(clep[t1],t(nnnn)) 
#    Clus$Silhouette[t1+1]=jhj[,3]
#  }else{
#    jhj = silhouette(Clus$level[which(!is.na(Clus$ClusterId))][2:length(Clus$level[which(!is.na(Clus$ClusterId))])],t(nnnn))
#    Clus$Silhouette[which(!is.na(Clus$ClusterId))[2:length(which(!is.na(Clus$ClusterId)))]]=jhj[,3]
#  }
  
#}
