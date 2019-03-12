#load_rData=F
#save_rData=F

if (load_rData){
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/x.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/df.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/Clus.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/clep.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/xN.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/listb.rData"))
  num_of_levels=ncol(x[colSums(!is.na(x)) > 0])-length(which(!str_detect(names(x),"level")))-1
  flagtic = F
  
  lastlist = listb
  # A list with the permat matrix for every branch
  perlist = lastlist$list
  persimlist = lastlist$listn
}

######################### Distances #########################
use_only_leafs=T
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

distances=computeDistances(xN,num_of_levels,df,Clus,use_only_leafs,sim,altsim)

################### graph
thrt = "finalID"
thr = 10  #percentage - remove all those connections that have a dissimilarity percentage > thr %
netyp = "Whole"
net_sil = FALSE
nt=Netw(xN,num_of_levels,thr,thrt,netyp,df,Clus,use_only_leafs, distances$ffg, distances$ffg2, distances$ffg3, distances$ffg2sim)

if (use_only_leafs) {u="use_only_leafs"} else u="use_all_nodes"

pal1 <- rainbow(6, alpha=1) 
colors=c("red","blue","green","pink","grey") #red:strong relationship

#whole
png(paste0("../",output_folder,slash_for_topics,data_topic,"/","whole_graph_thr ",thr," ",u,".png"), height = 2000, width = 2000)
plot(nt$net0.copy, edge.curved=.1, vertex.label.color = "black",edge.width=1,vertex.label.cex = 4)  #plot the network graph
legend("topleft", inset=c(0.1,0.2), colnames(df[str_which(names(df), "level.")])[1:(lev+1)], pch=21,
       col="#777777", pt.bg=unique(V(nt$net0.copy)$color), pt.cex=2, cex=.8, bty="n", ncol=1)
legend("topright", inset=c(0.1,0.2), c("4-5","3-4","2-3","1-3","0-1"), pch=21,
       col="#777777", pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1,title = "Relationship strength")
dev.off()

#partial
png(paste0("../",output_folder,slash_for_topics,data_topic,"/","partial_graph_thr ",thr," ",u,".png"), height = 2000, width = 2000)
plot(nt$net1.copy, edge.curved=.1, vertex.label.color = "black",edge.width=1, vertex.label.cex = 4)  #plot the network graph
legend("topleft", inset=c(0.1,0.2), paste(unique(c(0,clep[useful_nodes]))), pch=21,
       col="#777777", pt.bg=unique(V(nt$net0.copy)$color), pt.cex=2, cex=.8, bty="n", ncol=1)
legend("topright", inset=c(0.1,0.2), c("4-5","3-4","2-3","1-3","0-1"), pch=21,
       col="#777777", pt.bg=colors, pt.cex=2, cex=.8, bty="n", ncol=1,title = "Relationship strength")
dev.off()


matrix.heatmap(nt$thrtyp)

edge_table <- nt$tempor3 #the matrix with only those edges that have high connectivity 

# for each row find the nearest neighbor ie the colum with the biggest percentage!!!
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

#number of elements per node
n <- Clus$seqnum[useful_nodes+1]

for (i in 1:nrow(edges)){
  if (n[edges[i,1]]<n[edges[i,2]]){
    edges2[i,1]<-useful_nodes[edges[i,1]]
    edges2[i,2]<-useful_nodes[edges[i,2]]
  }else{
    edges2[i,1]<-useful_nodes[edges[i,2]]
    edges2[i,2]<-useful_nodes[edges[i,1]]
  }
}

################### merge
merge_list <- list()
for (i in 1:nrow(edges2)){
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
df_merged <- df[,1:3]
for (i in 1:length(merge_list_final)){
  df_merged[which(df_merged$clusters %in% merge_list_final[[i]]),3] <- as.numeric(names(merge_list_final)[i])
}

useful_nodes_merged <- unique(df_merged$clusters)
num_of_cl_merged <- length(useful_nodes_merged)

################### metrics ################### 
merged_mean_metrics=data.frame(matrix(0, nrow = 4, ncol = 9))
colnames(merged_mean_metrics)=c("num_of_cl", "Identity", "Identity_sd","Similarity_mean", "Similarity_median","Entropy_Id", "Entropy_Sim", "BS", "TopicSim" )

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
  matrix_merged=create_freq_matrix_for_merged_clusters(a,groups,let,sim)
  Clus_merged$seqnum[cl] = nrow(a)
  Clus_merged$Identity[cl]=100*length(which(matrix_merged$permat==100))/str_length(udata$AA.JUNCTION[1])
  Clus_merged$Similarity[cl]=100*length(which(matrix_merged$persim==100))/str_length(udata$AA.JUNCTION[1])
  Clus_merged$Entropy_Id[cl]=mean(matrix_merged$permat[nrow(matrix_merged$permat),])
  Clus_merged$Entropy_Sim[cl]=mean(matrix_merged$persim[nrow(matrix_merged$persim),])
  
  Clus_merged$BS[cl]=mean(colSums(computeSimilarityCl(matrix_merged$permat)[1:(nrow(matrix_merged$permat)-1),]))
  
  Clus_merged$Topic1[cl]=as.numeric(names(findTopicsPerCl(matrix_merged$udata,3))[1])
  Clus_merged$Topic2[cl]=as.numeric(names(findTopicsPerCl(matrix_merged$udata,3))[2])
  Clus_merged$Topic3[cl]=as.numeric(names(findTopicsPerCl(matrix_merged$udata,3))[3])
  Clus_merged$Topic1_notA[cl]=as.numeric(findTopicsPerCl(matrix_merged$udata,3)[1])
  Clus_merged$Topic2_notA[cl]=as.numeric(findTopicsPerCl(matrix_merged$udata,3)[2])
  Clus_merged$Topic3_notA[cl]=as.numeric(findTopicsPerCl(matrix_merged$udata,3)[3])
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
Clus_merged$TopicSim[which(temp<0.3)]=1
Clus_merged$TopicSim[which(Clus_merged$seqnum<5)]=NA

#mean values
mean(Clus_merged$Identity)
mean(Clus_merged$Similarity)
mean(Clus_merged$Entropy_Id)
mean(Clus_merged$Entropy_Sim)
mean(Clus_merged$BS)
mean(Clus_merged$TopicSim)

#Bar plot
Bar_plots_merged(df_merged,useful_nodes_merged, "_merged")

#or find the index of the node (merged_groups) that each node of the first column is related to
#merged_groups<- unique(edges2$col)
#merged_groups[edges2 %>% group_indices(col)]


######################### Silhouette #########################
if(net_sil == TRUE){
  use_only_leafs=F
  distances=computeDistances(xN,num_of_levels,df,lastlist,Clus,use_only_leafs)
  Clus$Silhouette=NA
  newffg = normalize(distances$ffg)
  newffg2 = normalize(distances$ffg2)
  newffg2sim = normalize(distances$ffg2sim)
  newffg3 = newffg * 0.5 + newffg2 * 0.25 + newffg2sim * 0.25  
  nnnn = normalize(newffg3)
  
  #this indecates how similar the clusters that belong to the same level are. when you go down the hierarchy of the tree the clusters 
  #that belong to each level become more similar, so they have a smaller silhouette value
  
  if (use_only_leafs){
    jhj = silhouette(clep[t1],t(nnnn)) 
    Clus$Silhouette[t1+1]=jhj[,3]
  }else{
    jhj = silhouette(Clus$level[which(!is.na(Clus$ClusterId))][2:length(Clus$level[which(!is.na(Clus$ClusterId))])],t(nnnn))
    Clus$Silhouette[which(!is.na(Clus$ClusterId))[2:length(which(!is.na(Clus$ClusterId)))]]=jhj[,3]
  }
  
}
