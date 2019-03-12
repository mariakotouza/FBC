load_rData=F
save_rData=T

if (load_rData){
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/x_initial.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/df_initial.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/Clus_initial.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/clep_initial.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/xN_initial.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/ff_initial.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/perlist.rData"))
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/persimlist.rData"))

  lev = max(na.omit(clep_initial))
  flagtic <- F
}

# real classes
if (file.exists(paste0("../",dataset,slash_for_topics,data_topic,"/true_classes.txt"))){
  actual_classes <- read.csv(paste0("../",dataset,slash_for_topics,data_topic,"/true_classes.txt"), stringsAsFactors = T, header = F)$V1
  if (is.numeric(actual_classes)){
    num_of_actual_cl <- length(unique(actual_classes))
    labels <- unique(actual_classes)
  }else{
    num_of_actual_cl <- length(labels)
    labels <- labels
  }
}


######################### Experiments for thresholds #########################
thr1=10
thrId=70 # change this to have more clusters and more clusters!!!
inner_thr=c(0.25,0.5,1,5)
thr_BS_initial=c(0,0,0,0)

exp_thresholds=data.frame(matrix(0, nrow = length(inner_thr), ncol = 16))
colnames(exp_thresholds)=c("Thr_BS_initial","Thr_Identity","Thr_entropy","thr_BS_initial","num_of_levels", "num_of_cl", "Identity", "Similarity_mean", "Entropy_Sim", "BS", "TopicSim", "Similarity_median","Entropy_Id","Identity_sd","Fscore", "Fscore_all_levels")

#without the condition having at least 1 non zero column and BS. All inner thr=c(0.25,0.5, 1,5)
for (exp in 1:length(inner_thr)){
  #invert when thresholds have been applyed at least once
  clep=clep_initial
  x=x_initial
  df=df_initial
  Clus=Clus_initial
  ff=ff_initial
  
  #Clus$TopicSim[which(Clus$Topic2_notA<0.3)]=NA
  #Clus$TopicSim[which(Clus$seqnum<10)]=NA
  
  threshold1 = inner_thr[exp]
  threshold2 = inner_thr[exp]
  threshold3 = inner_thr[exp]
  threshold4 = inner_thr[exp]
  enthr=T
  levcut=2
  
  if(enthr== TRUE){
    matches <- regmatches(unique(x_initial$pathString), gregexpr("[[:digit:]]+", unique(x_initial$pathString)))
    arxpath = unique(x_initial$pathString)
    for(i in 1:length(matches)){
      ggg = as.numeric(unlist(matches[i])) +1
      if(levcut <= (length(ggg)-1)){
        for(j in levcut:(length(ggg)-1)){
          if(algo == "Identity"){
            if (Clus_initial$Identity[ggg[j]]>thrId){
              #compute metrics !Add Sihlouette, topic similarity, ect
              members_metric = ((Clus_initial$seqnum[ggg[j]] - Clus_initial$seqnum[ggg[j+1]]) / Clus_initial$seqnum[ggg[j]]) * 100
              similarity_metric = (Clus_initial$Identity[ggg[j+1]] - Clus_initial$Identity[ggg[j]]) * 0.5 + (Clus_initial$Similarity[ggg[j+1]] - Clus_initial$Similarity[ggg[j]]) * 0.5
              silhouette_metric = Clus_initial$Silhouette[ggg[j+1]] - Clus_initial$Silhouette[ggg[j]]
              entropy_metric = abs(((Clus_initial$Entropy_Id[ggg[j+1]] - Clus_initial$Entropy_Id[ggg[j]]) * 0.5 + (Clus_initial$Entropy_Sim[ggg[j+1]] - Clus_initial$Entropy_Sim[ggg[j]]) * 0.5)*100)
              topicSim_metric = (Clus_initial$TopicSim[ggg[j+1]] - Clus_initial$TopicSim[ggg[j]])*100
              PerSim_metric = abs((Clus_initial$BS[ggg[j+1]] - Clus_initial$BS[ggg[j]]))
              
              #if the following condition is true -> branch break
              if((members_metric < threshold1 || similarity_metric < threshold2 || entropy_metric < threshold3 || PerSim_metric < threshold4)){
                deik = which(x_initial$pathString == arxpath[i])
                tem = str_locate(arxpath[i],as.character(ggg[j]-1))
                fftemp = str_sub(x_initial$pathString[deik[1]],1,tem[2])
                ff2 = str_sub(x_initial$pathString[deik[1]],tem[1],tem[2])
                ff3 = as.numeric(ff2)
                x$clusters[deik] = ff3
                ff4 = clep_initial[ff3]
                ff5 = str_which(names(x),sprintf("level.%d",ff4))
                ff6 = str_which(names(x), "level.")
                #replace the current cluster and its kids with NA
                if(is.na(ff5[1]) == FALSE){
                  x[deik,(ff5[1]+1):max(ff6)] = NA
                }
                x$pathString[deik] = fftemp
                break()
              }
            }
          }else if(algo == "Similarity"){
            if (Clus_initial$Similarity[ggg[j]]>thrId){
              members_metric = ((Clus_initial$seqnum[ggg[j]] - Clus_initial$seqnum[ggg[j+1]]) / Clus_initial$seqnum[ggg[j]]) * 100
              similarity_metric = Clus_initial$Similarity[ggg[j+1]] - Clus_initial$Similarity[ggg[j]]
              silhouette_metric = Clus_initial$Silhouette[ggg[j+1]] - Clus_initial$Silhouette[ggg[j]]
              entropy_metric = abs(Clus_initial$Entropy_Sim[ggg[j+1]] - Clus_initial$Entropy_Sim[ggg[j]])*100
              topicSim_metric = abs(Clus_initial$TopicSim[ggg[j+1]] - Clus_initial$TopicSim[ggg[j]])*100
              PerSim_metric = abs((Clus_initial$BS[ggg[j+1]] - Clus_initial$BS[ggg[j]]))
              
              #if the following condition is true -> branch break
              if((members_metric < threshold1 || similarity_metric < threshold2 || entropy_metric < threshold3 || PerSim_metric < threshold4)){
                deik = which(x_initial$pathString == arxpath[i])
                tem = str_locate(arxpath[i],as.character(ggg[j]-1))
                fftemp = str_sub(x_initial$pathString[deik[1]],1,tem[2])
                ff2 = str_sub(x_initial$pathString[deik[1]],tem[1],tem[2])
                ff3 = as.numeric(ff2)
                x$clusters[deik] = ff3
                ff4 = clep_initial[ff3]
                ff5 = str_which(names(x),sprintf("level.%d",ff4))
                ff6 = str_which(names(x), "level.")
                if(is.na(ff5[1]) == FALSE){
                  x[deik,(ff5[1]+1):max(ff6)] = NA
                }
                x$pathString[deik] = fftemp
                break()
              }
            }
          }
        }
      }
    }
    xN <- as.Node(x)
    bo = which(clep_initial == lev)
    levtel = clep_initial
    j = 1
    jfjf = vector()
    for(i in 1:max(bo)){
      temp1 = as.numeric(names(FindNode(xN,(sprintf("%d",i)))$children))
      if(length(temp1) == 1){
        if(temp1[1] %% 2 == 1){
          ff7 = temp1[1] + 1
        }else{
          ff7 = temp1[1] - 1
        }
        ff8 = clep_initial[ff7]
        ff9 = str_which(names(x),sprintf("level.%d",ff8))
        deik2 = which(df_initial[ff9] == ff7)
        x[deik2,ff9] = df_initial[deik2,ff9]
        tem = str_locate(df_initial$pathString[deik2[1]],as.character(ff7))
        fftemp = str_sub(df_initial$pathString[deik2[1]],1,tem[2])
        x$pathString[deik2] = fftemp
        jfjf[j] = ff7
        j = j+1
      }
      if(is.null(FindNode(xN,(sprintf("%d",i)))$isLeaf) && is.element(i, jfjf) == FALSE){
        levtel[i] = NA
        Clus[i+1,] = NA
        ff[i+1,] = NA
      }
    }
    xN <- as.Node(x)
  }
  #df=x
  df = x#[colSums(!is.na(x)) > 0]
  Clus[(which(levtel>ncol(x[colSums(!is.na(x)) > 0])-length(which(!str_detect(names(x),"level")))-1)+1),]=NA
  levtel[which(levtel>ncol(x[colSums(!is.na(x)) > 0])-length(which(!str_detect(names(x),"level")))-1)]=NA
  clep = levtel
  
  print(paste0("thr1=",threshold1,"  thr2=",threshold2))
  exp_thresholds$Thr_BS_initial[exp]=thr_BS_initial[exp]
  exp_thresholds$Thr_Identity[exp]=threshold2
  exp_thresholds$Thr_entropy[exp]=threshold3
  exp_thresholds$thr_BS_initial[exp]=threshold4
  num_of_levels=max(na.omit(clep))
  exp_thresholds$num_of_levels[exp]=num_of_levels
  num_of_cl=length(unique(df$clusters))
  exp_thresholds$num_of_cl[exp]=num_of_cl
  
  #Identity
  exp_thresholds$Identity[exp]=idenlev(num_of_levels,clep,Clus,"Identity",flagtic,logFile)
  #save min,max,mean,standard diviation of identity values of clusters per level
  t1 = which(clep == num_of_levels)
  # epipleon
  xm = as.numeric(as.data.frame(xN$leaves))
  orio = min(which(clep == num_of_levels))
  xm2 = sort(xm[xm < orio])
  t1 = sort(append(xm2,t1,after = length(xm2)))
  id_level=na.omit(Clus[t1+1,]$Identity)
  min_id=min(id_level)
  max_id=max(id_level)
  exp_thresholds$Identity_sd[exp]=sd(id_level)
  
  #Similarity
  exp_thresholds$Similarity_mean[exp]=mean(na.omit(Clus[t1+1,]$Similarity))
  exp_thresholds$Similarity_median[exp]=median(na.omit(Clus[t1+1,]$Similarity))
  
  #Entropy
  exp_thresholds$Entropy_Id[exp]=mean(na.omit(Clus[t1+1,]$Entropy_Id))
  exp_thresholds$Entropy_Sim[exp]=mean(na.omit(Clus[t1+1,]$Entropy_Sim))
  
  #Prelist_sim
  exp_thresholds$BS[exp]=mean(na.omit(Clus[t1+1,]$BS))
  
  #Topic similarity
  exp_thresholds$TopicSim[exp]=mean(na.omit(Clus[t1+1,]$TopicSim))
  
  #bar plot
  BarLev(num_of_levels,perlist,Clus,let,flagtic,logFile,paste0("_thr_BS_initial_",thr_BS_initial[exp],"_thr_inner_thr_",inner_thr[exp]))
  # png(paste0(getwd(),"/",output_folder,slash_for_topics,data_topic,"/","tree_algo_",algo,"_thr_BS_initial_",thr_BS_initial[exp],"_thr_inner_thr_",inner_thr[exp],".png"),width=1900, height=2900, res=100)
  # print({
  #   plot(xN)
  # })
  # dev.off()
  
  a <- x[,which(str_detect(colnames(x),"level"))]
  nodes <- c(na.omit(unique(unlist(a))))[which(c(na.omit(unique(unlist(a))))!=0)]
  
  if (file.exists(paste0("../",dataset,slash_for_topics,data_topic,"/true_classes.txt"))){
    Fscore <- 0
    F2 <- c()
    pres_all <- c()
    for (r in 1:length(labels)){
      # find the documents that belong to class r and save them into dataframe content 
      nr <- length(which(actual_classes==labels[r]))
      FLr <- c()
      pres_all_temp <- c()
      recall_all_temp <- c()
      for (i in  c(t1)){
        #find the documents that belong to cluster i
        content <- df[which(a==i, arr.ind = TRUE)[,1],]
        ni <- nrow(content)
        nri <- length(which(actual_classes[content$Sequence.ID]==labels[r]))
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
      Fscore <- Fscore + max(FLr, na.rm = T) * nr/nrow(a)
    }
    
    exp_thresholds$Fscore[exp] <- Fscore
    
    a <- x[,which(str_detect(colnames(x),"level"))]
    nodes <- c(na.omit(unique(unlist(a))))[which(c(na.omit(unique(unlist(a))))!=0)]
    
    Fscore <- 0
    F2 <- c()
    pres_all <- c()
    for (r in 1:length(labels)){
      # find the documents that belong to class r and save them into dataframe content 
      nr <- length(which(actual_classes==labels[r]))
      FLr <- c()
      pres_all_temp <- c()
      recall_all_temp <- c()
      for (i in  c(na.omit(unique(unlist(a))))){
        #find the documents that belong to cluster i
        content <- df[which(a==i, arr.ind = TRUE)[,1],]
        ni <- nrow(content)
        nri <- length(which(actual_classes[content$Sequence.ID]==labels[r]))
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
      Fscore <- Fscore + max(FLr, na.rm = T) * nr/nrow(a)
    }
    
    exp_thresholds$Fscore_all_levels[exp] <- Fscore
  }
  
  
}

# Write table
exp_thresholds[,2:ncol(exp_thresholds)]=apply(exp_thresholds[,2:ncol(exp_thresholds)], 2, function(x) formatC(x, format = "f", digits = 3))
write.table(exp_thresholds,paste0("../",output_folder,slash_for_topics,data_topic,"/exp_thresholds2_algo_",algo,".txt"),sep="\t",row.names = F)

#################################################################################

# save rData to disk
if (save_rData){
  out <- paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA") 
  if(!file.exists(out)){ 
    dir.create(out)
  }
  out <- paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo) 
  if(!file.exists(out)){ 
    dir.create(out)
  }
  save(x,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/x.rData"))
  save(df,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/df.rData"))
  save(Clus,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/Clus.rData"))
  save(clep,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/clep.rData"))
  save(xN,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/xN.rData"))
}


######## final clusters after applying thresholds ########
t1 = which(clep == num_of_levels)
# epipleon
xm = as.numeric(as.data.frame(xN$leaves))
orio = min(which(clep == num_of_levels))
xm2 = sort(xm[xm < orio])
t1 = sort(append(xm2,t1,after = length(xm2)))

#################### Compute num of clusters/levels ########
num_of_cl=length(t1) 
num_of_levels=ncol(x[colSums(!is.na(x)) > 0])-length(which(!str_detect(names(x),"level")))-1

############## Clus_only_final_level_after_thr
Clus_only_final_level_after_thr=Clus %>% filter(ClusterId %in% t1)
temp=Clus_only_final_level_after_thr$Topic2_notA

temp=Clus_only_final_level_after_thr$Topic1_notA
Clus_only_final_level_after_thr$TopicSim[which(temp<0.3)]=NA


#mean values
mean(Clus_only_final_level_after_thr$Identity)
mean(Clus_only_final_level_after_thr$Similarity)
mean(Clus_only_final_level_after_thr$Entropy_Id)
mean(Clus_only_final_level_after_thr$Entropy_Sim)
mean(Clus_only_final_level_after_thr$BS)
mean(na.omit(Clus_only_final_level_after_thr$TopicSim))

write.table(Clus_only_final_level_after_thr,paste0("../",output_folder,slash_for_topics,data_topic,"/Clus_only_final_level_after_thr_algo_",algo,".txt"),sep="\t",row.names = F)

#clusters with more than 2 elements
Clus_only_final_level_after_thr_more_than2Elem=Clus_only_final_level_after_thr %>% filter(seqnum>2)

#table with mean and sd values per level
metrics_per_level_after_thr=EmPin(Clus,dfsd,flagtic,logFile)
metrics_per_level_after_thr[,3:ncol(metrics_per_level_after_thr)]=apply(metrics_per_level_after_thr[,3:ncol(metrics_per_level_after_thr)], 2, function(x) formatC(x, format = "f", digits = 3))
write.table(metrics_per_level_after_thr, paste0("../",output_folder,slash_for_topics,data_topic,"/metrics_per_level_after_thr_whole_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)

