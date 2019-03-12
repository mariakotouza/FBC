############################### Initialization ############################### 
#change!!
endper = 100 # Set the percentage, which ends the programm 
#exclude_letters <- c()  #c("group_1")
#absoluteIdentity <- T

# An empty vector wich will store the level for every cluster
clep = vector('numeric')
br = 0  # Initial value of branch
cl = 0  # Initial value of new clusters
met = 0 # Initial value of level capacity counter
ep = 1  # Initial value of level
d = 0
nn = FALSE # Initial value for the condition sumper < endper%
listxx = list() # Initialize a list for saving the permat of all branches
listyy = list()
dfsum = data.frame(sumper = numeric(0),sumper2 = numeric(0),branch = numeric(0), len = numeric(0))

# Evaluation metrics -> depends on the domain
evaluation_metrics <- c("ff", "Num_of_clusters", "Average-Identity", "Average-Similarity", "Average-Entropy_Sim",
                        "Average-BS", "Average-TopicSim","SD-Similarity", "SD-Identity","Average-Entropy_Id", 
                        "SD-Entropy_Id", "SD-Entropy_Sim", "SD-BS", "SD-TopicSim")
as.data.frame(matrix(nrow = 0, ncol = (length(evaluation_metrics)-1), 0))
dfsd = data.frame(ff = character(0),as.data.frame(matrix(nrow = 0, ncol = (length(evaluation_metrics)-1), 0)),stringsAsFactors = FALSE)
colnames(dfsd) <- evaluation_metrics

ggdf = data.frame(branch = numeric(0), len = numeric(0))
last = 0
flagtic = TRUE
udata$AA.JUNCTION = as.character(udata$AA.JUNCTION)
udata$clusters = 0 # Initialiaze the column clusters with 0
udata$level.0 = 0 # Initialize the column of cl.0 with 0
udata$temp= 0 # Creating a temp column with 0
progend = FALSE
listax = list()
listax$temp = 1:length(udata$AA.JUNCTION)
names(listax)[length(listax)] = sprintf('cl.%d', 0) 
listq = list("listax" = listax,"ggdf" = ggdf, "dfsum" = dfsum,"list" = listxx, "listn" = listyy,"udata" = udata,"permat"= NA, "persim" = NA, "br" = br, "cl" = cl, "met" = met, "ep" = ep , "clep" = clep, "nn" = nn, "sumper" = NA,"sumper2" = NA, "selected_cell_id" = NA, "cel" = NA,"endper" = endper, "last" = last,"progend" = progend, "leaf" = FALSE,"leaf2" = FALSE)
algocol = 0
backcol = 0
backcolj6 = 0
listb = listq
levcut = 3
levcut =levcut +1
leaf2 = FALSE
enthr = FALSE

#logfile
out <- paste0("../",output_folder,slash_for_topics,data_topic,"/log")
if(!file.exists(paste0("../",output_folder))){ 
  dir.create(paste0("../",output_folder))
}
if(!file.exists(paste0("../",output_folder,slash_for_topics,data_topic))){ 
  dir.create(paste0(paste0("../",output_folder,slash_for_topics,data_topic)))
}
if(!file.exists(paste0(out))){ 
  dir.create(paste0(out))
}
logFile = paste0("../",output_folder,slash_for_topics,data_topic,"/log/log_file ",trunc(as.numeric(Sys.time())),".txt")
cat(paste0("Function","\t","Num of input rows","\t","Num of input columns","\t","Duration"), file=logFile, append=FALSE, sep = "\n")

if (load_rData){
  load(paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/listb.rData"))
}else{
  #################### Run ###########################
  tic()
  while (progend == FALSE){
    lista = Matrices(listb,let,sim,d,algo,algocol,backcol,backcolj6,logFile)
    progend = lista$progend
    if(progend == TRUE){
      print("ended up here")
      break
    }
    listb = Choice(lista,let,sim,d,algo,algocol,backcol,backcolj6,logFile)
    progend = listb$progend
  }
  en = toc(quiet = TRUE)
  cat(paste0("lastlist","\t",length(udata[udata$clusters == br,]$AA.JUNCTION),"\t",str_length(udata$AA.JUNCTION[1]),"\t",en$toc - en$tic), file=logFile, append=TRUE, sep = "\n")
  
}


# save rData to disk
if (save_rData){
  save(listb,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/listb.rData"))
}

############## Compute Metrics for each cluster #####################
lastlist = listb
lastlist2 = listb
lastlist3 = listb
# The final name of udata data frame 
df = lastlist$udata
# A list with the permat matrix for every branch
perlist = lastlist$list
persimlist = lastlist$listn
# Create a dataframe with all clusters and their identity and similarity percentage 
ff = lastlist$dfsum
ff$level[1]= 0
ff$level[2:(length(lastlist$clep[ff$branch]) +1 )] = lastlist$clep[ff$branch]  #without leaves
clep <- lastlist$clep

Clus = as.data.frame(matrix(100, ncol = 6, nrow = max(df$clusters)+1))
names(Clus) = c("ClusterId","Identity","Similarity","Entropy_Id","Entropy_Sim","BS")
Clus$ClusterId = 0:max(df$clusters)
Clus$seqnum[1] = nrow(df) 
Clus$seqnum[2:(max(df$clusters)+1)] = lastlist$ggdf$len
Clus$level = c(0,lastlist$clep)
Clus$Entropy_Id=0
Clus$Entropy_Sim=0
for(i in 1:length(ff$branch) ){
  ll = which(Clus$ClusterId == ff$branch[i])
  Clus$Identity[ll] = ff$sumper[i]
  Clus$Similarity[ll] = ff$sumper2[i]
  Clus$Entropy_Id[ll] = mean(perlist[[i]][nrow(perlist[[i]]),])
  Clus$Entropy_Sim[ll] = mean(persimlist[[i]][nrow(persimlist[[i]]),])
  Clus$BS[ll]=mean(colSums(computeSimilarityCl(perlist[[i]])[1:(nrow(perlist[[i]])-1),])) 
}

##############  Find the pathString for each node/cluster ############## 
lev = max(na.omit(lastlist$clep))
df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
df_args <- c(df[str_which(names(df), "level.")], sep="/")
if(lev == max(na.omit(lastlist$clep))){
  df$pathString<- do.call(paste, df_args)
  kk = df$pathString
  for(i in 1:length(kk)){
    temp = str_locate(kk[i],"/NA")
    if(is.na(temp[1]) == FALSE){
      temp2 = str_sub(kk[i], 1, temp[1]-1);
      kk[i] = temp2 
    }
  } 
}else{
  df$pathString<- do.call(paste, df_args)
  kk = df$pathString
  gg =as.data.frame(str_locate_all(kk,"/"))
  tem = 1
  for (i in 1:length(kk)){
    kk[i] = str_sub(kk[i],1,gg[,tem][lev+1]-1)
    tem = tem + 2
  }
  for(i in 1:length(kk)){
    temp = str_locate(kk[i],"/NA")
    if(is.na(temp[1]) == FALSE){
      temp2 = str_sub(kk[i], 1, temp[1]-1);
      kk[i] = temp2 
    }
  } 
}
df$pathString = kk
x <- ToDataFrameTree(df, "pathstring")
xN <- as.Node(x)


######################## Topic Similarity per cluster ######################## 
Clus$Topic1=NA
Clus$Topic1_notA=NA
Clus$Topic2=NA
Clus$Topic2_notA=NA
Clus$Topic3=NA
Clus$Topic3_notA=NA
Clus$TopicSim=NA
for (i in 1:nrow(Clus)){
  Clus$Topic1[i]=as.numeric(names(findTopicsPerCl(AminoCl(Clus$ClusterId[i],lastlist$clep,df,flagtic,logFile),3))[1])
  Clus$Topic2[i]=as.numeric(names(findTopicsPerCl(AminoCl(Clus$ClusterId[i],lastlist$clep,df,flagtic,logFile),3))[2])
  Clus$Topic3[i]=as.numeric(names(findTopicsPerCl(AminoCl(Clus$ClusterId[i],lastlist$clep,df,flagtic,logFile),3))[3])
  Clus$Topic1_notA[i]=as.numeric(findTopicsPerCl(AminoCl(Clus$ClusterId[i],lastlist$clep,df,flagtic,logFile),3)[1])
  Clus$Topic2_notA[i]=as.numeric(findTopicsPerCl(AminoCl(Clus$ClusterId[i],lastlist$clep,df,flagtic,logFile),3)[2])
  Clus$Topic3_notA[i]=as.numeric(findTopicsPerCl(AminoCl(Clus$ClusterId[i],lastlist$clep,df,flagtic,logFile),3)[3])
  if (!is.na(Clus$Topic2[i])){
    Clus$TopicSim[i]=topic_similarity[Clus$Topic1[i],Clus$Topic2[i]]
  }else{
    Clus$TopicSim[i]=1
  }
}

temp=Clus$Topic1_notA
Clus$TopicSim[which(temp<0.3)]=NA

temp=Clus$Topic2_notA
#Clus$Topic2_notA[which(temp<0.3)]=NA
#Clus$Topic3_notA[which(temp<0.3)]=NA
#Clus$Topic2[which(temp<0.3)]=NA
#Clus$TopicSim[which(Clus_merged$seqnum<10)]=NA
#Clus$Topic3[which(temp<0.3)]=NA
Clus$TopicSim[which(temp<0.3)]=1


#################### Compute num of clusters/levels ########
num_of_cl=length(unique(df$clusters)) 
num_of_levels=ncol(x[colSums(!is.na(x)) > 0])-length(which(!str_detect(names(x),"level")))-1

######## leaves - clusters ########
t1 = which(lastlist$clep == num_of_levels)
# epipleon
xm = as.numeric(as.data.frame(xN$leaves))
orio = min(which(lastlist$clep == num_of_levels))
xm2 = sort(xm[xm < orio])
t1 = sort(append(xm2,t1,after = length(xm2)))

Clus_only_final_level <- Clus %>% filter(ClusterId %in% t1)

#mean values
mean(Clus_only_final_level$Identity)
mean(Clus_only_final_level$Similarity)
mean(Clus_only_final_level$Entropy_Id)
mean(Clus_only_final_level$Entropy_Sim)
mean(Clus_only_final_level$BS)
mean(na.omit(Clus_only_final_level$TopicSim))

#table with mean and sd values per level
metrics_per_level=EmPin(Clus,dfsd,flagtic,logFile)
#format(metrics_per_level, digits=3)
#metrics_per_level[,2:ncol(metrics_per_level)]=round(metrics_per_level[,2:ncol(metrics_per_level)], 3)
metrics_per_level[,3:ncol(metrics_per_level)]=apply(metrics_per_level[,3:ncol(metrics_per_level)], 2, function(x) formatC(x, format = "f", digits = 3))
write.table(metrics_per_level, paste0("../",output_folder,slash_for_topics,data_topic,"/metrics_per_level_whole_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)

#### Save useful variables and matrices that have come from the initial tree #### 
x_initial=x
df_initial=df
Clus_initial=Clus
clep_initial=lastlist$clep
xN_initial <- xN
ff_initial <- ff

# save rData to disk
if (save_rData){
  out <- paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA")
  if(!file.exists(out)){ 
    dir.create(out)
  }
  save(x_initial,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/x_initial.rData"))
  save(df_initial,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/df_initial.rData"))
  save(Clus_initial,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/Clus_initial.rData"))
  save(clep_initial,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/clep_initial.rData"))
  save(xN_initial,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/xN_initial.rData"))
  save(ff_initial,file=paste0("../",output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/ff_initial.rData"))
}


#write.table(Clus_initial,paste0(output_folder,slash_for_topics,data_topic,"/Clus_initial_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)
#write.table(df_initial,paste0(output_folder,slash_for_topics,data_topic,"/df_initial_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)
#write.table(x_initial,paste0(output_folder,slash_for_topics,data_topic,"/x_initial_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)
#write.table(clep_initial,paste0(output_folder,slash_for_topics,data_topic,"/clep_initial_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)








