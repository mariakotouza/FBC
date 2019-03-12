## We assume that the data for item clustering have been created using the LDA script

###################################################################################################
################################# Load Packages ################################# 
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("parallel","plyr","dplyr","data.table","stringr","tidyr","entropy","ggplot2","ggseqlogo","gridExtra","cluster","seqinr","collapsibleTree","data.tree","DiagrammeR","stringdist","igraph","networkD3","plsgenomics","shinycssloaders","shiny","shinyFiles","shinyjs","shinyBS","DT","plotly","xtable","tictoc")
ipak(packages)

# if you only want to run binary tree construction 
#packages <- c("plyr","dplyr","data.table","stringr","tidyr","entropy","data.tree","stringdist","tictoc", "profvis")

###################################################################################################
source("functions.R")

save_rData=T
load_rData=F

########### Initialize bins ########### 
#algo = "Identity"
algo = "Similarity"
domain="items" #or "AA"
dataset_name <- as.character(args[2])
num_of_topics=20
dataset <- paste0("Data/",dataset_name)
output_folder <- paste0("output/", dataset_name)
data_topic=paste0(num_of_topics,"topics")
slash_for_topics <- "/"
num_of_bins=10
num_of_bins_per_group=2
overlaping=F
color_bar=c("lightyellow","honeydew2","lightgoldenrod2","lightgoldenrod4","lightpink","hotpink2","darkolivegreen1","green","deepskyblue","blue")

# Evaluation metrics -> depends on the domain
evaluation_metrics <- c("ff", "Num_of_clusters", "Average-Identity", "Average-Similarity", "Average-Entropy_Sim",
                        "Average-BS", "Average-TopicSim","SD-Similarity", "SD-Identity","Average-Entropy_Id", 
                        "SD-Entropy_Id", "SD-Entropy_Sim", "SD-BS", "SD-TopicSim")
dfsd = data.frame(ff = character(0),as.data.frame(matrix(nrow = 0, ncol = (length(evaluation_metrics)-1), 0)),stringsAsFactors = FALSE)
colnames(dfsd) <- evaluation_metrics

#group values using ranges 0:1 with step=1/num_of_bins
groups=seq(0,1,1/num_of_bins)
#to do: make it more general ~ do not use letters that are only 24, use combinations
let = LETTERS[1:num_of_bins] # The letters matrix
let=let[length(let):1]
sim=list()
color=matrix(0,nrow=1, ncol=length(let))
j=1
if (!overlaping){
  colors=colorRampPalette(c("blue", "green","yellow","red"))(num_of_bins/num_of_bins_per_group)
  for (i in 1:(num_of_bins/num_of_bins_per_group)){
    sim[[paste0(let[j],let[(j+num_of_bins_per_group-1)])]]=let[j:(j+num_of_bins_per_group-1)]
    color[j:(j+num_of_bins_per_group-1)]=colors[i]
    j=j+num_of_bins_per_group
  }
}else{
  colors=colorRampPalette(c("blue", "green","yellow","red"))(num_of_bins-1)
  for (i in 1:(num_of_bins-1)){
    sim[[paste0(let[j],let[(j+num_of_bins_per_group-1)])]]=let[j:(j+num_of_bins_per_group-1)]
    color[j:(j+num_of_bins_per_group-1)]=colors[i]
    j=j+1
  }
}

a <- list()
max_group_length <- c()
for (i in 1:num_of_topics){
  a[[i]] <- sim
  max_group_length <- c(max_group_length,length(a[[i]]))
}

sim <- a
max_group_length <- max(max_group_length)

#altsim = list("F","W",c("A","I","L","V"),c("M","C"),"P","G","Y",c("T","S"),c("H","K","R"),c("E","D"),c("Q","N"))
altsim <- sim
# Naming the group of similarities
for (i in 1:length(altsim)){
  names(altsim[[i]]) <- c("f","w","a","s","p","g","y","h","b","c","m")[1:(num_of_bins/num_of_bins_per_group)]
}

# Create custom colour scheme
cs1 = make_col_scheme(chars=let,
                      cols=c(color))

# Compute Letter Similarity 
letter_sim=as.data.frame(matrix(0,nrow = length(let),ncol = length(let)))
for (i in 1:length(let)){
  for (j in 1:length(let)){
    letter_sim[i,j]=1-abs(i-j)/length(let)
  }
}
colnames(letter_sim)=let
row.names(letter_sim)=let

# Topic similarity 
topic_similarity=read.csv(paste0("../",dataset,slash_for_topics,data_topic,"/Topic_Similarity.txt"),sep="\t",header = F)

############## Tranform Input Data ############## 
if (!file.exists(paste0("../",dataset,slash_for_topics,data_topic,"/udata_N_",as.numeric(args[1]),".txt"))){
  #read data
  if (!file.exists(paste0("../",dataset,slash_for_topics,data_topic,"/model-final.theta"))){
    temp <- read.csv(paste0("../",dataset,slash_for_topics,data_topic,"/", dataset_name, ".txt"),sep=",",header = F)
    actual_classes <- temp[,ncol(temp)]
    topics_folder <- paste0(num_of_topics,"topics")
    write.table(temp[,1:(ncol(temp)-1)],paste0("../","Data","/",dataset_name,"/",topics_folder,"/model-final.theta"),sep=" ",row.names = F,col.names = F)
    write.table(1:nrow(temp),paste0("../","Data","/",dataset_name,"/",topics_folder,"/docIDs.dat"),sep=" ",row.names = F, col.names = F)
    write.table(actual_classes,paste0("../","Data","/",dataset_name,"/",topics_folder,"/actual_classes.txt"),sep="\t",row.names = F, col.names = F)
    #write.table(ap_top_terms,paste0("Data","/",dataset_name,"/",topics_folder,"/ap_top_terms.txt"),sep="\t",row.names = F)
  }
  
  data=read.csv(paste0("../",dataset,slash_for_topics,data_topic,"/model-final.theta"),sep=" ",header = F)[,1:num_of_topics]
  
  if (file.exists(paste0("../",dataset,slash_for_topics,data_topic,"/docIDs.dat"))){
    data=cbind(sequence_id=read.csv(paste0("../",dataset,slash_for_topics,data_topic,"/docIDs.dat"),header = F)[,],data)
  }else{
    data=cbind(sequence_id=1:nrow(data),data)
  }
  
  data[,2:ncol(data)]=data[,2:ncol(data)]*10000/max(data[,2:ncol(data)]*10000)
  
  data_new=data
  data_discr=data
  
  for (i in 1:(length(groups)-1)){
    ids=which(data[,2:ncol(data)]>=groups[i] & data[,2:ncol(data)]<groups[i+1],arr.ind=TRUE)
    for (j in unique(ids[,2])){
      data_new[ids[which(ids[,2]==j),1],j+1]=let[i]
      data_discr[ids[which(ids[,2]==j),1],j+1]=groups[i]
    }
    if (i==length(groups)-1){
      ids=which(data[,2:ncol(data)]==groups[i+1],arr.ind=TRUE)
      for (j in unique(ids[,2])){
        data_new[ids[which(ids[,2]==j),1],j+1]= let[i] #let[i+1]
        data_discr[ids[which(ids[,2]==j),1],j+1]= groups[i]#let[i+1]
      }
    }
  }
  
  data_new$x <- apply( data_new[ , 2:ncol(data_new) ] , 1 , paste , collapse = "" )
  udata <- data_new[c(1,ncol(data_new))]
  colnames(udata)=c("Sequence.ID","AA.JUNCTION")
  
  
}else{
  udata <- read.csv(paste0("../",dataset,slash_for_topics,data_topic,"/udata_N_",as.numeric(args[1]),".txt"),sep="\t",header = T)
}
