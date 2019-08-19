run_FBC <- function(datapath, domain, algo, dataset_name, pipeline, endper_binary_tree, exclude_col_AA_binary_tree, exclude_col_end_AA_binary_tree, 
                    exclude_J6_AA_binary_tree,thr_branch_breaking, use_only_leafs_graph, thr_graph_connections, cores){
  
  save_rData=T
  load_rData=F
  
  endper <- endper_binary_tree
  algocol <- exclude_col_AA_binary_tree
  backcol <- exclude_col_end_AA_binary_tree
  backcolj6 <- exclude_J6_AA_binary_tree
  use_only_leafs <- use_only_leafs_graph

  # Read data
  udata <- read.csv(file = paste0(datapath, "/", dataset_name, "/udata.txt"), sep = "\t", stringsAsFactors = F)
  #print(head(udata))
  print(domain)
  if (domain == "items"){
    topic_similarity=read.csv(paste0(datapath, "/", dataset_name,"/Topic_Similarity.txt"),sep="\t",header = F)
    actual_classes=read.csv(paste0(datapath, "/", dataset_name,"/true_classes.txt"), stringsAsFactors = T, header = F)$V1
  }
  
  load(paste0(datapath, "/", dataset_name, "/params/let.rData"))
  load(paste0(datapath, "/", dataset_name, "/params/letter_sim.rData"))
  load(paste0(datapath, "/", dataset_name, "/params/cs1.rData"))
  load(paste0(datapath, "/", dataset_name, "/params/max_group_length.rData"))
  load(paste0(datapath, "/", dataset_name, "/params/sim.rData"))
  load(paste0(datapath, "/", dataset_name, "/params/altsim.rData"))

  num_of_topics <- (str_length(udata$AA.JUNCTION[1]))
  num_of_bins <- length(let)
  
  slash_for_topics <- ""
  data_topic <- ""
  output_folder <- paste0("output/", dataset_name)
  if(!file.exists(paste0("output"))){ 
    dir.create(paste0("output"))
  }
  
  # Evaluation metrics -> depends on the domain
  evaluation_metrics <- c("ff", "Num_of_clusters", "Average-Identity", "Average-Similarity", "Average-Entropy_Sim",
                          "Average-BS", "Average-TopicSim","SD-Similarity", "SD-Identity","Average-Entropy_Id", 
                          "SD-Entropy_Id", "SD-Entropy_Sim", "SD-BS", "SD-TopicSim")
  
  N <- nrow(udata)
  
  if (N > 10000){
    source('/home/analysis/binary_tree_construction_parallel_with_results.R', local=TRUE)
  }else{
    source('/home/analysis/binary_tree_construction.R',local=TRUE)
  }
  
  if (2 %in% pipeline){
    source('/home/analysis/branch_breaking_cwl.R', local=TRUE)
  }
  
  if (3 %in% pipeline){
    source('/home/analysis/distances_and_graph.R', local=TRUE)
  }
  
}

