#!/usr/bin/env Rscript

# Wrapper for the run_FBC_without_ui script. It lives as an R shell script in 
# order to be seamlessly integrated to the general variant calling pipeline 
# which uses only system interactions.
# To run the script using the command line using the default input arguments: Rscript --vanilla make_options_FBC.R

# Author: Maria Kotouza

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(
    opt_str=c("-a","--datapath"),
    action="store",
    default="../Data/Hitech",
    help=paste0(
      "The directory where the dataset is located."
    )
  ),
  make_option(
    opt_str=c("-b","--domain"),
    action="store",
    default="items",
    help=paste0(
      "Use:\n",
      "'AA' for gene sequence data\n",
      "'items' for document data\n",
      "'time-series' for time-series data\n",
      "'numeric' for arithmetic data\n",
      "'fashion' for clothing data\n",
      "(default %default)"
    )
  ),
  make_option(
    opt_str=c("-b","--algo"),
    action="store",
    default="Similarity",
    help="Supported algorithms: One of 'Similarity' (default) or 'Identity'."
  ),
  make_option(
    opt_str=c("-c","--dataset_name"),
    action="store",
    default="Hitech",
    help="The name of the dataset."
  ),
  make_option(
    opt_str=c("-d","--pipeline"),
    action="store",
    default="1,2,3",
    help=paste0(
      "Pipeline options:\n",
      "1. Binary Tree Construction \n",
      "2. Branch Breaking \n",
      "3. Graph Construction \n",
      "4. Evaluation \n",
      "Use comma to seperate the list of options (default %default)"
    )
  ),
  make_option(
    opt_str=c("-e","--endper_binary_tree"),
    action="store",
    default="100",
    help=paste0(
      "Set the maximum identity value for the leaf clusters of the binary tree.\n",
      "(default %default)"
    )
  ),
  make_option(
    opt_str=c("-f","--exclude_col_AA_binary_tree"),
    action="store",
    default="0",
    help=paste0(
      "This is applicable only at the AA domain. \n",
      "Exclude (default %default) columns from the beginning of the sequences.\n"
    )
  ),
  make_option(
    opt_str=c("-g","--exclude_col_end_AA_binary_tree"),
    action="store",
    default="0",
    help=paste0(
      "This is applicable only at the AA domain. \n",
      "Exclude (default %default) columns from the end of the sequences.\n"
    )
  ),
  make_option(
    opt_str=c("-h","--exclude_J6_AA_binary_tree"),
    action="store",
    default="0",
    help=paste0(
      "This is applicable only at the AA domain. \n",
      "Exclude gene J6 (default %default).\n"
    )
  ),
  make_option(
    opt_str=c("-i","--endper_binary_tree"),
    action="store",
    default="100",
    help=paste0(
      "Set the maximum identity value for the leaf clusters of the binary tree.\n",
      "(default %default)"
    )
  ),
  make_option(
    opt_str=c("-j","--thr_branch_breaking"),
    action="store",
    default="0.01,0.01,0.01,0.01",
    help=paste0(
      "Select thresholds for the branch breaking algorithm for the metrics: \n",
      "Number of clusters' elements, Identity, Entropy, Bin Similarity, seperated by comma (default %default). \n",
      "Description: for each branch of the tree, the appropri-ate level to be cut is examined \n",
      "by recursively comparing the parent cluster with its two children using the values \n",
      "of the evaluation metrics and those thresholds"
    )
  ),
  make_option(
    opt_str=c("-k","--use_only_leafs_graph"),
    action="store",
    default="Yes",
    help=paste0(
      "Select 'Yes' if you want to create a graph only with the leaf clusters of the tree, or otherwise select 'No'\n",
      "(default %default)"
    )
  ),
  make_option(
    opt_str=c("-l","--thr_graph_connections"),
    action="store",
    default="40",
    help=paste0(
      "remove all those connections that have a dissimilarity percentage > thr % \n",
      "(default %default)"
    )
  ),
  make_option(
    opt_str=c("-m","--cores"),
    action="store",
    default="4",
    help=paste0(
      "Number of cores to be used. \n",
      "(default %default)"
    )
  )
);


opt <- parse_args(OptionParser(option_list=option_list));

if (opt$use_only_leafs_graph == "Yes"){
  use_only_leafs_graph <- T
}else{
  use_only_leafs_graph <- F
}

source("run_FBC_without_ui.R")
run_FBC(datapath=opt$datapath, domain=opt$domain, algo=opt$algo, dataset_name=opt$dataset_name, pipeline=as.numeric(strsplit(opt$pipeline,",")[[1]]), 
        endper_binary_tree=as.numeric(opt$endper_binary_tree), exclude_col_AA_binary_tree = as.numeric(opt$exclude_col_AA_binary_tree),
        exclude_col_end_AA_binary_tree=as.numeric(opt$exclude_col_end_AA_binary_tree), exclude_J6_AA_binary_tree=as.numeric(opt$exclude_J6_AA_binary_tree),
        thr_branch_breaking=as.numeric(strsplit(opt$thr_branch_breaking,",")[[1]]), 
        use_only_leafs_graph=use_only_leafs_graph, thr_graph_connections=as.numeric(opt$thr_graph_connections), cores=opt$cores)

# Print help
# print_help( OptionParser(option_list=option_list))

