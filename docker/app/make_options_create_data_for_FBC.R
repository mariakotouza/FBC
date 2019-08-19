#!/usr/bin/env Rscript

# Wrapper for the run_discretization script. It lives as an R shell script in 
# order to be seamlessly integrated to the general variant calling pipeline 
# which uses only system interactions.
# To run the script using the command line using the default input arguments: Rscript --vanilla make_options_discretization.R

# Author: Maria Kotouza

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option(
    opt_str=c("-a","--datapath"),
    action="store",
    default="Data/Hitech",
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
    opt_str=c("-c","--dataset_name"),
    action="store",
    default="Hitech",
    help="The name of the dataset."
  ),
  make_option(
    opt_str=c("-d","--num_of_bins"),
    action="store",
    default="10",
    help=paste0(
      "The number of bins selected for discretization (default %default) \n"
    )
  ),
  make_option(
    opt_str=c("-e","--num_of_bins_per_group"),
    action="store",
    default="2",
    help=paste0(
      "The number of consecutive bins that will be assigned to a group for the Similarity-based algorithm (default %default)"
    )
  ),
  make_option(
    opt_str=c("-f","--num_of_topics"),
    action="store",
    default="20",
    help=paste0(
      "The length of the vectors (default %default)"
    )
  )
);

opt <- parse_args(OptionParser(option_list=option_list));

source("run_create_data_for_FBC.R")
run_create_data_for_FBC(datapath=opt$datapath, domain=opt$domain, dataset_name=opt$dataset_name, 
                   num_of_bins=as.numeric(opt$num_of_bins), num_of_bins_per_group=as.numeric(opt$num_of_bins_per_group),
                   num_of_topics=opt$num_of_topics)

# Print help
# print_help( OptionParser(option_list=option_list))

