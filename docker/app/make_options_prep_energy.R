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
    default="Datasets/01_Baseline",
    help=paste0(
      "The directory where the dataset is located."
    )
  ),
  make_option(
    opt_str=c("-b","--dataset_name"),
    action="store",
    default="01_Baseline",
    help="The name of the dataset."
  ),
  make_option(
    opt_str=c("-c","--pipeline"),
    action="store",
    default="1,2,3,4,6",
    help=paste0(
      "Pipeline options:\n",
      "1. Covert negative values to positive ones (Reversed polarization of smart meters). \n",
      "2. Delete the days for which we have no measurements \n",
      "3. Delete the days that have only one load value (Stuck metering equipment) \n",
      "4. Find which days are not close to the min,max,average (outlier values) \n",
      "5. Continues same values (Stuck metering equipment) \n",
      "6. Find big peaks (outlier values)\n",
      "Use comma to seperate the list of options (default %default)"
    )
  ),
  make_option(
    opt_str=c("-d","--large_peak"),
    action="store",
    default="10000",
    help=paste0(
      "Big peak values.\n",
      "(default %default)"
    )
  )
);

opt <- parse_args(OptionParser(option_list=option_list));

source("run_prep_energy.R")
run_prep_energy(datapath=opt$datapath, dataset_name=opt$dataset_name, pipeline=as.numeric(strsplit(opt$pipeline,",")[[1]]), 
                large_peak=opt$large_peak)

# Print help
# print_help( OptionParser(option_list=option_list))

