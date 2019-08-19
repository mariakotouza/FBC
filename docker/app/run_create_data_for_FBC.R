run_create_data_for_FBC <- function(datapath, domain, dataset_name, num_of_bins, num_of_bins_per_group,num_of_topics){
  num_of_topics <- as.numeric(num_of_topics)
  dir.create("output_create_udata")
  dir.create(paste0("output_create_udata/", dataset_name))
  dir.create(paste0("output_create_udata/",dataset_name,"/params"))
  
  if (domain %in%  c('items', 'time-series', 'numeric')){
    #group values using ranges 0:1 with step=1/num_of_bins
    groups=seq(0,1,1/num_of_bins)
    #to do: make it more general ~ do not use letters that are only 24, use combinations
    let = LETTERS[1:num_of_bins] # The letters matrix
    let=let[length(let):1]
    sim=list()
    color=matrix(0,nrow=1, ncol=length(let))
    j=1
    overlaping <- F
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
    
    write.table(x = read.csv(paste0(datapath, "/", dataset_name, "/udata.txt"), sep="\t"), file = paste0("output_create_udata/",dataset_name,"/udata.txt"), sep = "\t", row.names = F)
    write.table(x = read.csv(paste0(datapath, "/", dataset_name,"/Topic_Similarity.txt"), sep="\t"), file = paste0("output_create_udata/", dataset_name,"/Topic_Similarity.txt"), sep = "\t", row.names = F)
    write.table(x = read.csv(paste0(datapath, "/", dataset_name,"/true_classes.txt"), sep="\t"), file = paste0("output_create_udata/", dataset_name,"/true_classes.txt"), sep = "\t", row.names = F)
    
    #file.copy(paste0(datapath, "/", dataset_name, "/udata.txt"), paste0("output_create_data/", dataset_name), overwrite = TRUE)
  }
  
  if (domain == 'AA'){
    sim = list("F","W",c("A","I","L","V"),c("M","C"),"P","G","Y",c("T","S"),c("H","K","R"),c("E","D"),c("Q","N"))
    # Naming the group of similarities
    names(sim) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
    
    # A table with the letters
    let = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y") # The letters matrix
    
    # Create custom colour scheme
    cs1 = make_col_scheme(chars=c("F","W","A","I","L","V","M","C","P","G","Y","T","S","H","K","R","E","D","Q","N"),
                          cols=c("#1E90FF", "#BA55D3", "#0000FF", "#0000FF", "#0000FF", "#0000FF", "#C6E2FF", "#C6E2FF", "#FFD700", "#00EE00", "#C1FFC1", "#54FF9F", "#54FF9F", "#FF0000", "#FF0000", "#FF0000", "#FFD700", "#FFD700", "#ED9121", "#ED9121"))
    
    a <- list()
    max_group_length <- c()
    for (i in 1:num_of_topics){
      a[[i]] <- sim
      max_group_length <- c(max_group_length,length(a[[i]]))
    }
    
    sim <- a
    max_group_length <- max(max_group_length)
    
    altsim <- sim
    # Naming the group of similarities
    for (i in 1:length(altsim)){
      names(altsim[[i]]) <- c("f","w","a","s","p","g","y","h","b","c","m")
    }
    
    # Compute Letter Similarity 
    letter_sim=as.data.frame(matrix(0,nrow = length(let),ncol = length(let)))
    for (i in 1:length(let)){
      for (j in 1:length(let)){
        letter_sim[i,j]=0.5
      }
    }
    colnames(letter_sim)=let
    row.names(letter_sim)=let
    
    num_of_bins <- length(let)
    
    udata <- read.csv(paste0(datapath, "/filter_in_", dataset_name, ".txt"), sep = "\t", stringsAsFactors = F) %>% 
      filter(Summary.CDR3.IMGT.length == (num_of_topics - 2)) %>%
      select(Summary.Sequence.ID, Summary.AA.JUNCTION)
    colnames(udata) <- c('Sequence.ID', 'AA.JUNCTION')
    
    write.table(udata, file = paste0("output_create_udata/",dataset_name,"/udata.txt"), sep = "\t", row.names = F)
  }
  
  
  tmp_path <- getwd()
  
  save(let,file=paste0("output_create_udata/",dataset_name,"/params/let.rData"))
  save(letter_sim,file=paste0("output_create_udata/",dataset_name,"/params/letter_sim.rData"))
  save(cs1,file=paste0("output_create_udata/",dataset_name,"/params/cs1.rData"))
  save(max_group_length,file=paste0("output_create_udata/",dataset_name,"/params/max_group_length.rData"))
  save(sim,file=paste0("output_create_udata/",dataset_name,"/params/sim.rData"))
  save(altsim,file=paste0("output_create_udata/",dataset_name,"/params/altsim.rData")) 
  
}