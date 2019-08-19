#args <- commandArgs(TRUE)
#args <- c("100", "Hitech", "1", "0")
#source("preprocessing_items_terminal.R")

library("profvis") 
library(inline)
includes <- '#include <sys/wait.h>'
code <- 'int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};'
wait <- cfunction(body=code, includes=includes, convention='.C')

#remove(data, data_discr)
save_rData <- F

if (Sys.info()[1] == "Windows"){
  cores <- 1
}else{
  if (as.numeric(args[4])==0){
    cores <- detectCores(all.tests = FALSE, logical = TRUE)
  }else{
    cores <- as.numeric(args[4])
  }
}

############################### Initialization ############################### 
#change!!
#endper = 100 # Set the percentage, which ends the programm 
#exclude_letters <- c()  #c("group_1")
#absoluteIdentity <- T

d = 0

udata$AA.JUNCTION = as.character(udata$AA.JUNCTION)
udata$clusters = 0 # Initialiaze the column clusters with 0
udata$level.0 = 0 # Initialize the column of cl.0 with 0
ff <- data.frame()
perlist <- list()
persimlist <- list()

#logfile
out <- paste0(output_folder,slash_for_topics,data_topic,"/log")
if(!file.exists(paste0(output_folder))){ 
  dir.create(paste0(output_folder))
}
if(!file.exists(paste0(output_folder,slash_for_topics,data_topic))){ 
  dir.create(paste0(paste0(output_folder,slash_for_topics,data_topic)))
}
if(!file.exists(paste0(out))){ 
  dir.create(paste0(out))
}
logFile = paste0(output_folder,slash_for_topics,data_topic,"/log/log_file ",trunc(as.numeric(Sys.time())),"_cores_",cores,".txt")
cat(paste0("Level","\t","Cluster","\t","N","\t","Duration","\t","Memory"), file=logFile, append=FALSE, sep = "\n")

if(!file.exists(paste0(output_folder,slash_for_topics,data_topic,"/RDATA"))){ 
  dir.create(paste0(output_folder,slash_for_topics,data_topic,"/RDATA"))
}

if(!file.exists(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo))){ 
  dir.create(paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo))
}

#################### Run ###########################
#prof_FHC <- profvis({
en1 = tic()
level = -1
cl_at_level <- 0 #contains the clusters ids as they are saved in the structure
all_leafs <- F
while (!all_leafs) {
  level <- level + 1
  next_level_cl <- c()
  all_leafs <- T
  templ2=mclapply(cl_at_level,mc.cores = cores, mc.preschedule = TRUE, function(cluster_id) {
                    #lista <- tree_construction_cluster_i(cluster_id,level,let,sim,d,algo,algocol,backcol,backcolj6,logFile)
                    br <- cluster_id
                    ids_br <- which(udata$clusters == br)
                    vv <- length(ids_br)
                    ep <- level
                    tic()
                    #print(paste0("cl ",br," has ",length(which(udata$clusters == br))," sequences"))
                    if (domain!="AA"){
                      mymat = matrix(0,nrow=length(let), ncol=num_of_topics)
                      permat = matrix(0,nrow=length(let) + 1, ncol=num_of_topics)
                      simmat = matrix(0,nrow = max_group_length,ncol =num_of_topics)
                      persim = matrix(0,nrow = max_group_length+1,ncol =num_of_topics)
                      rownames(mymat) = let
                      rownames(permat) = c(let,"Entropy")
                      rownames(simmat) = c(paste0("group_",1:max_group_length))
                      rownames(persim) = c(paste0("group_",1:max_group_length),"Entropy")
                    }else{
                      # Find the sequences with gene J6
                      ind1 = str_which(udata$J.GENE.and.allele[ids_br],"J6")
                      qw = 1:length(udata$AA.JUNCTION[ids_br])
                      if(length(ind1) == 0){
                        ind2 = qw
                      }else{
                        ind2 = qw[-ind1]
                      }
                      
                      # Initialize matrices for J6 sequences
                      mymat1 = matrix(0,nrow=length(let), ncol=num_of_topics)
                      simmat1 = matrix(0,nrow = max_group_length,ncol =num_of_topics)
                      rownames(mymat1) = let
                      rownames(simmat1) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
                      
                      # Initialize matrices for NO J6 sequences
                      mymat2 = matrix(0,nrow=length(let), ncol=num_of_topics)
                      simmat2 = matrix(0,nrow = max_group_length,ncol =num_of_topics)
                      rownames(mymat2) = let
                      rownames(simmat2) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
                      
                      # Initialize matrices for all sequences
                      permat = matrix(0,nrow=length(let) + 1, ncol=num_of_topics)
                      rownames(permat) = c(let,"Entropy")
                      persim = matrix(0,nrow = max_group_length+1,ncol =num_of_topics)
                      rownames(persim) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am","Entropy")
                      
                      # Compute matrices for J6 sequences
                      if(length(ind1) != 0 ){
                        trimudata = strsplit(udata[ids_br,]$AA.JUNCTION[ind1],"")
                        align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
                        for(i in (algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6)){
                          temptab = plyr::count(align[i],vars = colnames(align)[i])
                          names(temptab)[1] = "X1"
                          mymat1[which(is.na(match(names(mymat1[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
                          temptab1 = temptab
                          gg = as.vector(unlist(temptab1[1]))
                          temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),gg[x]))])))
                          #temptab1[1] = names(sim[[i]])[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim,x)))]
                          temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
                          simmat1[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim[[i]]) == x)),i] = as.vector(unlist(temptab1[2]))
                        }
                      }
                      
                      # Compute matrices for NO J6 sequences
                      if(length(ind2) != 0 ){
                        trimudata = strsplit(udata[ids_br,]$AA.JUNCTION[ind2],"")
                        align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
                        for(i in (algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcol)){
                          temptab = plyr::count(align[i],vars = colnames(align)[i])
                          names(temptab)[1] = "X1"
                          mymat2[which(is.na(match(names(mymat2[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
                          permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
                          temptab1 = temptab
                          gg = as.vector(unlist(temptab1[1]))
                          temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim[[i]])[which(str_detect(unlist(sim[[i]]),gg[x]))]))) #delete digits from the names of similarity groups
                          temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
                          simmat2[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim[[i]]) == x)),i] = as.vector(unlist(temptab1[2]))
                          persim[max_group_length+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
                        }
                      }
                      
                      # Combine matrices for No J6 and J6 sequences
                      mymat = mymat1 + mymat2
                      simmat = simmat1 + simmat2
                      meg = length(ind1) + length(ind2)
                      if(length(ind2) == 0 || backcol == backcolj6){
                        permat[1:length(let),(algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6)] = (mymat[,(algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
                        persim[1:max_group_length,(algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6)] = (simmat[,(algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
                      }else{
                        permat[1:length(let),(algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6)] = (mymat[,(algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
                        persim[1:max_group_length,(algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6)] = (simmat[,(algocol+1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6)] / meg) * 100
                        if(backcolj6 != 0){
                          permat[1:length(let),(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcol)] = (mymat[,(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcol)] / length(ind2)) * 100
                          persim[1:max_group_length,(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcol)] = (simmat[,(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcolj6 + 1):(str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcol)] / length(ind2)) * 100
                        }
                      }
                    }
                    
                    
                    # Find the Entropy for all the sequences
                    align = strsplit(udata$AA.JUNCTION[ids_br],"")
                    align = data.frame(matrix(unlist(align), nrow=length(align), byrow=T))
                    
                    for(i in 1:num_of_topics){
                      temptab = plyr::count(align[i],vars = colnames(align)[i])
                      names(temptab)[1] = "X1"
                      match_elements=match(names(mymat[,i]),as.vector(unlist(temptab[1])))
                      ids=as.vector(unlist(temptab[1]))[na.omit(match_elements)]
                      mymat[ids,i] = as.vector(unlist(temptab[2]))[na.omit(match_elements)]
                      #mymat[which(is.na(match(names(mymat[,i]),as.vector(unlist(temptab[1])))) == FALSE),i] = as.vector(unlist(temptab[2]))
                      permat[length(let) + 1,i] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
                      temptab1 = temptab
                      temptab1[1] = names(sim[[i]])[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim[[i]],x)))]
                      temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
                      
                      simmat[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim[[i]]) == x)),i] = as.vector(unlist(temptab1[2]))
                      persim[max_group_length+1,i] = entropy(simmat[,i],base = exp(1))
                    }
                    permat[1:length(let),] = (mymat / vv) * 100
                    persim[1:max_group_length,] = (simmat / vv) * 100
                    
                    if (domain=="AA"){
                      if(algocol != 0 && backcol!=0){
                        # Initialize matrices
                        mymat3 = matrix(0,nrow=length(let), ncol=num_of_topics)
                        simmat3 = matrix(0,nrow = max_group_length,ncol =num_of_topics)
                        rownames(mymat3) = let
                        rownames(simmat3) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am")
                        permat3 = matrix(0,nrow=length(let) + 1, ncol=num_of_topics)
                        rownames(permat3) = c(let,"Entropy")
                        persim3 = matrix(0,nrow = max_group_length+1,ncol =num_of_topics)
                        rownames(persim3) = c("F","W","Al","Su","P","G","Y","Hy","Ba","Ac","Am","Entropy")
                        
                        trimudata = strsplit(udata[ids_br,]$AA.JUNCTION,"")
                        align = data.frame(matrix(unlist(trimudata), nrow=length(trimudata), byrow=T))
                        exc = 1:algocol
                        exc2 = (str_length(udata[ids_br,]$AA.JUNCTION[1]) - backcol):str_length(udata[ids_br,]$AA.JUNCTION[1])
                        exc = c(exc,exc2)
                        for(i in 1:length(exc)){
                          temptab = plyr::count(align[exc[i]],vars = colnames(align)[exc[i]])
                          names(temptab)[1] = "X1"
                          mymat3[which(is.na(match(names(mymat3[,exc[i]]),as.vector(unlist(temptab[1])))) == FALSE),exc[i]] = as.vector(unlist(temptab[2]))
                          permat3[length(let) + 1,exc[i]] = entropy(xtabs(freq ~ X1, temptab),base=exp(1))
                          temptab1 = temptab
                          gg = as.vector(unlist(temptab1[1]))
                          temptab1[1] = sapply(1:length(gg),function(x) gsub('[[:digit:]]+', '', names(unlist(sim)[which(str_detect(unlist(sim),gg[x]))])))
                          #temptab1[1] = names(sim[[i]])[sapply(as.vector(unlist(temptab1[1])),function(x) which(str_detect(sim,x)))]
                          temptab1 = aggregate(temptab1[2], by=temptab1[1], FUN=sum)
                          simmat3[sapply(as.vector(unlist(temptab1[1])),function(x) which(names(sim[[i]]) == x)),exc[i]] = as.vector(unlist(temptab1[2]))
                          persim3[max_group_length+1,i] = entropy(xtabs(freq ~ X1, temptab1),base = exp(1))
                        }
                        permat3[1:length(let),] = (mymat3 / length(udata[ids_br,]$AA.JUNCTION)) * 100
                        persim3[1:max_group_length,] = (simmat3 / length(udata[ids_br,]$AA.JUNCTION)) * 100
                        permat[1:length(let),exc] = permat3[1:length(let),exc] 
                        persim[1:max_group_length,exc] = persim3[1:max_group_length,exc] 
                      }
                    }
                    
                    permat_m <- permat
                    
                    ################################ Finish ###########################
                    t1 =which(permat[1:(nrow(permat)-1),] == 100,arr.ind = TRUE)
                    sumper = (length(as.numeric(t1[,2]))* 100) / num_of_topics
                    #find the number of non-zero groups. if this number is = 2 than sim=100
                    t2 =which(persim[1:(nrow(persim)-1),] >= 100,arr.ind = TRUE)
                    
                    sumper2 = (length(as.numeric(t2[,2]))* 100) / num_of_topics
                    
                    #dfsum[nrow(dfsum) + 1,] = c(sumper,sumper2,br,vv)
                    leaf2 <- F
                    if(algo == "Identity"){
                      if ((sumper >= endper) | (vv==1)){
                        nn = TRUE  # When nn = TRUE the percentage of sumper < endper%
                        leaf2 <- T
                        en = toc(quiet = TRUE)
                        cat(paste0(level,"\t",br,"\t",vv,"\t",en$toc - en$tic,"\t",(pryr::mem_used())/10^6), file=logFile, append=TRUE, sep = "\n")
                      }
                    }else{
                      if ((sumper2 >= endper) | (vv==1)){
                        leaf2 = TRUE	
                        en = toc(quiet = TRUE)
                        cat(paste0(level,"\t",br,"\t",vv,"\t",en$toc - en$tic,"\t",(pryr::mem_used())/10^6), file=logFile, append=TRUE, sep = "\n")
                      }
                    }
                  
                    #result1 = list("permat"= permat, "persim" = persim,"sumper" = sumper,"sumper2" = sumper2, "leaf2" = leaf2)
                    #result1 = list("leaf2" = leaf2)
                    if (!leaf2) {
                      
                      ####################  Choice ####################
                      
                      
                      ################## Find desired cell for division ################
                      if(algo == "Similarity" ){
                        permat = persim # If sumper < endper% we want to check only the persim matrix
                      }
                      
                      if (domain!="AA"){
                        permat_temp <- permat
                      }else{
                        permat_temp <- permat[,(algocol + 1):(ncol(permat)-backcol)]
                      }
                      
                      cel = which(permat_temp == max(permat_temp), arr.ind = TRUE)
                      selected_cell_id = 1
                      poss = max(permat_temp) 
                      # We exclude the 100 % from the max values
                      if (max(permat_temp) == 100){ 
                        cel = which(permat_temp == max(permat_temp[permat_temp!=max(permat_temp)]), arr.ind = TRUE) # The desired cell
                        poss = max(permat_temp[permat_temp!=max(permat_temp)])
                      }
                      
                      # Check if clusters division is stopped or not
                      if(poss != 0){  
                        # if cel contain more than one cells, find the best cell matching some criteria
                        selected_cell_id = 1
                        if(algo == "Identity"){
                          if ((length(cel)/2) > 1){
                            dddff = min(permat_temp[nrow(permat_temp),cel[,2]])
                            dddff2 = which(permat_temp[nrow(permat_temp),cel[,2]] == dddff)
                            selected_cell_id = dddff2[1]
                            if (length(dddff2)>1 && nn == FALSE){ # If the vector has 2 or more numbers means that we have columns with the same entropy and nn = FALSE in order not to double check the persim
                              dddff3 = persim[,(algocol + 1):(ncol(permat)-backcol)][sapply(1:length(dddff2), function (x){ str_which(names(sim[[selected_cell_id]]),gsub('[[:digit:]]+', '', names(unlist(sim[[selected_cell_id]])[which(str_detect(unlist(sim[[selected_cell_id]]),let[cel[x,1]]))])))}),cel[,2]]
                              dddff3 = max(diag(dddff3))
                              dddff4 = which(diag( persim[,(algocol + 1):(ncol(permat)-backcol)][sapply(1:length(dddff2), function (x){str_which(names(sim[[selected_cell_id]]),gsub('[[:digit:]]+', '', names(unlist(sim[[selected_cell_id]])[which(str_detect(unlist(sim[[selected_cell_id]]),let[cel[x,1]]))])))}),cel[,2]]) == dddff3)
                              selected_cell_id = dddff4[1]
                              if(length(dddff4)> 1){
                                dddff5 = min(persim[,(algocol + 1):(ncol(permat)-backcol)][nrow(persim[,(algocol + 1):(ncol(permat)-backcol)]),cel[dddff4,2]])
                                dddff6 = which(persim[,(algocol + 1):(ncol(permat)-backcol)][nrow(persim[,(algocol + 1):(ncol(permat)-backcol)]),cel[dddff4,2]] == dddff5)
                                selected_cell_id = dddff6[1] 
                              }
                            }
                          }
                        }else{
                          if ((length(cel)/2) > 1){
                            dddff = min(permat_temp[nrow(permat_temp),cel[,2]])
                            dddff2 = which(permat_temp[nrow(permat_temp),cel[,2]] == dddff)
                            selected_cell_id = dddff2[1]
                          }
                        }
                        
                        ################### Divide #####################
                        # If we need a new level, then we create a new column in udata dataframe with its name (level.ep)
                        
                        if (algo == "Identity"){
                          # Find sequences that contain the cell's letter in cell's position 
                          x1 = str_which(str_detect(str_sub(udata$AA.JUNCTION[ids_br],(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), let[cel[selected_cell_id,1]]), "TRUE")
                          length_cl1 = length(x1)
                          cl1 = (2*(br+1)-1)
                          
                          # The other sequences of the cluster
                          x2 = str_which(str_detect(str_sub(udata$AA.JUNCTION[ids_br],(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), let[cel[selected_cell_id,1]]), "FALSE")
                          length_cl2 = length(x2)
                          cl2 = (2*(br+1)-1) + 1
                        }else{
                          # Find sequences contains the cell's similarity group in cell's position 
                          strings.to.find = unlist(sim[[selected_cell_id]][cel[selected_cell_id,1]])
                          x1 = str_which(str_detect(str_sub(udata$AA.JUNCTION[ids_br],(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), str_c(strings.to.find, collapse="|")), "TRUE")
                          length_cl1 = length(x1)
                          cl1 = (2*(br+1)-1)
                          
                          # The other sequences of the cluster
                          x2 = str_which(str_detect(str_sub(udata$AA.JUNCTION[ids_br],(cel[selected_cell_id,2]+algocol),(cel[selected_cell_id,2]+algocol)), str_c(strings.to.find, collapse="|")), "FALSE")
                          length_cl2 = length(x2)
                          cl2 = (2*(br+1)-1) + 1
                        }
                        # Find which new sub-cluster has more sequences and give it first cluster name
                        if(length_cl1 < length_cl2){
                          temp2 = x1
                          x1 = x2
                          x2 = temp2
                        }
                        
                      }
                      en = toc(quiet = TRUE)
                      cat(paste0(level,"\t",br,"\t",vv,"\t",en$toc - en$tic,"\t",(pryr::mem_used())/10^6), file=logFile, append=TRUE, sep = "\n")
                    }
                    
                    if (!leaf2) {
                      next_level_cl <- c((2*(cluster_id+1)-1),2*(cluster_id+1))
                      cl <- c(cl1,cl2)
                      #result <- list("next_level_cl"=next_level_cl,"leaf2"=leaf2,"x1"=x1, "x2"=x2, "cl"=cl)
                      result <- list("next_level_cl"=next_level_cl,"leaf2"=leaf2,"x1"=x1, "x2"=x2, "cl"=cl,"ff"=c(sumper,sumper2,cluster_id,vv,level),
                                     "persim"=persim, "permat_m"=permat_m)
                    }else{
                      #result <- list("next_level_cl"=c(),"leaf2"=leaf2)
                      result <- list("next_level_cl"=c(),"leaf2"=leaf2,"ff"=c(sumper,sumper2,cluster_id,vv,level),
                                     "persim"=persim, "permat_m"=permat_m)
                    }
                    
  })
  #wait()
  all_leafs <- all(laply(templ2,.fun = function(x)  x$leaf2))
  if (!all_leafs){
    j=0
    for (cluster_id in cl_at_level){
      j=j+1
      if (!templ2[[j]]$leaf2){
        udata[udata[udata$clusters == cluster_id,]$Sequence.ID[templ2[[j]]$x1],sprintf('level.%d', level+1)] = templ2[[j]]$cl[1]
        udata[udata[udata$clusters == cluster_id,]$Sequence.ID[templ2[[j]]$x2],sprintf('level.%d', level+1)] = templ2[[j]]$cl[2]
        #cat(paste0(cluster_id,"\t",templ2[[j]]$cl[1], "\t", templ2[[j]]$cl[2]), file=logFile, append=TRUE, sep = "\n")
      }
    }
    udata$clusters <- apply(udata[,4:(ncol(udata))],1, function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))
  }
  ff <- rbind(ff,ldply(templ2,.fun = function(x)  x$ff))
  per_temp <- llply(templ2,.fun = function(x)  x$permat_m) 
  names(per_temp) <- sprintf('permat_br.%d', cl_at_level) # Save the permat with this format
  perlist$temp <- per_temp
  names(perlist)[length(perlist)] <- sprintf('permat_br.%d&', cluster_id) 
  per_temp <- llply(templ2,.fun = function(x)  x$persim)
  names(per_temp) = sprintf('persim_br.%d', cl_at_level)
  persimlist$temp <- per_temp
  names(persimlist)[length(persimlist)] <- sprintf('persim_br.%d&', cluster_id) 
  
  if (!all_leafs){ 
    cl_at_level <- unlist(llply(templ2,.fun = function(x)  x$next_level_cl))
    }
  
}

en2 = toc(quiet = TRUE)
cat(paste0((en2$toc - en1)), file=logFile, append=TRUE, sep = "\n")
cat(paste0("lastlist","\t",(en2$toc - en1)), file=logFile, append=TRUE, sep = "\n")

perlist <- unlist(perlist, recursive=FALSE)
persimlist <- unlist(persimlist, recursive=FALSE)

names(persimlist) <- ldply(strsplit(names(persimlist), '&.'))$V2
names(perlist) <- ldply(strsplit(names(perlist), '&.'))$V2

#})

#mem_max_FHC <- max(prof_FHC$x$message$prof$memalloc)
#mem_mean_FHC <- mean(prof_FHC$x$message$prof$memalloc)
#time_FHC <- max(prof_FHC$x$message$prof$time)
#num_of_items_FHC <- N
#cat(paste(max(prof_FHC$x$message$prof$memalloc), mean(prof_FHC$x$message$prof$memalloc),
#          max(prof_FHC$x$message$prof$time,(en2$toc - en1)), sep="\t"), file=out2, append=T, sep = "\n")

mem <- read.csv(logFile,sep="\t")$Memory
#mem <- mem %>% group_by(Level) %>% summarize(n=sum(Memory))
#print(mem$n)
stat_FHC2 <- cbind(N, mem_max = mean(na.omit(mem)), max_mem = max(na.omit(mem)),Time=(en2$toc - en1))
write.table(stat_FHC2,paste0(output_folder,slash_for_topics,data_topic,"/stat_FHC_parallel_N_",N,"_",as.numeric(args[3]),"_cores_",cores,".txt"),sep="\t",row.names = F)

udata[udata==0] <- NA
udata$level.0 <- 0
lev <- as.numeric(strsplit(names(udata)[ncol(udata)],"\\.")[[1]][2])

colnames(ff) <- c("sumper", "sumper2", "branch",  "len", "level")
initial_cl_ids <- ff$branch
ff$branch <- 0:(nrow(ff)-1)

# Change cluster ids to continues numbers at udata
udata_temp <- udata[,c(1,3)]
for (i in 1:(lev+1)){
  udata_temp <- rbind(udata_temp,data.frame(Sequence.ID=(udata$Sequence.ID+i*nrow(udata)),clusters=udata[[i+3]]))
}

udata_temp <- udata_temp[order(udata_temp$clusters),]
udata_temp$clusters <- as.numeric(as.factor(udata_temp$clusters)) - 1
udata_temp <- udata_temp[order(udata_temp$Sequence.ID),]

udata2 <- udata
udata2[[3]] <- udata_temp[1:nrow(udata),2]
for (i in 1:(lev+1)){
  udata2[[i+3]] <- udata_temp[(i*nrow(udata)+1):(nrow(udata)*(i+1)),2]
}

udata <- udata2
rm(udata2)
rm(udata_temp)

clep <- c()
for (i in 0:lev){
  temp <- length(unique(na.omit(udata[[sprintf('level.%d', i)]])))
  clep <- c(clep,rep(i,times=temp))
}

# save rData to disk
if (save_rData){
  save(listb,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/listb.rData"))
}

############## Compute Metrics for each cluster #####################
#lastlist = listb
#lastlist2 = listb
#lastlist3 = listb
# The final name of udata data frame 
df = udata
# A list with the permat matrix for every branch
#perlist = lastlist$list
#persimlist = lastlist$listn
# Create a dataframe with all clusters and their identity and similarity percentage 
#ff = lastlist$dfsum
ff$level <- clep  #without leaves
#clep <- clep

Clus = as.data.frame(matrix(100, ncol = 6, nrow = nrow(ff)))
names(Clus) = c("ClusterId","Identity","Similarity","Entropy_Id","Entropy_Sim","BS")
Clus$ClusterId = ff$branch
Clus$seqnum = ff$len
Clus$level = c(na.omit(clep))
Clus$Entropy_Id=0
Clus$Entropy_Sim=0

for(i in 1:length(ff$branch) ){
  ll = which(Clus$ClusterId == ff$branch[i])
  Clus$Identity[ll] = ff$sumper[i]
  Clus$Similarity[ll] = ff$sumper2[i]
  Clus$Entropy_Id[ll] = mean(perlist_temp[[i]][nrow(perlist_temp[[i]]),])
  Clus$Entropy_Sim[ll] = mean(persimlist_temp[[i]][nrow(persimlist_temp[[i]]),])
  Clus$BS[ll]=mean(colSums(computeSimilarityCl(perlist_temp[[i]],letter_sim,let)[1:(nrow(perlist_temp[[i]])-1),])) 
}

rm(perlist_temp, persimlist_temp)

##############  Find the pathString for each node/cluster ############## 
lev = max(na.omit(clep))
df = df[ do.call( order , df[ , match(  colnames(df[str_which(names(df), "level.")]) , names(df) ) ]  ) , ]
df_args <- c(df[str_which(names(df), "level.")], sep="/")
if(lev == max(na.omit(clep))){
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
  Clus$Topic1[i]=as.numeric(names(findTopicsPerCl(AminoCl(Clus$ClusterId[i],clep[2:length(clep)],df,flagtic,logFile),3,let))[1])
  Clus$Topic2[i]=as.numeric(names(findTopicsPerCl(AminoCl(Clus$ClusterId[i],clep[2:length(clep)],df,flagtic,logFile),3,let))[2])
  Clus$Topic3[i]=as.numeric(names(findTopicsPerCl(AminoCl(Clus$ClusterId[i],clep[2:length(clep)],df,flagtic,logFile),3,let))[3])
  Clus$Topic1_notA[i]=as.numeric(findTopicsPerCl(AminoCl(Clus$ClusterId[i],clep[2:length(clep)],df,flagtic,logFile),3,let)[1])
  Clus$Topic2_notA[i]=as.numeric(findTopicsPerCl(AminoCl(Clus$ClusterId[i],clep[2:length(clep)],df,flagtic,logFile),3,let)[2])
  Clus$Topic3_notA[i]=as.numeric(findTopicsPerCl(AminoCl(Clus$ClusterId[i],clep[2:length(clep)],df,flagtic,logFile),3,let)[3])
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
t1 = which(clep == num_of_levels)
# epipleon
xm = as.numeric(as.data.frame(xN$leaves))
orio = min(which(clep == num_of_levels))
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
flagtic <- T
metrics_per_level=EmPin(Clus,clep,xN,dfsd,flagtic,logFile)
#format(metrics_per_level, digits=3)
#metrics_per_level[,2:ncol(metrics_per_level)]=round(metrics_per_level[,2:ncol(metrics_per_level)], 3)
metrics_per_level[,3:ncol(metrics_per_level)]=apply(metrics_per_level[,3:ncol(metrics_per_level)], 2, function(x) formatC(x, format = "f", digits = 3))
write.table(metrics_per_level, paste0(output_folder,slash_for_topics,data_topic,"/metrics_per_level_whole_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)

#### Save useful variables and matrices that have come from the initial tree #### 
x_initial=x
df_initial=df
Clus_initial=Clus
clep_initial=clep
xN_initial <- xN
ff_initial <- ff

# save rData to disk
if (save_rData){
  out <- paste0(output_folder,slash_for_topics,data_topic,"/RDATA")
  if(!file.exists(out)){ 
    dir.create(out)
  }
  save(x_initial,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/x_initial.rData"))
  save(df_initial,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/df_initial.rData"))
  save(Clus_initial,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/Clus_initial.rData"))
  save(clep_initial,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/clep_initial.rData"))
  save(xN_initial,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/xN_initial.rData"))
  save(ff_initial,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/ff_initial.rData"))
  save(perlist,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/perlist.rData"))
  save(persimlist,file=paste0(output_folder,slash_for_topics,data_topic,"/RDATA/",algo,"/persimlist.rData"))
}


#write.table(Clus_initial,paste0(output_folder,slash_for_topics,data_topic,"/Clus_initial_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)
#write.table(df_initial,paste0(output_folder,slash_for_topics,data_topic,"/df_initial_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)
#write.table(x_initial,paste0(output_folder,slash_for_topics,data_topic,"/x_initial_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)
#write.table(clep_initial,paste0(output_folder,slash_for_topics,data_topic,"/clep_initial_algo_",algo,".txt"),sep="\t",row.names = FALSE, col.names = TRUE)




