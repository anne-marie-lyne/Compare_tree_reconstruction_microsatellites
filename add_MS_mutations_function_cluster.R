add_MS_mutations_function = function(M,mut_rate,tree_string,tr,mut_rate_spr,ini_len){
  #stepwise mutation model
  load(paste("true_tree_rands_gil",tree_string[tr],".RData",sep=""))
  load(paste("phytree_variables",tree_string[tr],".RData",sep=""))
  
  #parameters for microsatellite mutations
  ###########################################################################
  l = sample(seq(10,20,1), M, replace=T, prob = 0.3*0.7^seq(10,20,1))
  MS = array(0,dim=c(sum(reac_rands<=3)*2+1,M))
  MS[1,] = l
  MS_site = 1:M
  count_mut = 0
  mut_rate_sep = 10^{log10(mut_rate)*matrix(1,nrow=1,ncol=M)+runif(M,min=-mut_rate_spr,max=mut_rate_spr)}
  #mut_rate_sep = mut_rate + runif(M,min=-mut_rate,max=mut_rate)
  #mut_rate_sep[mut_rate_sep<mut_rate] = 0
  ###########################################################################
  ix=1
  for(i in 1:t_max) {
    #use inputted tree
    reac = reac_rands[i]
    #depending on reaction, choose a node to update
    node = node_rands[i]
    #if there was a division, add a mutation with prob mut_rate
    if(reac<=3){
      MS[(2*ix):(2*ix+1),] = MS[rep(node,2),]
      #then sim a mutation
      mut_site = MS_site[runif(M)<mut_rate_sep]
      if(length(mut_site)!=0){
        #SHOULD REALLY CHECK HERE WHETHER CURRENT LENGTH IS 0
        p_or_m = 2*(runif(length(mut_site))<0.5)-1
        cell = (runif(length(mut_site))<0.5)
        MS[cbind((2*ix)+cell,mut_site)] = MS[cbind((2*ix)+cell,mut_site)] + p_or_m
      }
      count_mut = count_mut + length(mut_site)
      ix = ix+1
    }
  } 

  #then take these rows from the MS
  MS_to_cluster = MS[alive_leaves,]
  #add stutter error at N-1 with probability 0.0024
  rnds = matrix(runif(ncol(MS_to_cluster)*nrow(MS_to_cluster)),nrow = nrow(MS_to_cluster))
  MS_to_cluster[rnds<0.0024] = MS_to_cluster[rnds<0.0024] - 1
  rownames(MS_to_cluster) = alive_leaves
  op = list()
  op$MS_to_cluster = MS_to_cluster
  op$count_mut = count_mut
  return(op)
}


add_more_complex_MS_mutations = function(M,mut_rate,tree_string,tr,mut_rate_spr,ini_len){
  load(paste("true_tree_rands_gil",tree_string[tr],".RData",sep=""))
  load(paste("phytree_variables",tree_string[tr],".RData",sep=""))
  
  #parameters for microsatellite mutations
  ###########################################################################
  #sample from geometric distribution on 10:20 to get initial lengths
  l = sample(seq(10,20,1), M, replace=T, prob = 0.3*0.7^seq(10,20,1))
  MS = array(0,dim=c(sum(reac_rands<=3)*2+1,M))
  MS[1,] = l
  MS_site = 1:M
  count_mut = 0
  mut_rate_sep = 10^{log10(mut_rate)*matrix(1,nrow=1,ncol=M)+runif(M,min=-mut_rate_spr,max=mut_rate_spr)}
  p = 0.8 #probability that length change is of size one (given mutation occurs)
  q = 0.8 #geometric parameter given that length change is bigger than 1
  len_pdf = c(q,q*(1-q),q*(1-q)^2)/(1-(1-q)^3)
  len_cdf = c(len_pdf[1],len_pdf[1]+len_pdf[2],1)
  #mut_rate_sep = mut_rate + runif(M,min=-mut_rate,max=mut_rate)
  #mut_rate_sep[mut_rate_sep<mut_rate] = 0
  ###########################################################################
  ix=1
  for(i in 1:t_max) {
    #use inputted tree
    reac = reac_rands[i]
    #depending on reaction, choose a node to update
    node = node_rands[i]
    #if there was a division, add a mutation with prob mut_rate
    if(reac<=3){
      MS[(2*ix):(2*ix+1),] = MS[rep(node,2),]
      #then sim a mutation
      mut_site = MS_site[runif(M)<mut_rate_sep]
      if(length(mut_site)!=0){
        #SHOULD REALLY CHECK HERE WHETHER CURRENT LENGTH IS 0
        #first check if mutation is plus or minus
        p_or_m = 2*(runif(length(mut_site))<0.5)-1
        #then check which cell inherits change
        cell = (runif(length(mut_site))<0.5)
        #then check if mutation is size one or not
        si = (runif(length(mut_site))<p)
        if(sum(si==0)!=0){
          #then simulate size change if length change is greater than one
          rands = runif(sum(si==0))
          mat = 
            matrix(rands,nrow=length(len_cdf),ncol=length(rands),byrow=TRUE)<matrix(len_cdf,nrow=length(len_cdf),ncol=length(rands),byrow=FALSE)
          si[si==0] = apply(mat,2,function(x) min(which(x==1))) + 1
        }
        ch = si*p_or_m 
        MS[cbind((2*ix)+cell,mut_site)] = MS[cbind((2*ix)+cell,mut_site)] + ch
      }
      count_mut = count_mut + length(mut_site)
      ix = ix+1
    }
  } 
  
  #then take these rows from the MS
  MS_to_cluster = MS[alive_leaves,]
  #add stutter error at N-1 with probability 0.0024
  rnds = matrix(runif(ncol(MS_to_cluster)*nrow(MS_to_cluster)),nrow = nrow(MS_to_cluster))
  MS_to_cluster[rnds<0.0024] = MS_to_cluster[rnds<0.0024] - 1
  rownames(MS_to_cluster) = alive_leaves
  op = list()
  op$MS_to_cluster = MS_to_cluster
  op$count_mut = count_mut
  return(op)
}


add_more_complex_MS_mutations_autosome = function(M,mut_rate,tree_string,tr,mut_rate_spr,ini_len){
  load(paste("true_tree_rands_gil",tree_string[tr],".RData",sep=""))
  load(paste("phytree_variables",tree_string[tr],".RData",sep=""))
  
  #parameters for microsatellite mutations
  ###########################################################################
  #sample from geometric distribution on 10:20 to get initial lengths
  l1 = sample(seq(10,20,1), M/2, replace=T, prob = 0.3*0.7^seq(10,20,1))
  l2 = c(l1[1:floor(length(l1)/2)], (l1[(floor(length(l1)/2)+1):length(l1)] + 
           rgeom(length(l1[(floor(length(l1)/2)+1):length(l1)]),0.4) +1))
  l = c(matrix(c(l1,l2), byrow = T, nrow = 2))
  MS = array(0,dim=c(sum(reac_rands<=3)*2+1,M))
  MS[1,] = l
  MS_site = 1:M
  count_mut = 0
  mut_rate_sep = 10^{log10(mut_rate)*matrix(1,nrow=1,ncol=M)+runif(M,min=-mut_rate_spr,max=mut_rate_spr)}
  p = 0.8 #probability that length change is of size one (given mutation occurs)
  q = 0.8 #geometric parameter given that length change is bigger than 1
  len_pdf = c(q,q*(1-q),q*(1-q)^2)/(1-(1-q)^3)
  len_cdf = c(len_pdf[1],len_pdf[1]+len_pdf[2],1)
  #mut_rate_sep = mut_rate + runif(M,min=-mut_rate,max=mut_rate)
  #mut_rate_sep[mut_rate_sep<mut_rate] = 0
  ###########################################################################
  ix=1
  for(i in 1:t_max) {
    #use inputted tree
    reac = reac_rands[i]
    #depending on reaction, choose a node to update
    node = node_rands[i]
    #if there was a division, add a mutation with prob mut_rate
    if(reac<=3){
      MS[(2*ix):(2*ix+1),] = MS[rep(node,2),]
      #then sim a mutation
      mut_site = MS_site[runif(M)<mut_rate_sep]
      if(length(mut_site)!=0){
        #SHOULD REALLY CHECK HERE WHETHER CURRENT LENGTH IS 0
        #first check if mutation is plus or minus
        p_or_m = 2*(runif(length(mut_site))<0.5)-1
        #then check which cell inherits change
        cell = (runif(length(mut_site))<0.5)
        #then check if mutation is size one or not
        si = (runif(length(mut_site))<p)
        if(sum(si==0)!=0){
          #then simulate size change if length change is greater than one
          rands = runif(sum(si==0))
          mat = 
            matrix(rands,nrow=length(len_cdf),ncol=length(rands),byrow=TRUE)<matrix(len_cdf,
                                                                                    nrow=length(len_cdf),
                                                                                    ncol=length(rands),byrow=FALSE)
          si[si==0] = apply(mat,2,function(x) min(which(x==1))) + 1
        }
        ch = si*p_or_m 
        MS[cbind((2*ix)+cell,mut_site)] = MS[cbind((2*ix)+cell,mut_site)] + ch
      }
      count_mut = count_mut + length(mut_site)
      ix = ix+1
    }
  } 
  
  #then take these rows from the MS
  MS_to_cluster = MS[alive_leaves,]
  #add stutter error at N-1 with probability 0.0024
  rnds = matrix(runif(ncol(MS_to_cluster)*nrow(MS_to_cluster)),nrow = nrow(MS_to_cluster))
  MS_to_cluster[rnds<0.0024] = MS_to_cluster[rnds<0.0024] - 1
  #for each pair of alleles, put shortest in first position
  MS_to_cluster = do.call("cbind", lapply(seq(ncol(MS_to_cluster)/2), function(x) {
    rel.cols = MS_to_cluster[,((x-1)*2+1):(x*2)]
    mins = apply(rel.cols, 1, function(y) min(y))
    maxs = apply(rel.cols, 1, function(y) max(y))
    return(cbind(mins,maxs))
  }))
  rownames(MS_to_cluster) = alive_leaves
  op = list()
  op$MS_to_cluster = MS_to_cluster
  op$count_mut = count_mut
  return(op)
}



add_asymmetric_MS_mutations = function(M,mut_rate,tree_string,tr,mut_rate_spr,ini_len){
  load(paste("true_tree_rands_gil",tree_string[tr],".RData",sep=""))
  load(paste("phytree_variables",tree_string[tr],".RData",sep=""))
  
  #parameters for microsatellite mutations
  ###########################################################################
  l = sample(seq(10,20,1), M, replace=T, prob = 0.3*0.7^seq(10,20,1)) #initial number of tandem repeats in each MS
  MS = array(0,dim=c(sum(reac_rands<=3)*2+1,M))
  MS[1,] = l
  MS_site = 1:M
  count_mut = 0
  mut_rate_sep = 10^{log10(mut_rate)*matrix(1,nrow=1,ncol=M)+runif(M,min=-mut_rate_spr,max=mut_rate_spr)}
  p = 0.8 #probability that length change is of size one (given mutation occurs)
  q = 0.8 #geometric parameter given that length change is bigger than 1
  len_pdf = c(q,q*(1-q),q*(1-q)^2)/(1-(1-q)^3)
  len_cdf = c(len_pdf[1],len_pdf[1]+len_pdf[2],1)
  #mut_rate_sep = mut_rate + runif(M,min=-mut_rate,max=mut_rate)
  #mut_rate_sep[mut_rate_sep<mut_rate] = 0
  ###########################################################################
  ix=1
  for(i in 1:t_max) {
    #use inputted tree
    reac = reac_rands[i]
    #depending on reaction, choose a node to update
    node = node_rands[i]
    #if there was a division, add a mutation with prob mut_rate
    if(reac<=3){
      MS[(2*ix):(2*ix+1),] = MS[rep(node,2),]
      #then sim a mutation
      mut_site = MS_site[runif(M)<mut_rate_sep]
      if(length(mut_site)!=0){
        #SHOULD REALLY CHECK HERE WHETHER CURRENT LENGTH IS 0
        #first check if mutation is plus or minus
        p_or_m = 2*(runif(length(mut_site))<0.4)-1
        #then check which cell inherits change
        cell = (runif(length(mut_site))<0.5)
        #then check if mutation is size one or not
        si = (runif(length(mut_site))<p)
        if(sum(si==0)!=0){
          #then simulate size change if length change is greater than one
          rands = runif(sum(si==0))
          mat = 
            matrix(rands,nrow=length(len_cdf),ncol=length(rands),byrow=TRUE)<matrix(len_cdf,nrow=length(len_cdf),ncol=length(rands),byrow=FALSE)
          si[si==0] = apply(mat,2,function(x) min(which(x==1))) + 1
        }
        ch = si*p_or_m 
        MS[cbind((2*ix)+cell,mut_site)] = MS[cbind((2*ix)+cell,mut_site)] + ch
      }
      count_mut = count_mut + length(mut_site)
      ix = ix+1
    }
  } 
  
  #then take these rows from the MS
  MS_to_cluster = MS[alive_leaves,]
  #add stutter error at N-1 with probability 0.0024
  rnds = matrix(runif(ncol(MS_to_cluster)*nrow(MS_to_cluster)),nrow = nrow(MS_to_cluster))
  MS_to_cluster[rnds<0.0024] = MS_to_cluster[rnds<0.0024] - 1
  rownames(MS_to_cluster) = alive_leaves
  op = list()
  op$MS_to_cluster = MS_to_cluster
  op$count_mut = count_mut
  return(op)
}









