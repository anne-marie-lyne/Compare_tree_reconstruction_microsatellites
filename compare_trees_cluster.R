robinson_foulds = function(treeRec,phyTree){
  rf = treedist(treeRec,phyTree,check.labels = TRUE)
  return(rf[1])
}

compare_clades = function(treeRec,children_alive){
  #extract all clades by going through internal nodes
  Ntip = length(treeRec$tip.label)
  alltips = as.numeric(treeRec$tip.label)
  children_rec = list()
  j = 1
  for(i in (Ntip+2):(Ntip+treeRec$Nnode)){
    children_rec[[j]] = as.numeric(tips(treeRec,i)) #NB this is in original numbers
    children_rec[[j+1]] = setdiff(alltips,children_rec[[j]])
    j = j+2
  }
  for(i in 1:Ntip){
    children_rec[[j]] = setdiff(alltips,alltips[i])
    j = j+1
  }
  #then go through each clade in true tree and see if it is present in rec tree
  sizeClades = numeric(length=(length(children_alive)-1))
  cladeFound = vector(length=(length(children_alive)-1))
  for(i in 2:length(children_alive)){
    sizeClades[i-1] = length(children_alive[[i]])
    cladeFound[i-1] = any(sum(sapply(children_rec,
           function(x) length(x)==length(children_alive[[i]]) & 
             length(intersect(x,children_alive[[i]]))==length(children_alive[[i]]))))
  }
  
  #find percentage of each size reconstructed
  sizeUnique = sort(unique(sizeClades))
  numRec = sapply(sizeUnique,function(x) sum(cladeFound[sizeClades==x]))
  percRec = sapply(sizeUnique,function(x) numRec[sizeUnique==x]/sum(sizeClades==x))
  op = list()
  op$percRecSize = percRec
  op$sizeUnique = sizeUnique
  op$percRec = 1 - sum(numRec)/(length(children_alive)-1)
  op$numRec = numRec
  op$numSize = sapply(sizeUnique,function(x) sum(sizeClades==x))
  return(op)
}