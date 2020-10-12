#script with function to make phylo object

make_phytree_array = function(tree_string,tr){
  
  load(paste("true_tree_rands_gil",tree_string[tr],".RData",sep=""))
  utri = tree
  utri[lower.tri(utri,diag=TRUE)] = 0
  Nleaves = sum(rowSums(utri)==0)
  Nnodes = ncol(tree)
  
  #create orig_node_names, element i is original name of node i
  new_node_names = 1:ncol(tree)
  orig_node_names = c()
  orig_node_names[1:Nleaves] = new_node_names[rowSums(utri)==0]
  internal_nodes = new_node_names[rowSums(utri)!=0]
  orig_node_names[(Nleaves+1):Nnodes] = internal_nodes
  
  phyTree = list()
  phyTree$edge = matrix(0,nrow=2*(Nnodes-Nleaves),ncol=2)
  count = 1
  for(i in 1:(Nnodes-Nleaves)){
    #find row internal_nodes[i]
    #find children (these are off-diagonal columns which aren't empty)
    #find new number of these element (which elements of orig_node_names)
    ind = which(tree[internal_nodes[i],]!=0)
    phyTree$edge[count,] = c(which(orig_node_names==internal_nodes[i]),
                             which(orig_node_names==ind[2]))
    count = count + 1
    phyTree$edge[count,] = c(which(orig_node_names==internal_nodes[i]),
                             which(orig_node_names==ind[3]))
    count = count + 1
  }
  phyTree$tip.label = orig_node_names[1:Nleaves]
  phyTree$node.label = orig_node_names[(Nleaves+1):Nnodes]
  phyTree$Nnode = Nnodes-Nleaves
  class(phyTree) = "phylo"
  
  #find which nodes are alive
  alive_leaves = new_node_names[rowSums(utri)==0 & diag(tree)<=2]
  alive_leaves_new = sapply(alive_leaves,function(x) which(orig_node_names==x))
  
  #find children of each node in new numbering
  children_new = list()
  for(i in 1:length(children)){
    if(length(children[[i]])!=0){
      children_new[[i]] = 
        sapply(children[[i]],function(x) which(orig_node_names==x))
    }
  }
  
  #keep only the leaves
  children_new_alive = 
    lapply(children_new,function(x) intersect(x,alive_leaves_new))
  children_new_alive = 
    children_new_alive[sapply(children_new_alive,function(x) length(x)!=0)]
  children_new_alive = lapply(children_new_alive,sort)
  children_alive = lapply(children,function(x) intersect(x,alive_leaves))
  children_alive = 
    children_alive[sapply(children_alive,function(x) length(x)>1)]
  children_alive = lapply(children_alive,sort)
  
  #find dead leaves
  dead_leaves = new_node_names[diag(tree)==3]
  dead_leaves_new = sapply(dead_leaves,function(x) which(orig_node_names==x))
  
  if(length(dead_leaves)!=0){
    #prune tree
    phyTree = drop.tip(phyTree, tip=dead_leaves_new, trim.internal = TRUE, subtree = FALSE)
  }
  
  #save outputs to a file
  save(alive_leaves,children_alive,phyTree,orig_node_names,
       file=paste("phytree_variables",tree_string[tr],".RData",sep=""))
  return(phyTree)
}