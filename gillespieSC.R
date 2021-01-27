#script to carry out Gillespie algorithm to stochastically evolve a pool
#of stem cells and differentiated cells

#this version starts with one SC, initially expands the pool of stem cells via self renewal,
#and then allows them to differentiate or self renew, in two different phases
#random numbers are saved so that same tree can be resimulated for use
#with microsatellite simulations
#NB script was originally more general and allowed to start from any number of SC
#can't guarantee there won't be errors if you try this, particularly in saving...
#plus it will be really slow and input to microsatellite part assumes a tree...

#inputs are
# iniNum: vector [N_SC N_Diff] where N_SC and N_Diff are initial numbers of SC and diff cells
# simTime: scalar giving simulation time

#stem cells can either divide to give two stem cells, or to give two
#differentiated cells
#differentiated cells either divide to give two differentiated cells or die

rm(list=ls(all=TRUE))
library(ape)
library(phangorn)
library(phytools)
library(geiger)

#reactions
# SC -> SC + SC
# SC -> Diff + Diff
# Diff -> Diff + Diff
# Diff -> 0

#SET WORKING DIRECTOR HERE
setwd("PATH_TO/gillespie_output")

#Simulate a number of different trees with different params, 
#give each a different name
#output for tree 5 is provided on github
tr=5
tree_string = c("","2","_death01","_2death01","_death02","_2death02",
                "_death03","_2death03","_test")

#parameters for Gillespie algorithm determining which cells divide and when
###########################################################################
#code SC = 1, Diff = 2, Dead = 3
iniNum = c(1,0)
Ntree = sum(iniNum)
simTime = 25
stoich_mat = matrix(c(1,0,0,-1,2,0,0,1,0,0,-1,1),
                    nrow=4,ncol=3,byrow=TRUE)
#define probabilities of reaction for each equation
#initial phase, expansion of SC, no differentiation
probSC = c(1,0) #probabilities of each SC reaction
probDiff = c(0.5,0.5) #probabilities of each diff reaction
rateSC = 1 #SC reaction rate
rateDiff = 0 #diff reaction rate
cSC = probSC * rateSC
cDiff = probDiff * rateDiff
#second phase, diff cells divide x100 more frequently than SC
#for diff cells, symm division more likely than death
probSC2 = c(0.5,0.5)
probDiff2 = c(0.7,0.3)
rateSC2 = 0.01
rateDiff2 = 1
cSC2 = probSC2 * rateSC2
cDiff2 = probDiff2 * rateDiff2
#third phase, diff cells divide x100 more frequently than SC
#for diff cells, symm division and death equally likely
probSC3 = c(0.5,0.5)
probDiff3 = c(0.5,0.5)
rateSC3 = 0.01
rateDiff3 = 1
cSC3 = probSC3 * rateSC3
cDiff3 = probDiff3 * rateDiff3
#times for end of each phase
#careful, cell numbers will explode if t_1 and t_2 are too large
t_1 = 3.5
t_2 = 13
###########################################################################

#adjacency matrix with the type of each node and the nodes each turns into
Ndiv = vector("list", Ntree)
trees = list()
for(i in 1:Ntree){
  if(i<=iniNum[1]){
    trees[[i]] = matrix(1,1,1)
  } else{
    trees[[i]] = matrix(2,1,1)
  }
}

#save total number of each species at each iter, the children of each node, 
Ntype = matrix(0,nrow=sum(iniNum),ncol=2)
if(iniNum[1]>0){
  Ntype[c(1:iniNum[1]),1] = c(1:iniNum[1])
  if(iniNum[2]>0){
    Ntype[(iniNum[1]+1):sum(iniNum),1] = iniNum[1]
    Ntype[(iniNum[1]+1):sum(iniNum),2] = 1:iniNum[2]
  }
} else{
  Ntype[(iniNum[1]+1):sum(iniNum),2] = 1:iniNum[2]
}
spNum = matrix(c(iniNum,0),nrow=1,ncol=3)
spNum_save = spNum
children = list()
reac_rands = c()
node_rands = c()

t = 0
t_save = t
while(t<simTime){
  print(t)
  print(Ntype)
  #sim time tau and which reaction
  if(t<t_1){
    #phase 1
    lambda = sum(cSC*spNum[1] + cDiff*spNum[2])
    tau = rexp(1,rate=lambda)
    reac = which(rmultinom(1,1,
                           c(cSC*spNum[1],cDiff*spNum[2]))==1)
    reac_rands = append(reac_rands,reac)
  }else if(t>=t_1 && t<t_2){
    #phase 2
    lambda = sum(cSC2*spNum[1] + cDiff2*spNum[2])
    tau = rexp(1,rate=lambda)
    reac = which(rmultinom(1,1,
                           c(cSC2*spNum[1],cDiff2*spNum[2]))==1)
    reac_rands = append(reac_rands,reac)
  }else{
    #phase 3
    lambda = sum(cSC3*spNum[1] + cDiff3*spNum[2])
    tau = rexp(1,rate=lambda)
    reac = which(rmultinom(1,1,
                           c(cSC3*spNum[1],cDiff3*spNum[2]))==1)
    reac_rands = append(reac_rands,reac)
  }
  #choose which node to update
  node_overall = floor(runif(1,min=1,
                             max=((reac<3)*spNum[1]+(reac>=3)*spNum[2])+1))
  #update overall cell numbers
  spNum = spNum + stoich_mat[reac,]
  spNum_save = rbind(spNum_save,spNum)
  
  #find which tree this node is in, and which node within tree
  #if reac<3 then need to count SC=1, if reac>=3 then need to count Diff=2 cells
  #each time a tree is changed, need to update that and all subsequent trees
  sp_to_count = (reac<3)*1 + (reac>=3)*2
  inds = which(((node_overall - Ntype[,sp_to_count])<=0)!=0,arr.ind=T)
  j = inds[1]
  ind = 1:nrow(trees[[j]])
  if(j==1){
    if(nrow(trees[[j]])==1){
      node_tree = 1
    } else{
      n = node_overall
      utri = trees[[j]]
      utri[lower.tri(utri,diag=TRUE)] = 0
      rel_ind = ind[diag(trees[[j]])==sp_to_count & 
                      rowSums(utri)==0]
      node_tree = rel_ind[n]      
    }
  } else{
    if(nrow(trees[[j]])==1){
      node_tree = 1
    } else{    
      n = node_overall - Ntype[j-1,sp_to_count]
      utri = trees[[j]]
      utri[lower.tri(utri,diag=TRUE)] = 0
      rel_ind = ind[diag(trees[[j]])==sp_to_count & 
                      rowSums(utri)==0]
      node_tree = rel_ind[n]
    }
  }
  node_rands = c(node_rands,node_tree)
  Ntype[j:Ntree,] = Ntype[j:Ntree,] + stoich_mat[rep(reac,Ntree-j+1),1:2] 
  
  #update adjacency matrix for relevant tree
  if(reac<=3){
    #copy matrix (if not reaction 4) and make new bigger matrix
    tree_copy = trees[[j]]
    trees[[j]] = matrix(0,nrow=nrow(tree_copy)+2,ncol=ncol(tree_copy)+2)
    trees[[j]][1:nrow(tree_copy),1:ncol(tree_copy)] = tree_copy   
  }
  #add rows for new daughter cells, and info about children
  sp_created = (reac==1)*1 + (reac==2 || reac==3)*2
  if(reac<=3){
    trees[[j]][nrow(tree_copy)+1,nrow(tree_copy)+1] = sp_created
    trees[[j]][nrow(tree_copy)+2,nrow(tree_copy)+2] = sp_created
    trees[[j]][node_tree,(nrow(tree_copy)+1):(nrow(tree_copy)+2)] = rep(sp_created,2)
    Ndiv[[j]] = append(Ndiv[[j]],rep(Ndiv[[j]][node_tree]+1,2))
    #also update children
    children[[node_tree]] = c(nrow(trees[[j]])-1,nrow(trees[[j]]))
    #find previous nodes which have 'node' as a child and add new children
    logi = sapply(children,function(x) sum(x==node_tree)>0)
    newchildren= lapply(children[logi],
                        function(x) append(x,c(nrow(trees[[j]])-1,nrow(trees[[j]]))))
    children[logi] = newchildren   
  } else{
    trees[[j]][node_tree,node_tree] = 3;
  }
  #update time
  t = t + tau
  t_save = append(t_save,t)
  #end if there are no SC or Diff cells left
  if(sum(spNum[1:2])==0){
    break
  }
}
tree = trees[[1]]
t_max = length(reac_rands)

#save output
save(tree,reac_rands,node_rands,children,t_max,
     file=paste("true_tree_rands_gil",tree_string[tr],".RData",sep=""))


# make phylo object required for later parts of simulation
# source("PATH_TO/make_phytree_array.R")
# phyTree = make_phytree_array(tree_string,tr)
# load(paste("phytree_variables",tree_string[tr],".RData",sep=""))
# plot(phyTree)

