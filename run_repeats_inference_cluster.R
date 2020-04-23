#script to run simulations
#this script takes one tree produced by gillespieSC.R and 
#simulates random mutations 
#a variety of phylogenetic algorithms can then be used to reconstruct the tree
#the reconstruction is then assessed
rm(list = ls())

options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
tsk = as.integer(args[1])
tr = as.integer(args[2])
i = as.integer(args[3])
j = as.integer(args[4])
md.scen = as.integer(args[5])
wd = as.character(args[6])
set.seed(tsk)
setwd(wd)

library(stringr)
library(ape)
library(phangorn)
library(phytools)
library(geiger)
#in paper, results are shown for trees: 1, 5 and 7 and 9
tree_string = c("","2","_death01","_2death01","_death02","_2death02",
                "_death03","_2death03","_diff_dep")
load(paste("true_tree_rands_gil",tree_string[tr],".RData",sep=""))

#run function to produce all required variables
source("make_phytree_array_cluster.R")
source("add_MS_mutations_function_cluster.R")
source("reconstruct_tree_cluster.R")
source("compare_trees_cluster.R")
source("make_external_data_files.R")
load(paste("phytree_variables",tree_string[tr],".RData",sep=""))

mut_rates = c(0.001,0.0001)
#mut_rates = (10^(log10(mut_rates)+1) - 10^(log10(mut_rates)-1))/(2*log(10))
spr = 0 #spread of mutation rates
Ms = c(50,100,500,1000,5000,100000)
ini_len = 10

#sample final leaves
samp_rate = 1
if(samp_rate!=1){
  Ntips = length(phyTree$tip.label)
  tips = 1:Ntips
  tips_discard = phyTree$tip.label[sample(tips,ceiling((1-samp_rate)*Ntips),replace=FALSE)]
  phyTree_sample = drop.tip(phyTree,tip=as.character(tips_discard),trim.internal=TRUE,subtree=FALSE)
  children_alive = lapply(children_alive,function(x) setdiff(x,tips_discard))
  children_alive = unique(children_alive[sapply(children_alive,function(x) length(x)>1)])
}
true.tree = read.tree(file=paste("true_tree",tree_string[tr],".new",sep=""))
#outputs for each repeat
if( (i==1 && j==1) || (i==1 && md.scen==2) ){
  if(md.scen>1){
    num_mut = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.pars = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.njeu = list()
    perRec.bal = list()
    perRec.bayesFK = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.mlFK = array(0,dim=c(length(mut_rates),length(Ms)))
    recTree = list()
    timings = list()
    for(x in 1:2){
      recTree[[x]] = list()
      timings[[x]] = list()
      for(y in 1:5){
        recTree[[x]][[y]] = list()
        timings[[x]][[y]] = list()
      }
    }
    for(k in 1:6){
      perRec.njeu[[k]] = array(0,dim=c(length(mut_rates),length(Ms)))
      perRec.bal[[k]] = array(0,dim=c(length(mut_rates),length(Ms)))
    }
  }else{
    num_mut = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.pars = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.njeu = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.njma = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.njcos = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.baleu = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.balma = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.balcos = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.bayesFK = array(0,dim=c(length(mut_rates),length(Ms)))
    perRec.mlFK = array(0,dim=c(length(mut_rates),length(Ms)))
    recTree = list()
    timings = list()
    for(x in 1:2){
      recTree[[x]] = list()
      timings[[x]] = list()
      for(y in 1:5){
        recTree[[x]][[y]] = list()
        timings[[x]][[y]] = list()
      }
    }
  }
}else{
  load(paste("results",tree_string[tr],"_task_",tsk,"_md_",md.scen,".RData",sep=""))
}

#simulate MS length changes
add_MS_op = add_more_complex_MS_mutations(Ms[j],mut_rates[i],tree_string,tr,spr,ini_len)
num_mut[i,j] = add_MS_op$count_mut
MS_to_cluster = add_MS_op$MS_to_cluster
if(exists("tips_discard")){
  MS_to_cluster =
    MS_to_cluster[which(!(rownames(MS_to_cluster) %in% as.character(tips_discard))),]
}

#add missing data if required
if(md.scen>=3){
  #for each cell, simulate number of MS to be unobserved, then randomly set to NA
  missing.prop = (md.scen-2)/10
  Nna = ceiling(missing.prop*ncol(MS_to_cluster))*rep(1,nrow(MS_to_cluster))
  NAcells = lapply(Nna,function(x) sample(ncol(MS_to_cluster),x,replace=FALSE))
  for(k in 1:nrow(MS_to_cluster)){
    MS_to_cluster[k,NAcells[[k]]] = NA
  }
}

#check if there are any changes
if(any(apply(MS_to_cluster,2,function(x) length(unique(x[!is.na(x)]))>1))){
  trees = dist_recons(MS_to_cluster,ini_len,md.scen,tsk,i,j)
  phyTree$node.label = NULL
  phyTree$edge.length = runif(length(phyTree$tip.label)*2-2)
  timings[[i]][[j]][[2]] = trees$t.NJ
  timings[[i]][[j]][[3]] = trees$t.bal
  if(md.scen>1){
    recTree[[i]][[j]][[2]] = list()
    recTree[[i]][[j]][[3]] = list()
    for(k in 1:6){
      njtree = paste("treeNJeu",k,sep="")
      baltree = paste("treeBal",k,sep="")
      recTree[[i]][[j]][[2]][[k]] = trees[[njtree]]
      recTree[[i]][[j]][[3]][[k]] = trees[[baltree]]
      op.njeu = compare_clades(trees[[njtree]],children_alive)
      perRec.njeu[[k]][i,j] = op.njeu$percRec
      
      op.bal = compare_clades(trees[[baltree]],children_alive)
      perRec.bal[[k]][i,j] = op.bal$percRec
    }
  }else{
    recTree[[i]][[j]][[2]] = trees$treeNJeu
    recTree[[i]][[j]][[3]] = trees$treeNJma
    recTree[[i]][[j]][[4]] = trees$treeNJcos
    recTree[[i]][[j]][[5]] = trees$treeBaleu
    recTree[[i]][[j]][[6]] = trees$treeBalma
    recTree[[i]][[j]][[7]] = trees$treeBalcos
    op.njeu = compare_clades(trees$treeNJeu,children_alive)
    perRec.njeu[i,j] = op.njeu$percRec
    op.njma = compare_clades(trees$treeNJma,children_alive)
    perRec.njma[i,j] = op.njma$percRec
    op.njcos = compare_clades(trees$treeNJcos,children_alive)
    perRec.njcos[i,j] = op.njcos$percRec
    
    op.baleu = compare_clades(trees$treeBaleu,children_alive)
    perRec.baleu[i,j] = op.baleu$percRec
    op.balma = compare_clades(trees$treeBalma,children_alive)
    perRec.balma[i,j] = op.balma$percRec
    op.balcos = compare_clades(trees$treeBalcos,children_alive)
    perRec.balcos[i,j] = op.balcos$percRec
  }
  #make files for revbayes and paup
  makerevbayes(MS_to_cluster,tsk,md.scen)
  makepaup(MS_to_cluster,tsk,md.scen)
  #makerbcont(MS_to_cluster,tsk,md.scen)
}else{
  timings[[i]][[j]][[2]] = NA
  timings[[i]][[j]][[3]] = NA
  if(md.scen>1){
    for(k in 1:6){
      perRec.njeu[[k]][i,j] = NA
      perRec.bal[[k]][i,j] = NA
    }
  }else{
    perRec.njeu[i,j] = NA
    perRec.baleu[i,j] = NA
    perRec.njma[i,j] = NA
    perRec.balma[i,j] = NA
    perRec.njcos[i,j] = NA
    perRec.balcos[i,j] = NA
  }
  write.table(0,file=paste("nst_",tsk,"_md_",md.scen,".txt",sep=""),quote=FALSE,
              row.names = FALSE,col.names = FALSE)
}

if(md.scen>1){
  save(num_mut,recTree,timings,mut_rates,Ms,md.scen,
       perRec.pars,perRec.njeu,perRec.bal,perRec.bayesFK,perRec.mlFK,
       file=paste("results",tree_string[tr],"_task_",tsk,"_md_",md.scen,".RData",sep=""))
}else{
  save(num_mut,recTree,timings,mut_rates,Ms,md.scen,
       perRec.pars,perRec.njeu,perRec.njma,perRec.njcos,perRec.baleu,perRec.balma,perRec.balcos,
       perRec.bayesFK,perRec.mlFK,
       file=paste("results",tree_string[tr],"_task_",tsk,"_md_",md.scen,".RData",sep=""))
}

quit("no")

