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
f_treesFK = as.character(args[5])
md.scen = as.integer(args[6])
wd = as.character(args[7])
set.seed(tsk)

library(stringr)
library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(methods)
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

#outputs for each repeat
load(paste("results",tree_string[tr],"_task_",tsk,"_md_",md.scen,".RData",sep=""))

#read in mcmc trees GTK
mcmc_treesFK = read.table(f_treesFK,colClasses = c("character"),header=TRUE)
#find which column tree is in
c = which(colnames(mcmc_treesFK)=="psi")
write.table(mcmc_treesFK[which.max(as.numeric(mcmc_treesFK[,3])),c],
            file=paste("ml_tree_",tsk,"_md_",md.scen,".txt",sep=""),
            quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(mcmc_treesFK[which.max(as.numeric(mcmc_treesFK[,2])),c],
            file=paste("map_tree_",tsk,"_md_",md.scen,".txt",sep=""),
            quote=FALSE,col.names=FALSE,row.names=FALSE)
treeML = read.tree(paste("ml_tree_",tsk,"_md_",md.scen,".txt",sep=""))
treeML$root.edge = NULL
treeBayes = read.tree(paste("map_tree_",tsk,"_md_",md.scen,".txt",sep=""))
treeBayes$root.edge = NULL
write.tree(treeML, file = paste("ml_tree_",tsk,"_md_",md.scen,".new",sep=""))
write.tree(treeBayes, file = paste("map_tree_",tsk,"_md_",md.scen,".new",sep=""))

#compute distances
op.mlFK = compare_clades(treeML,children_alive)
perRec.mlFK[i,j] = op.mlFK$percRec


op.bayesFK = compare_clades(treeBayes,children_alive)
perRec.bayesFK[i,j] = op.bayesFK$percRec


#read in time taken for MCMC
timings[[i]][[j]][[4]] = read.table(paste0(wd,"/timing_",tsk,"_md_",md.scen,".txt"), 
                                    header = F, colClasses = "numeric")[1,1]

if(md.scen>1){
  recTree[[i]][[j]][[4]] = treeBayes
  recTree[[i]][[j]][[5]] = treeML
  save(num_mut,recTree,timings,
       perRec.pars,perRec.njeu,perRec.bal,perRec.bayesFK,perRec.mlFK,
       mut_rates,Ms,md.scen,file=paste("results",tree_string[tr],"_task_",tsk,"_md_",md.scen,".RData",sep=""))
}else{
  recTree[[i]][[j]][[8]] = treeBayes
  recTree[[i]][[j]][[9]] = treeML
  save(num_mut,recTree,timings,
       perRec.pars,perRec.njeu,perRec.njma,perRec.njcos,perRec.baleu,perRec.balma,perRec.balcos,
       perRec.bayesFK,perRec.mlFK,
       mut_rates,Ms,md.scen,file=paste("results",tree_string[tr],"_task_",tsk,"_md_",md.scen,".RData",sep=""))
}
quit("no")


