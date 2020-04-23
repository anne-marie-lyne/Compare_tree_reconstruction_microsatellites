#script to compute distances for paup parsimony tree
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

library(stringr)
library(ape)
library(phangorn)
library(phytools)
library(geiger)
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

#read in tree outputted by paup
t.pars.paup = read.nexus(paste("pars_paup_output",tsk,"_md_",md.scen,".nex",sep=""))


#use first tree
#remove 't's from taxon names
tiplabels = t.pars.paup$tip.label
t.pars.paup$tip.label = substring(tiplabels,2)
recTree[[i]][[j]][[1]] = t.pars.paup
rf.tree = multi2di(t.pars.paup,random=TRUE)
write.tree(rf.tree, file = paste("pars_paup_output",tsk,"_md_",md.scen,".new",sep=""))
op.pars = compare_clades(rf.tree,children_alive)
perRec.pars[i,j] = op.pars$percRec


timings[[i]][[j]][[1]] = read.table(paste0(wd,"/timing_",tsk,"_md_",md.scen,".txt"), 
                            header = F, colClasses = "numeric")[1,1]


if(md.scen>1){
  save(num_mut,recTree,timings,
       perRec.pars,perRec.njeu,perRec.bal,perRec.bayesFK,perRec.mlFK,
       mut_rates,Ms,md.scen,file=paste("results",tree_string[tr],"_task_",tsk,"_md_",md.scen,".RData",sep=""))
}else{
  save(num_mut,recTree,timings,
       perRec.pars,perRec.njeu,perRec.njma,perRec.njcos,perRec.baleu,perRec.balma,perRec.balcos,
       perRec.bayesFK,perRec.mlFK,
       mut_rates,Ms,md.scen,file=paste("results",tree_string[tr],"_task_",tsk,"_md_",md.scen,".RData",sep=""))
}
quit("no")




