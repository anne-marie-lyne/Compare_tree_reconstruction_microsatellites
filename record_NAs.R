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
set.seed(tsk)
tree_string = c("","2","_death01","_2death01","_death02","_2death02",
                "_death03","_2death03","_diff_dep")
load(paste("results",tree_string[tr],"_task_",tsk,"_md_",md.scen,".RData",sep=""))

perRec.bayesFK[i,j] = NA
perRec.mlFK[i,j] = NA
perRec.pars[i,j] = NA

timings[[i]][[j]][[1]] = NA
timings[[i]][[j]][[4]] = NA
  
  
if(md.scen>1){
  save(num_mut,recTree,timings,
       perRec.pars,perRec.njeu,perRec.bal,perRec.bayesFK,perRec.mlFK,
       mut_rates,Ms,md.scen,file=paste("results",tree_string[tr],"_task_",tsk,
                                       "_md_",md.scen,".RData",sep=""))
}else{
  save(num_mut,recTree,timings,
       perRec.pars,perRec.njeu,perRec.njma,perRec.njcos,perRec.baleu,perRec.balma,perRec.balcos,
       perRec.bayesFK,perRec.mlFK,
       mut_rates,Ms,md.scen,file=paste("results",tree_string[tr],"_task_",tsk,
                                       "_md_",md.scen,".RData",sep=""))
}
quit("no")
