#script to go through each results file and compile results
Ntsk = 60
tr=1
md.scen=1
Model = "_tree1"
setwd(paste0("/data/tmp/run_repeats_tree_inference",md.scen,Model))

tree_string = c("","2","_death01","_2death01","_death02","_2death02",
                "_death03","_2death03","_diff_dep")
load(paste("results",tree_string[tr],"_task_1_md_",md.scen,".RData",sep=""))

#make matrices to hold results
RecTree = list()
if(md.scen>1){
  PerRec.njeu = list()
  PerRec.bal = list()
  for(i in 1:6){
    PerRec.njeu[[i]] = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
    PerRec.bal[[i]] = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
  }
}else{
  PerRec.njeu = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
  PerRec.njma = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
  PerRec.njcos = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
  PerRec.baleu = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
  PerRec.balma = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
  PerRec.balcos = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
}
Num_mut = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
PerRec.pars = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
PerRec.bayesFK = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
PerRec.mlFK = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))

Timings = list()
for(i in 1:4){
  Timings[[i]] = array(0,dim=c(length(mut_rates),length(Ms),Ntsk))
}

if(md.scen>1){
  for(i in 1:Ntsk){
    load(paste("results",tree_string[tr],"_task_",i,"_md_",md.scen,".RData",sep=""))
    Num_mut[,,i] = num_mut
    PerRec.pars[,,i] = perRec.pars
    PerRec.bayesFK[,,i] = perRec.bayesFK
    PerRec.mlFK[,,i] = perRec.mlFK
    for(k in 1:6){
      PerRec.njeu[[k]][,,i] = perRec.njeu[[k]]
      PerRec.bal[[k]][,,i] = perRec.bal[[k]]
    }
    RecTree[[i]] = recTree
  }
}else{
  for(i in 1:Ntsk){
    load(paste("results",tree_string[tr],"_task_",i,"_md_",md.scen,".RData",sep=""))
    Num_mut[,,i] = num_mut
    PerRec.pars[,,i] = perRec.pars
    PerRec.njeu[,,i] = perRec.njeu
    PerRec.njma[,,i] = perRec.njma
    PerRec.njcos[,,i] = perRec.njcos
    PerRec.baleu[,,i] = perRec.baleu
    PerRec.balma[,,i] = perRec.balma
    PerRec.balcos[,,i] = perRec.balcos
    PerRec.bayesFK[,,i] = perRec.bayesFK
    PerRec.mlFK[,,i] = perRec.mlFK
    RecTree[[i]] = recTree
  }
}
print(length(RecTree))
print(length(RecTree[[1]]))
print(length(RecTree[[1]][[1]]))
print(length(RecTree[[1]][[1]][[1]]))
print(length(RecTree[[1]][[1]][[1]][[1]]))

for(i in 1:Ntsk){
  load(paste("results",tree_string[tr],"_task_",i,"_md_",md.scen,".RData",sep=""))
  for(j in 1:length(mut_rates)){
    for(k in 1:(length(Ms)-1)){
      for(l in 1:4){
        Timings[[l]][j,k,i] = timings[[j]][[k]][[l]]
      }
    }
  }
}

Model = ""
if(md.scen>1){
  save(Num_mut,
       PerRec.pars,PerRec.njeu,PerRec.bal,PerRec.bayesFK,PerRec.mlFK,
       mut_rates,Ms,Ntsk,md.scen,RecTree,Timings,
       file=paste("full_results",tree_string[tr],"_bayes_md_",md.scen,Model,".RData",sep=""))
}else{
  save(Num_mut,
       PerRec.pars,PerRec.njeu,PerRec.njma,PerRec.njcos,PerRec.baleu,PerRec.balma,PerRec.balcos,
       PerRec.bayesFK,PerRec.mlFK,
       mut_rates,Ms,Ntsk,md.scen,RecTree,Timings,
       file=paste("full_results",tree_string[tr],"_bayes_md_",md.scen,Model,".RData",sep=""))
}

quit("no")