makerevbayes = function(MS_to_cluster,tsk,md.scen){
  modefunc <- function(x){
    tabresult <- tabulate(x)
    themode <- which(tabresult == max(tabresult))
    if(sum(tabresult == max(tabresult))>1) themode <- NA
    return(themode)
  }
  #produce nexus file for input to RevBayes
  #assume mode is the original and centre around this value
  modes = apply(MS_to_cluster, 2, modefunc)
  #replace NAs in each column with this value
  # MS_to_cluster = do.call("cbind",lapply(seq(length(modes)),
  #                                        function(x) {MS_to_cluster[which(is.na(MS_to_cluster[,x])),x]=modes[x]
  #                                        return(MS_to_cluster[,x])} ))
  MS_to_cluster2 = t(apply(MS_to_cluster,1,function(x) x-modes+1))
  minMS = min(MS_to_cluster2,na.rm = TRUE)
  MS_to_cluster2 = MS_to_cluster2 - minMS
  #take only the MS which show a change
  write.table((max(MS_to_cluster2,na.rm=TRUE)+1),file=paste("nst_",tsk,"_md_",md.scen,".txt",sep=""),
              quote=FALSE,row.names = FALSE,col.names = FALSE)
  logi = apply(MS_to_cluster2,2,function(x) length(unique(x[!is.na(x)]))>1)
  write.table(MS_to_cluster2[,logi],paste("ms_sim_data_rb_",tsk,"_md_",md.scen,".tsv",sep=""),
              sep="\t",quote=FALSE,col.names=FALSE,row.names=TRUE,na="?")
}

makepaup = function(MS_to_cluster,tsk,md.scen){
  modefunc <- function(x){
    tabresult <- tabulate(x)
    themode <- which(tabresult == max(tabresult))
    if(sum(tabresult == max(tabresult))>1) themode <- NA
    return(themode)
  }
  #produce nexus file for input to RevBayes
  modes = apply(MS_to_cluster, 2, modefunc)
  #replace NAs with modal value
  # MS_to_cluster = do.call("cbind",lapply(seq(length(modes)),
  #                                        function(x) {MS_to_cluster[which(is.na(MS_to_cluster[,x])),x]=modes[x]
  #                                        return(MS_to_cluster[,x])} ))
  MS_to_cluster2 = t(apply(MS_to_cluster,1,function(x) x-modes+1))
  minMS = min(MS_to_cluster2,na.rm = TRUE)
  MS_to_cluster2 = MS_to_cluster2 - minMS + 1
  #take only the MS which show a change
  logi = apply(MS_to_cluster2,2,function(x) length(unique(x[!is.na(x)]))>1)
  #change to letters
  if(sum(logi)>1){
    MS_to_cluster3 = matrix(LETTERS[MS_to_cluster2[,logi]],nrow=nrow(MS_to_cluster2[,logi]),
                            byrow=FALSE)
  }else{
    MS_to_cluster3 = matrix(LETTERS[MS_to_cluster2[,logi]],nrow=length(MS_to_cluster2[,logi]),
                            byrow=FALSE)
  }
  rownames(MS_to_cluster3) = paste("t",rownames(MS_to_cluster),sep="")
  MS_to_cluster3[is.na(MS_to_cluster3)] = "?"
  write.nexus.data(MS_to_cluster3,paste("ms_sim_data_paup_",tsk,"_md_",md.scen,".nex",sep=""),format="protein")
}


makerbcont = function(MS_to_cluster,tsk,md.scen){
  modefunc <- function(x){
    tabresult <- tabulate(x)
    themode <- which(tabresult == max(tabresult))
    if(sum(tabresult == max(tabresult))>1) themode <- NA
    return(themode)
  }
  #produce nexus file for input to RevBayes
  #assume mode is the original and centre around this value
  modes = apply(MS_to_cluster, 2, modefunc)
  MS_to_cluster2 = t(apply(MS_to_cluster,1,function(x) x-modes+1))
  minMS = min(MS_to_cluster2,na.rm = TRUE)
  MS_to_cluster2 = MS_to_cluster2 - minMS + 1
  #take only the MS which show a change
  write.table(max(MS_to_cluster2,na.rm=TRUE)+2,file=paste("nst_",tsk,"_md_",md.scen,".txt",sep=""),quote=FALSE,
              row.names = FALSE,col.names = FALSE)
  logi = apply(MS_to_cluster2,2,function(x) length(unique(x[!is.na(x)]))>1)
  write.table(MS_to_cluster2[,logi],paste("ms_sim_data_rb_",tsk,"_md_",md.scen,".nex",sep=""),
              sep="\t",quote=FALSE,col.names=FALSE,row.names=TRUE,na="?")
}
