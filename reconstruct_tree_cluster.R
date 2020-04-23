dist_recons = function(data,ini_len,md.scen,tsk,i,j){
  trees = list()
  Mode <- function(x) {
    ux <- na.omit(unique(x) )
    tab <- tabulate(match(x, ux)); ux[tab == max(tab) ][1]
  }
  distimputemean <- function(x) {
    #set NA to mode value for that col
    means=apply(x,2,function(x) mean(x,na.rm=TRUE))
    for(i in 1:ncol(x)){
      x[is.na(x[,i]),i] = means[i]
    }
    di = dist(x,method="euclidean")
    return(di)
  }
  distimputemode <- function(x) {
    #set NA to mode value for that col
    modes=apply(x,2,function(x) Mode(x))
    for(i in 1:ncol(x)){
      x[is.na(x[,i]),i] = modes[i]
    }
    di = dist(x,method="euclidean")
    return(di)
  }
  distuseobs1 <- function(x){
    d1 = matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(i in 1:(nrow(x)-1)){
      for(j in (i+1):nrow(x)){
        n.both.obs = sum(apply(!is.na(rbind(x[i,],x[j,])),2,all))
        if(n.both.obs>0){
          d1[i,j] = sum(abs(x[i,]-x[j,]),na.rm=TRUE)
          if(is.na(d1[i,j])){d1[i,j]=0}
          d1[j,i] = d1[i,j]
        }else{
          d1[i,j] = 0
          d1[j,i] = 0
        }
      }
    }
    rownames(d1) = rownames(x)
    di = as.dist(d1)
    return(di)
  }
  distuseobs2 <- function(x){
    d1 = matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(i in 1:(nrow(x)-1)){
      for(j in (i+1):nrow(x)){
        n.both.obs = sum(apply(!is.na(rbind(x[i,],x[j,])),2,all))
        if(n.both.obs>0){
          d1[i,j] = sqrt(sum((x[i,]-x[j,])^2,na.rm=TRUE))
          d1[j,i] = d1[i,j]
        }else{
         d1[i,j] = 0
         d1[j,i] = 0
        }
      }
    }
    rownames(d1) = rownames(x)
    di = as.dist(d1)
    return(di)
  }
  distuseobs1scale <- function(x){
    d1 = matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(i in 1:(nrow(x)-1)){
      for(j in (i+1):nrow(x)){
        n.both.obs = sum(apply(!is.na(rbind(x[i,],x[j,])),2,all))
        if(n.both.obs>0){
          d1[i,j] = ncol(x)*sum(abs(x[i,]-x[j,]),na.rm=TRUE)/n.both.obs
          d1[j,i] = d1[i,j]
        }else{
          d1[i,j]=0
          d1[j,i]=0
        }
      }
    }
    rownames(d1) = rownames(x)
    di = as.dist(d1)
    return(di)
  }
  distuseobs2scale <- function(x){
    d1 = matrix(0,nrow=nrow(x),ncol=nrow(x))
    for(i in 1:(nrow(x)-1)){
      for(j in (i+1):nrow(x)){
        n.both.obs = sum(apply(!is.na(rbind(x[i,],x[j,])),2,all))
        if(n.both.obs>0){
          d1[i,j] = sqrt(sum((ncol(x)*(x[i,]-x[j,])/n.both.obs)^2,na.rm=TRUE))
          d1[j,i] = d1[i,j]
        }else{
          d1[i,j]=0
          d1[j,i] = d1[i,j]
        }
      }
    }
    rownames(d1) = rownames(x)
    di = as.dist(d1)
    return(di)
  }
  if(md.scen>1){
    dm1 = distimputemean(data)
    dm2 = distimputemode(data)
    dm3 = distuseobs1(data)
    dm4 = distuseobs2(data)
    dm5 = distuseobs1scale(data)
    dm6 = distuseobs2scale(data)
    trees$treeNJeu1 = NJ(dm1)
    write.tree(trees$treeNJeu1,paste("nj_tree1_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeBal1 = fastme.bal(dm1, nni = TRUE, spr = TRUE, tbr = TRUE)
    write.tree(trees$treeBal1,paste("bal_tree1_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeNJeu2 = NJ(dm2)
    write.tree(trees$treeNJeu2,paste("nj_tree2_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeBal2 = fastme.bal(dm2, nni = TRUE, spr = TRUE, tbr = TRUE)
    write.tree(trees$treeBal2,paste("bal_tree2_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeNJeu3 = NJ(dm3)
    write.tree(trees$treeNJeu3,paste("nj_tree3_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeBal3 = fastme.bal(dm3, nni = TRUE, spr = TRUE, tbr = TRUE)
    write.tree(trees$treeBal3,paste("bal_tree3_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeNJeu4 = NJ(dm4)
    write.tree(trees$treeNJeu4,paste("nj_tree4_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeBal4 = fastme.bal(dm4, nni = TRUE, spr = TRUE, tbr = TRUE)
    write.tree(trees$treeBal4,paste("bal_tree4_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeNJeu5 = NJ(dm5)
    write.tree(trees$treeNJeu5,paste("nj_tree5_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeBal5 = fastme.bal(dm5, nni = TRUE, spr = TRUE, tbr = TRUE)
    write.tree(trees$treeBal5,paste("bal_tree5_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeNJeu6 = NJ(dm6)
    write.tree(trees$treeNJeu6,paste("nj_tree6_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeBal6 = fastme.bal(dm6, nni = TRUE, spr = TRUE, tbr = TRUE)
    write.tree(trees$treeBal6,paste("bal_tree6_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    data[is.na(data)] <- 0
    dmat = matrix(0,nrow=nrow(data),ncol=nrow(data) )
    for(l in 1:(nrow(data)-1)){
      for(k in (l+1):nrow(data)){
        dmat[l,k] = 1 - sum(data[l,]*data[k,])/sqrt(sum(data[l,]^2)*sum(data[k,]^2))
      }
    }
    dmat[lower.tri(dmat)] <- t(dmat)[lower.tri(dmat)]
    rownames(dmat) = rownames(data)
    colnames(dmat) = rownames(data)
    dm = as.dist(dmat)
    ptm <- proc.time()
    trees$treeNJcos = NJ(dm)
    t.NJ = (proc.time() - ptm)[3]
    write.tree(trees$treeNJcos,paste("njcos_tree_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    ptm <- proc.time()
    trees$treeBalcos = fastme.bal(dm, nni = TRUE, spr = TRUE, tbr = TRUE)
    t.bal = (proc.time() - ptm)[3]
    write.tree(trees$treeBalcos,paste("balcos_tree_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$t.NJ = t.NJ
    trees$t.bal = t.bal
  }else{
    dm.ma = dist(data,method = "manhattan")
    dm.eu = dist(data,method = "euclidean")
    trees$treeNJeu = NJ(dm.eu)
    write.tree(trees$treeNJeu,paste("njeu_tree_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeBaleu = fastme.bal(dm.eu, nni = TRUE, spr = TRUE, tbr = TRUE)
    write.tree(trees$treeBaleu,paste("baleu_tree_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeNJma = NJ(dm.ma)
    write.tree(trees$treeNJma,paste("njma_tree_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$treeBalma = fastme.bal(dm.ma, nni = TRUE, spr = TRUE, tbr = TRUE)
    write.tree(trees$treeBalma,paste("balma_tree_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    dmat = matrix(0,nrow=nrow(data),ncol=nrow(data) )
    for(l in 1:(nrow(data)-1)){
      for(k in (l+1):nrow(data)){
        dmat[l,k] = 1 - sum(data[l,]*data[k,])/sqrt(sum(data[l,]^2)*sum(data[k,]^2))
      }
    }
    dmat[lower.tri(dmat)] <- t(dmat)[lower.tri(dmat)]
    rownames(dmat) = rownames(data)
    colnames(dmat) = rownames(data)
    dm = as.dist(dmat)
    ptm <- proc.time()
    trees$treeNJcos = NJ(dm)
    t.NJ = (proc.time() - ptm)[3]
    write.tree(trees$treeNJcos,paste("njcos_tree_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    ptm <- proc.time()
    trees$treeBalcos = fastme.bal(dm, nni = TRUE, spr = TRUE, tbr = TRUE)
    t.bal = (proc.time() - ptm)[3]
    write.tree(trees$treeBalcos,paste("balcos_tree_",tsk,"_",i,"_",j,"_",md.scen,".new",sep=""))
    trees$t.NJ = t.NJ
    trees$t.bal = t.bal
  }
  return(trees)
}
