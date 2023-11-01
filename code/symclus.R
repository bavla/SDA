# clustering mixed symbolic data
# 30. October 2023
# Vladimir Batagelj

dSkip <- function(Y,p,q){return(0)}

dMembers <- function(Y,p,q){
  P <- Y[[p]]$R; Q <- Y[[q]]$R
  R <- P+Q; M <- max(R); t <- as.integer(2*R >= M)
  return(sum((1-t)*R + t*(M-R)))
}

dIntSq <- function(Y,p,q){
  P <- Y[[p]]$L; Q <- Y[[q]]$L
  wp <- Y[[p]]$s; wq <- Y[[q]]$s
  return(wp*wq*sum((P-Q)**2)/(wp+wq))
}

# computes dissimilarity between SOs
# global: alpha, dSel, nSel
distSO <- function(U,p,q){  
  D <- numeric(nSel)
  for(i in 1:nSel) {
    X <- dSel[[i]]; d <- X$d; Y <- U[[i]]
    D[i] <- d(Y,p,q) 
  }
  dis <- as.numeric(D %*% alpha)
  if (is.na(dis)) dis <- Inf
  return(dis)
}

updateL <- function(U,dSel,j,ip,iq){
  dt <- dSel[[j]]$dType; Y <- U[[j]]
  if(dt == "skip"){
    return(list(L="",R="",s=0))
  } else if(dt == "membersR"){
    P <- Y[[ip]]$R; Q <- Y[[iq]]$R
    R <- P+Q; M <- max(R); t <- as.integer(2*R >= M)
    s <- Y[[ip]]$s + Y[[iq]]$s
    return(list(L=t,R=R,s=Y[[ip]]$s+Y[[iq]]$s))
  } else if(dt == "intervalSq"){
    P <- Y[[ip]]$L; Q <- Y[[iq]]$L
    Pr <- Y[[ip]]$R; Qr <- Y[[iq]]$R
    wp <- Y[[ip]]$s; wq <- Y[[iq]]$s
    t <- (wp*P+wq*Q)/(wp+wq); R <- c(min(Pr,Qr),max(Pr,Qr))
    return(list(L=t,R=R,s=wp+wq))
  } else cat(j,ip,iq, "Error\n")
}

hclustSO <- function(SD,dSel){
  orDendro <- function(i){if(i<0) return(-i)
    return(c(orDendro(m[i,1]),orDendro(m[i,2])))}

  nUnits <- SD$head$nUnits; nmUnits <- nUnits-1; nSel <- length(dSel)
  npUnits <- nUnits+1; n2mUnits <- nUnits+nmUnits
  w <- rep(1,nUnits)
  alpha <<- vars <- rep(NA,nSel)
  for(i in 1:nSel) {X <- dSel[[i]]; vars[i] <- X$var; alpha[i] <<- X$alpha }
  H <- SD$SDF[,vars]; U <- H
  for(i in 1:nSel) for(j in 1:nUnits) U[[i]][[j]] <- list(L=H[[i]][[j]],R=H[[i]][[j]],s=1)
  D <- matrix(nrow=nUnits,ncol=nUnits)
  for(p in 1:nmUnits) for(q in (p+1):nUnits) { 
    D[q,p] <- D[p,q] <- distSO(U,p,q)
  }
  diag(D) <- Inf
  active <- 1:nUnits; m <- matrix(nrow=nmUnits,ncol=2)
  node <- rep(0,nUnits); h <- numeric(nmUnits)
  for(j in npUnits:n2mUnits) { U[nrow(U)+1,] <- vector("list",nSel)
    for(i in 1:nSel) U[[i]][[j]] <- list(L=NA,R=NA,s=NA)}
  rownames(U)[npUnits:n2mUnits] <- paste("L",1:nmUnits,sep="")
  for(k in 1:nmUnits){
    ind <- active[sapply(active,function(i) which.min(D[i,active]))]
    dd <- sapply(active,function(i) min(D[i,active]))
    pq <- which.min(dd)
    p<-active[pq]; q <- ind[pq]; h[k] <- D[p,q]
    if(node[p]==0){m[k,1] <- -p; ip <- p
    } else {m[k,1] <- node[p]; ip <- node[p]}
    if(node[q]==0){m[k,2] <- -q; iq <- q
    } else {m[k,2] <- node[q]; iq <- node[q]}
    ik <- nUnits + k
    for(j in 1:nSel) U[[j]][[ik]] <- updateL(U,dSel,j,ip,iq)
    active <- setdiff(active,p)
    for(s in setdiff(active,q)){
      is <- ifelse(node[s]==0,s,node[s])
      D[s,q] <- D[q,s] <- distSO(U,ik,is)
    }
    node[[q]] <- ik
  }
  for(i in 1:nmUnits) for(j in 1:2) if(m[i,j]>nUnits) m[i,j] <- m[i,j]-nUnits 
  hc <- list(merge=m,height=h,order=orDendro(nmUnits),labels=rownames(SD$SDF),
    method=NULL,call=NULL,dist.method=NULL,leaders=U[npUnits:n2mUnits,])
  class(hc) <- "hclust"
  return(hc)
}

