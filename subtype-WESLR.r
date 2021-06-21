# ##############################functions###################################################
get.knn <- function(v, k){
  ind <- order(v, decreasing = T)
  return(ind[1:k])
}

get.knn.graph <- function(S, k){
  n <- nrow(S)
  S.knn <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    ind <- get.knn(S[i,], k)
    S.knn[i,ind] <- 1
    S.knn[ind,i] <- 1
  }
  return(S.knn)
}

subtypeWESLR <- function(Xs, Ls,Y_or,mclusters,mx, ml, mu,lam, gam1,gam2,sigma,delta,max.iter,eps){
  
  n <- ncol(Xs[[1]])
  c <- mclusters
  Gs <- vector("list", mx)
  e1ds <- vector("list", mx)
  ds <- rep(0, mx)
  for(i in 1:mx){
    ds[i] <- nrow(Xs[[i]])
    set.seed(1)
    Gs[[i]] <- matrix(runif(ds[i]*c), nrow = ds[i], ncol = c) 
    e1ds[[i]] <- rep(1, ds[i])
  }
  set.seed(2)
  F.mat=matrix(runif(n*c), nrow =n , ncol = c)
  check.step <- 20;
  Gs.old <- vector("list", mx)
  As <- vector("list", mx)
  As.pos <- vector("list", mx)
  As.neg <- vector("list", mx)
  Bs <- vector("list", mx)
  Bs.pos <- vector("list", mx)
  Bs.neg <- vector("list", mx) 
  
  
  Q=length(Y_or)
  LQ=vector("list", length(Y_or))
  for(i in 1:Q){
    LQ[[i]]<- diag(colSums(Y_or[[i]])) - Y_or[[i]]
  }
  
  alpha1=rep(1,mx)/mx
  alpha2=rep(1,ml)/ml
  L <- matrix(0, nrow = n, ncol = n)
  for(i in 1:mx){
    L <- L + alpha1[i]^gam1*Ls[[i]]
  }
  for(i in 1:ml){
    L <-L +delta*alpha2[i]^gam2* LQ[[i]]
  }
  
  t <- 0
  while(t < max.iter){
    t <- t + 1
    L.pos <- (abs(L) + L)/2
    L.neg <- (abs(L) - L)/2
    for(i in 1:mx){
      Gs.old[[i]] <- Gs[[i]]
      B <-  mu*t(Xs[[i]])%*%Gs[[i]]
    }
    B.pos <- (abs(B) + B)/2
    B.neg <- (abs(B) - B)/2
    F.mat=F.mat*sqrt((L.neg%*%F.mat+B.pos+sigma*F.mat)/(L.pos%*%F.mat+mx*mu*F.mat+B.neg+sigma*F.mat%*%t(F.mat)%*%F.mat))
    
    for(i in 1:mx){
      As[[i]] <- mu*Xs[[i]]%*%t(Xs[[i]]) + lam*(e1ds[[i]]%*%t(e1ds[[i]]))
      As.pos[[i]] <- (As[[i]] + abs(As[[i]]))/2
      As.neg[[i]] <- (abs(As[[i]]) - As[[i]])/2
      Bs[[i]] <- matrix(0,nrow(Xs[[i]]),c)
      Bs[[i]] <-mu*Xs[[i]]%*%F.mat
      Bs.pos[[i]] <- (Bs[[i]] + abs(Bs[[i]]))/2
      Bs.neg[[i]] <- (abs(Bs[[i]]) - Bs[[i]])/2
    }
    
    for(i in 1:mx){
      Gs[[i]] <- Gs[[i]]*sqrt((Bs.pos[[i]] + As.neg[[i]]%*%Gs[[i]])/(Bs.neg[[i]] + As.pos[[i]]%*%Gs[[i]]))
    }
    
    for(i in 1:mx){
      alpha1[i] <- (1/sum(diag(t(F.mat)%*%Ls[[i]]%*%F.mat)))^(1/(gam1 - 1))
    }
    alpha1 <- alpha1/sum(alpha1)
    for(i in 1:ml){
      alpha2[i] <- (1/sum(diag(t(F.mat)%*%LQ[[i]]%*%F.mat)))^(1/(gam2 - 1))
    }
    alpha2 <- alpha2/sum(alpha2)
    
    diff.G <- rep(0, mx)
    for(i in 1:mx){
      diff.G[i] <- norm(Gs[[i]] - Gs.old[[i]], "f")/norm(Gs.old[[i]], "f")
    }
    
    if(t%%check.step == 0){
      mesg <- sprintf("t = %i, diffG mean = %.4e", t, mean(diff.G))
      print(mesg)
    }
    if(mean(diff.G) < eps)
      break
  }
  return(list(Gs = Gs, F.mat = F.mat, alpha1 = alpha1,alpha2 = alpha2, diff.G = diff.G, t = t))
}

###################################  subtype-WESLR  ###########################################
## load the synthetic data and cluster indexs of base clustering
load("./subtype-WESLR.RData")


## initialize the hyperparameters
mx <- 3  #the type of omic data   #
K=20
theta =0.5
mclusters=4 #
eps=1e-5
max.iter=2000
mu=0.001
lambda =1
gamma1 =2
gamma2 = 2
sigma=100
delta=0.1
ml=length(clusters) #  kinds of base clustering
n=nrow(dataL[[1]])
truelabel = c(matrix(1,50,1),matrix(2,50,1),matrix(3,50,1),matrix(4,50,1));
Y_or =list(matrix(0,n,n),matrix(0,n,n),matrix(0,n,n),matrix(0,n,n))

## compute the graph Laplacian matrix of each data type 
library(SNFtool)
distL = lapply(dataL, function(x) (dist2(x, x))^(1/2))
affinityL = lapply(distL, function(x) affinityMatrix(x, K, theta))
S1.knn <- get.knn.graph(affinityL[[1]], K)
S2.knn <- get.knn.graph(affinityL[[2]], K)
S3.knn <- get.knn.graph(affinityL[[3]], K)
D1 <- diag(colSums(S1.knn))
D2 <- diag(colSums(S2.knn))
D3 <- diag(colSums(S3.knn))
L1 <- D1 - S1.knn
L2 <- D2 - S2.knn
L3 <- D3 - S3.knn


## compute the indicator matrix of each base clustering
for(i in 1:length(clusters)){
  for(j in 1:n){
    for(k in 1:n){
      if(clusters[[i]][j]==clusters[[i]][k]){
        Y_or[[i]][j,k]=1
        Y_or[[i]][k,j]=1
      }
    }
  }
}
## run subtype-WESLR
Xs=list(t(dataL[[1]]),t(dataL[[2]]),t(dataL[[3]]))
Ls<-list(L1, L2, L3)
train.res <-subtypeWESLR(Xs, Ls,Y_or,mclusters,mx, ml, mu,lambda,gamma1,gamma2,sigma,delta,max.iter,eps)
clusters_subtype-WESLR=max.col(train.res$F.mat)
NMI_subtype-WESLR = calNMI(clusters_subtype-WESLR,truelabel);
NMI_moCluster = calNMI( clusters[[1]],truelabel);
NMI_iClusterPlus = calNMI(clusters[[2]],truelabel);
NMI_LRAcluster= calNMI(clusters[[3]],truelabel);
NMI_SNF = calNMI(clusters[[4]],truelabel);
