## subtype-WESLR

subtype-WESLR.r is the R file of the subtype-WESLR method. There are three functions including get.knn, get.knn.graph,subtypeWESLR. A graph model can be conscructed by get.knn and get.knn.graph. subtypeWESLR is the main function. Synthetic data of three data types and base clustering of moCluster,iClusterPlus,LRAcluster, and SNF are contained in subtype-WESLR.RData. There are four steps running subtype-WELSR in the demo.

### 1. Load subtype-WESLR.RData.

```
load("./subtype-WESLR.RData")
```

### 2. Initialize the hyperparameters.

For example

```
mx <- 3  #the type of omic data
K=20
theta =0.5
mclusters=4
eps=1e-5
max.iter=2000
mu=0.001
lambda =1
gamma1 =2
gamma2 = 2
sigma=100
delta=0.1
```

### 3. Compute the graph Laplacian matrix of each data type and the indicator matrix of each base clustering.

```
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
```

```
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
```

### 4. Run subtype-WESLR.

```
Xs=list(t(dataL[[1]]),t(dataL[[2]]),t(dataL[[3]]))
Ls<-list(L1, L2, L3)
train.res <-subtypeWESLR(Xs, Ls,Y_or,mclusters,mx, ml, mu,lambda,gamma1,gamma2,sigma,delta,max.iter,eps)
```

