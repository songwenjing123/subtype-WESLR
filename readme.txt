## subtype-WESLR

#subtype-WESLR.r is the R file of the subtype-WESLR method.
#There are three functions including get.knn, get.knn.graph,subtypeWESLR.
#A graph model can be conscructed by get.knn and get.knn.graph.
#subtypeWESLR is the main function.
#Synthetic data of three data types and base clustering of moCluster,iClusterPlus,LRAcluster, and SNF
# are contained in subtype-WESLR.RData.



# #Firstly, subtype-WESLR.RData needs to be loaded.
# load("./subtype-WESLR.RData")



##Secondly, initialize the hyperparameters.
##For example
# mx <- 3  #the type of omic data   #
# K=20
# theta =0.5
# mclusters=4 #
# eps=1e-5
# max.iter=2000
# mu=0.001
# lambda =1
# gamma1 =2
# gamma2 = 2
# sigma=100
# delta=0.1


##Thirdly, compute the graph Laplacian matrix of each data type and the indicator matrix of each base clustering.


##Finally,run subtype-WESLR
