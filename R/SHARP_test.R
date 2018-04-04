# rm(list=ls())#remove all 
# library(gplots)
# #library(ggplot2)
# library(reshape2)
# library(pvclust)
# library(RColorBrewer);
# library(stringr)
# library(clues)
# library(Matrix)
# library(Rtsne)
# library (cluster)
#This is to implement SHARP (Single-cell RNA-seq data Hyper-fast and Accurate via ensemble Random Projection) on different parameters (including ensemble size, reduced dimensions and different runs)

#by Shibiao WAN, Feb-27-2018, Philadelphia, PA, USA


#sub-function to be called


#' Obtain a random matrix
#'
#' This function is to generate a random matrix according to a sparse-version Achlioptas distribution (see our paper for details), given the input matrix and the dimension we want to reduce to.
#'
#' @param scdata input single-cell expression matrix
#' @param p the dimension to be reduced to
#' 
#' @examples
#' newE = RPmat(E, p)
#'
#' @import Matrix
#'
#' @export
RPmat <- function(scdata, p){#for random projection; note scdata: m*n, m is the feature dimensions, n is the sample number; p is the reduced dimension
  m = nrow(scdata)#the number of features
  n = ncol(scdata)#number of samples/cells
#   set.seed(123)#a flag to fix the random number
  s = sqrt(m)#according to the paper "Very Sparse Random Projection"
  x0 = sample(c(sqrt(s),0, -sqrt(s)), size= m*p, replace=TRUE, prob=c(1/(2*s), 1 - 1/s, 1/(2*s)))
#   return(x)
  x <- Matrix(x0, nrow = m, byrow = TRUE, sparse = TRUE)#reshape a vector to a sparse matrix
  
  scmat = Matrix(data.matrix(scdata), sparse = TRUE)#convert it to a sparse matrix
  projmat = 1/sqrt(p)*t(x)%*%scmat#the projected matrix by random projection;matrix multiplication %*%
#   print("It works!")
#   rmat = ranmat(m, p)#the random matrix

#   if(n <10000){
#     x <- matrix(x0, nrow = m, byrow = TRUE)#reshape a vector to a matrix
# #     scmat = data.matrix(scdata)#convert the data frame to a matrix
#     projmat = 1/sqrt(p)*t(x)%*%scdata#the projected matrix by random projection;matrix multiplication %*%
#   }else{
#     x2 <- Matrix(x0, nrow = m, byrow = TRUE, sparse = TRUE)#reshape a vector to a matrix
#     projmat = 1/sqrt(p)*t(x2)%*%scdata#the projected matrix by random projection;matrix multiplication %*%
#   }

#   projmat = t(projmat)
  list1 = list()
  list1$R = t(x)
  list1$projmat = projmat
  return(list1)#the same format, feature*sample
}

# getbA <- function(){
#   r = which(rowColor %in% k)#for the data in the same cluster with value k
#   ind = expand.grid(r,r)#all the potential combinations
#   tmp = sparseMatrix(i = ind[,1], j = ind[,2], x= 1, dims = c(N,N))#non-zero entries
# }

#' Obtain the co-location matrix
#'
#' This function is to obtain the weighted co-location matrix from the clustering result for following ensemble clustering.
#'
#' @param rowColor the clustering results by some clustering algorithms
#' 
#' @examples
#' AA = getA(rowColor)
#'
#' @import Matrix
#'
#' @export
getA <- function(rowColor){
#This is to obtain the weighted co-association matrix for clustering solution rowColor
  N = length(rowColor)#number of points

#   -------------------------------------------------------
#   #convert a (categorical) vector to a numeric vector 
#   t1 = factor(rowColor)#rowColor is a vector like "red red purple brown ..."
#   levels(t1) = c(1:length(levels(t1)))#convert the categorical to numeric
#   t1 = as.numeric(as.character(t1))
#   -------------------------------------------------------

# t1 = factor(rowColor)#rowColor is a vector like "red red purple brown ..."
# levels(t1) = c(1:length(levels(t1)))#convert the categorical to numeric
# t1 = as.numeric(as.character(t1))

  L = levels(factor(rowColor))
#   A = matrix(0, N, N)#declare a 0-value matrix
#   A = Matrix(0, nrow = N, ncol = N, sparse = T)#use a parse matrix for efficient storage
# ind = which(rowColor[, c] %in% k)#the same cluster

  #find indices for each cluster, then all combinations of indices
  tmp = sapply(L, function(k){r = which(rowColor %in% k);expand.grid(r,r)})
  #reshape to the indices
  allind = matrix(unlist(t(tmp)), ncol = 2, byrow= F)#need transpose
  A = sparseMatrix(i = allind[,1], j = allind[,2], x= 1, dims = c(N,N))#non-zero entries
  
#   for (k in L) {
#     r = which(rowColor %in% k)#for the data in the same cluster with value k
#     ind = expand.grid(r,r)#all the potential combinations
#     tmp = sparseMatrix(i = ind[,1], j = ind[,2], x= 1, dims = c(N,N))#non-zero entries
#     A = A + tmp
# #     A[which(rowColor %in% k), which(rowColor %in% k)] = 1
#   }
  
#   r = which(rowColor %in% k); ind = expand.grid(r,r)
#   apply(L, function(k) A[which(rowColor %in% k), which(rowColor %in% k)] = 1)
  return(A)
}

getnewk <- function(k, R, x, N){#This is to get the original index of the sample 
  k1 = which(x %in% R[k])#find samples with k-th cluster
  d1 = unlist(strsplit(R[k], "_"))#the name contains only two parts; get the numbering part
  d = as.numeric(tail(d1, n = 1))#the last element of the split arrays
  newk1 = k1 - (d - 1)*N#the index
  return(newk1)
}

#' Obtain the co-location matrix
#'
#' This function is to obtain the weighted co-location matrix from the clustering result for following ensemble clustering.
#'
#' @param rowColor the clustering results by some clustering algorithms
#' 
#' @examples
#'  x = as.vector(sapply(1:C, function(i){paste(nC[,i], "_", i, sep = "")}))#convert the matrix (N*C) to vector (concatenating them)
#'  R = unique(x)#all unique labels
#'  alls = apply(cb, 2, getss, R = R, x = x, w1 = w1)#calculate the weight s for all combinations 
#'
#' @export
getss<- function(pind, R, x, w1){#This is to get the element of S
  pairk = lapply(pind, getnewk, R = R, x = x, N = length(w1))#run for two indices
  
  intset = intersect(unlist(pairk[1]), unlist(pairk[2]))#set intersection
  
  ss = 0
  if (length(intset) != 0){
    uset = union(unlist(pairk[1]),unlist(pairk[2]))#set union
    ss = sum(w1[intset])/sum(w1[uset])
  }
  return(ss)
}

#' weighted ensemble clustering
#'
#' This function is to do weighted ensemble clustering for meta-clustering.
#'
#' @param nC a m*n matrix, where m is the number of cells, n is the number of clustering algorithms, and the (i,j)-element of nC represents the cluter for the i-th cell by the j-th clutering predictor.
#' 
#' @examples
#' finalrowColor = wMetaC(enrp)
#'
#' @import cluster
#'
#' @export
wMetaC <- function(nC){
  #This is to obtain the weight matrix for each cluster solution for following meta-clustering
  N = nrow(nC)#number of points
  C = ncol(nC)#number of clustering methods/times; or K
#   AA = Matrix(0, nrow = N, ncol = N, sparse = T)#initialization
#   for (i in 1:C){
#     AA = AA + getA(nC[,i])
#   }
  
  AA = Reduce('+', apply(nC, 2, getA))#sum over obtained matrices; execute along the column and then matrix sum
  AA = AA/C
#   AA = sum(apply(nC, 2, getA))/C
  indA = which(AA!=0,arr.ind = T)#find non-zero indices of AA
  nd = vapply(AA[indA], function(x) x*(1-x), numeric(1))
  newAA = sparseMatrix(i = indA[,1], j = indA[,2], x= nd, dims = c(N,N))
#   AA[] <- vapply(AA, function(x) x*(1-x), numeric(1))#apply to every element
  w0 = 4/N*rowSums(newAA)#the weight for each point
#   e = 0.001
  e = 0.01
  w1 = (w0 + e)/(1 + e)#adjusted point weight
  
#   #revise the labels for differentiation
#   x <- vector(mode="character", length=0)
#   for (i in 1:C){
# #     t1 = factor(nC[,i])#rowColor is a vector like "red red purple brown ..."
# #     levels(t1) = paste(levels(t1), "_", i, sep = "")#make each cluster in each column different
#     t1 = paste(nC[,i], "_", i, sep = "")
#     x = c(x, t1)#reshape the matrix to a vector
#   }
#   
#   #Cluster weight value
#   
#   newnC <- matrix(x, nrow = N, byrow = FALSE)#reshape a vector to a matrix; by column
  
  #revise the labels for differentiation
# #   newnC = matrix('0', nrow = N, ncol = C)#revised cluster name
#   x = character(length = 0)#initialization
#   for (i in 1:C){
#       newnC = paste(nC[,i], "_", i, sep = "")#distinguish the cluster labels for each column/solution
#       x = c(x, newnC)#extend the matrix to a vector for ease of processing
# #     newnC[,i] = paste(nC[,i], "_", i, sep = "")#distinguish the cluster labels for each column/solution
# #     x = c(x, newnC[,i])#extend the matrix to a vector for ease of processing
#   }
  x = as.vector(sapply(1:C, function(i){paste(nC[,i], "_", i, sep = "")}))#convert the matrix (N*C) to vector (concatenating them)
  newnC <- matrix(x, nrow = N, byrow = FALSE)#reshape a vector to a matrix; by column
  
  R = unique(x)#all unique labels
  allC = length(R)#number of all unique labels
#   S = matrix(0, allC, allC)
  
  cb = combn(allC,2)#all possible combinations (n*(n-1)/2)
  alls = apply(cb, 2, getss, R = R, x = x, w1 = w1)#calculate the weight s for all combinations
  S0 = sparseMatrix(i = cb[1,], j = cb[2,], x= alls, dims = c(allC,allC))#triangle part of the S
  S = S0 + t(S0) + diag(allC)
#   for (k in c(1:(allC-1))){
#     k1 = which(x %in% R[k])#find samples with k-th cluster
#     d1 = as.numeric(unlist(strsplit(R[k], "_"))[2])#the name contains only two parts; get the numbering part
#     #convert d1 to an integer
#     newk1 = k1 - (d1 - 1)*N#the index
#     for (j in c((k + 1):allC)){
#       k2 = which(x %in% R[j]) #%% N#modulus over N
#       d2 = as.numeric(unlist(strsplit(R[j], "_"))[2])#the name contains only two parts; get the numbering part
#       newk2 = k2 - (d2 - 1)*N#the index
#       
#       intset = intersect(newk1, newk2)#set intersection
#       uset = union(newk1,newk2)#set union
#       if (length(intset) != 0){#if 0, no need to calculate
# 	S[k, j] = sum(w1[intset])/sum(w1[uset])
# 	S[j, k] = S[k, j]#symmetric
#       }
#      }
#      S[k, k] = 1
#   }


#       if ((length(k1) != 0) && (length(k2) != 0)){
# # 	if (k1 == 0){
# # 	  k1 = N
# # 	}
# 	k1[k1 == 0] = N#the final one
# 	k2[k2 == 0] = N
#       
# # 	if (k2 == 0){
# # 	  k2 = N
# # 	}
# 	intset = intersect(k1, k2)
# 	uset = union(k1,k2)
# 	if (length(intset) != 0){
# 	  S[k, j] = sum(w1[intset])/sum(w1[uset])
# # 	  S[j, k] = S[k, j]#symmetric
# 	}
#       }
#     }
    
  d=as.dist(1-S)
#   metah = hclust(d, method="ward.D")
  metah = hclust(d, method="ward.D")
  
  ###################determine the optimal number of clusters###################
  maxc = min(40, dim(S)[1]-1)#maximum number of clusters; at most (maxc-1) clusters because maxc clusters mean no needing to cluster
  v = matrix(0, nrow = dim(S)[1], ncol = maxc-1)#for all different numbers of clusters
#   print(dim(v))
  msil = rep(0, maxc-1)#declare a vector of zeros
  
#   mdunn = rep(0, maxc-1)#for dunn index
#   mdb = rep(0, maxc-1)#for DB index
#   wss = rep(0, maxc-2)#within-cluster sum of squares
#   CHind = rep(0, maxc-1)#CH index
  

  allc = c(2:maxc)
  for(i in 1:length(allc)){
    v[, i] = cutree(metah, k = allc[i])#for different numbers of clusters
#     print(v[,i])
    sil = silhouette(v[, i], d)#calculate the silhouette index 
    
#     msil[i-1] = mean(sil[,3])#the mean value of the index
    msil[i] = median(sil[,3])#the mean value of the index
    
#     mdunn[i-1] = dunn(d, v[,i-1])
    
#     db = index.DB(d, cl = v[, i-1])
#     mdb[i-1] = db$DB

#     spl <- split(d, v[, i-2])
#     wss[i-2] <- sum(sapply(spl, wss))#within-cluster sum of squares
#     print("OK")
#     CHind[i] = 
	  (S, v[, i], disMethod = "1-corr")
#     print(CHind[i])
  }

#   oind = which.max(msil)#the index corresponding to the max sil index
  tmp = which(msil == max(msil))#in case there are more than one maximum
  if(length(tmp)>1){oind = tmp[ceiling(length(tmp)/2)]}else{oind = tmp}
#   print(msil)#the average sil index for different numbers of clusters
  
#   print("OK!")
#   print(CHind)
  
#   if(max(msil)<=0.2){
#     difCH = diff(CHind)
#     x1 = which(difCH<0)
#     if(length(x1)>0){
#       oind = min(x1)#local maximum
#     }else{
#       oind = 1
#     }
# #   oind = which.max(CH)
#   }

#   #if the silhouette index is not reliable, we use the gap statistics
#   if(max(msil)<=0.1){
#     
#     
#     gskmn <- clusGap(S, FUN = kmeans, nstart = 20, K.max = min(maxc, 15), B = min(50, nrow(S)))
#     g = gskmn$Tab
#     gap = g[, "gap"]#the gap info
#     print(gap)
#     
# #     oind = which.max(gap)#maximum gap
#     
#     sesim = g[, "SE.sim"]#standard error of the gap
#     print(sesim)
#     oind = maxSE(gap, sesim)#maximize the gap with parsimony of the model
#     
# #     if(oind >=floor(maxc*0.8)){#if the gap stastic keeps increasing until very big number of cluster, we use the first SE-rule (Gap(k) >= Gap(k+1) - SE)
# #       sesim = g[, "SE.sim"]#standard error of the gap
# #       print(sesim)
# #       oind = maxSE(gap, sesim)#maximize the gap with parsimony of the model
# #     }
#   }

#   print(wss)
#   print(mdunn)
#   print(mdb)
#   oind = 10
  tf = v[,oind]#the optimal clustering results
  ############################################################################
  
#   tf=cutree(metah,k=N.cluster)
  newnC[] <- vapply(newnC, function(x) tf[match(x, R)], numeric(1))#apply to every element
  finalC = apply(newnC, 1, function(d) names(sort(table(d),decreasing=TRUE)[1])) #find the most repeated elements for each row
  
  N.cluster = length(unique(finalC))#note that the number of clusters for meta-clustering is not determined by previous selection, but by the unique number in the final round.
  
  print(paste("The optimal number of clusters for ensemble clustering is: ", N.cluster, sep = ""), quote = FALSE)
  
  return(finalC)
}

#' hierarchical clustering with number of clusters determined automatically
#'
#' This function is to do hierarchical clustering, with the number of clusters determined by a strategy combining Silhouette index, CH index and the heights after hierarchical clustering.
#'
#' @param E the feature matrix whose columns represent the features and whose rows represent data/cells.
#' @param tag the tag for each independent random projection
#' @param colorL the color set for different clusters
#'
#' @examples
#' rowColor= getrowColor(E1, tag, colorL)
#'
#' @import cluster
#'
#' @import clues
#'
#' @export
# getrowColor <- function(E, tag, outdir, colorL){#hierarchical clustering
getrowColor <- function(E, tag, colorL){#hierarchical clustering
	
	
	repnum = ncol(E)

	#mExp = E[,c("WT.0h_1","WT.0h_2","WT.8h_1","WT.8h_2","KO.0h_1","KO.0h_2","KO.8h_1","KO.8h_2")]
# 	mExp= E[,c(3,4,1,2,7,8,5,6)]

# 	mExp= E[,c(1:40)]#the number of replicates
	mExp= E[,c(1:repnum)]#change this one######
	#For each experiment, note to change the number of replicates########
	
	#mlabel = E[,c("FC","FC")]#mExp = E[,]
	#mlabel[mlabel[,1]<0,1] <-  0
	#mlabel[mlabel[,2]>0,2] <-  0
	maxval = 2
	minval = -2
	my = mExp
	my = my[apply(my,1,function(x) sd(x)!=0),]#remove those rows whose elements are the same
	my <- t(scale(t(my)))

# 	cluster=vector(mode="character",length=nrow(mExp))

	d=as.dist(1-cor(t(my)))#cor function is to calculate the column-wise correlation
# 	d = dist(my, method = "euclidean")
	h=hclust(d, method="ward.D")#ward to ward.D
# 	h=hclust(d, method="ward.D")#ward to ward.D
# 	dend = as.dendrogram(h)
# 	lownum = 1126+1144+176 #(orange, purple, red)
# 	wGreen = lownum+1097
# 	myorder = c((wGreen+1):4998, (lownum+1):wGreen,1:lownum)
# 	myorder = c((wGreen+1):4998, 1:lownum,  (lownum+1):wGreen)
# 	dend <- reorder(dend,myorder,agglo.FUN=mean)
# 	#dend <- reorder(dend,1:4998)
# 
# 	d2 =as.dist(1-cor((my)))
# # 	h2 = hclust(d2,method="ward.D")
# 	h2=hclust(d2, method="ward.D")#ward to ward.D
# 	dend2 = as.dendrogram(h2)
# 	
# 	bk <- seq(-2, 2, by=0.1)
# 	data.mat=as.matrix(my)
  
#  	fname = unlist(strsplit(outdir, "/"))[-1]#part of the names
# 	pngnamedist = paste(outdir, fname, tag, ".dist.png",sep="")
# # 	unlink(pngnamedist)#remove the old figure
# 	png(pngnamedist,width=500,height=500)
# 	par(cex=1.2)
# 	plot(h2,main="distance",cex.main=1)
# 	dev.off()
# 
# 	pngname=paste(outdir,fname, tag,".Dendro.png",sep="")
# 	png(pngname,width=1000,height=1200)
# 	#par(mfrow=c(1,3))
# 	par(mar=c(10,1,1,1))

        
        ###################determine the optimal number of clusters###################
	nc = 2:40
	v = matrix(0, nrow = nrow(my), ncol = length(nc))#for all different numbers of clusters
	msil = rep(0, length(nc))#declare a vector of zeros
# 	wss = rep(0, length(nc))#within-cluster sum of squares
# 	mdunn = rep(0, 39)#for dunn index
# 	mdb = rep(0, 39)#for dunn index
	CHind = rep(0, length(nc))
	
# 	print(paste("The height for the top 10 are: ", tail(h$height, n = 10), sep = ""))
	
	my1 = as.matrix(my)#convert to full matrix
	my = my1
	for(i in 1:length(nc)){
	  v[,i] = cutree(h, k = nc[i])#for different numbers of clusters
	  sil = silhouette(v[,i], d)#calculate the silhouette index 
# 	  msil[i] = mean(sil[,3])#the mean value of the index
	  msil[i] = median(sil[,3])#the mean value of the index
	  
# 	  mdunn[i] = dunn(d, v[,i])
	  
# 	  db = index.DB(d, cl = v[, i])
# 	  mdb[i] = db$DB

# 	  #within-cluster sum of squares
# 	  spl <- split(d, v[,i])
# 	  wss[i] <- sum(sapply(spl, wss))
	  CHind[i] = get_CH(my, v[,i], disMethod = "1-corr")
	}

# 	print(msil)#the average sil index for all cases
# 	print(CHind)
# 	print(wss)
# 	print(mdunn)
# 	print(mdb)
# 	oind = which.max(msil)#the index corresponding to the max sil index
	tmp = which(msil == max(msil))#in case there are more than one maximum
	if(length(tmp)>1){oind = tmp[ceiling(length(tmp)/2)]}else{oind = tmp}
	
	if(max(msil)<=0.35){
	  oind = which.max(CHind)
	  if(oind ==1){#it's likely that the CH index is not reliable either
	    tmp = tail(h$height, n =10)#the height
	    diftmp = diff(tmp)
	    flag = diftmp > tmp[1:(length(tmp)-1)]#require the height is more than 2 times of the immediate consecutive one
	    
	    if(any(flag)){#if any satifies the condition; make sure at least one satisfied
	      pind = which.max(flag)
	      opth = (tmp[pind] + tmp[pind+1])/2#the optimal height to cut
	      optv = cutree(h, h = opth)#using the appropriate height to cut
	      oind = length(unique(optv)) - 1#for consistency
	    }
	  }
# 	  difCH = diff(CHind)
# 	  x1 = which(difCH<0)
# 	  if(length(x1)>0){
# 	    oind = min(x1)#local maximum
# 	  }else{
# 	    oind = 1
# 	  }
	  
	}
# 	#if the silhouette index is not reliable, we use the gap statistics
# 	if(max(msil)<=0.1){
# 	  gskmn <- clusGap(my, FUN = kmeans, nstart = 20, K.max = 15, B = min(nrow(my), 50))
# 	  g = gskmn$Tab
# 	  gap = g[, "gap"]#the gap info
# 	  print(gap)
#     
# # 	  oind = which.max(gap)#maximum gap
# 	  
# 	  sesim = g[, "SE.sim"]#standard error of the gap
# 	  print(sesim)
# 	  oind = maxSE(gap, sesim)#maximize the gap with parsimony of the model
#     
# # 	  if(oind >=floor(maxc*0.8)){#if the gap stastic keeps increasing until very big number of cluster, we use the first SE-rule (Gap(k) >= Gap(k+1) - SE)
# # 	    sesim = g[, "SE.sim"]#standard error of the gap
# # 	    print(sesim)
# # 	    oind = maxSE(gap, sesim)#maximize the gap with parsimony of the model
# # 	  }
# 	}
	
# 	oind = 5
	f = v[,oind]#the optimal clustering results
	N.cluster = nc[oind]#the optimal number of clusters
# 	N.cluster = oind + 1#the optimal number of clusters
	print(paste("The optimal number of clusters for individual RP is: ", N.cluster, sep = ""), quote = FALSE)
	############################################################################
	
	rowColor=vector(mode="character",length=nrow(my))
	rowColor2=vector(mode="character",length=nrow(my))
	cluster=vector(mode="character",length=nrow(my))
# 	unlink(paste(outdir,"Cluster_",str_replace(fname,".txt",""), tag, "*.txt", sep = ""))#delete the old files
	for(j in 1:N.cluster){
		sub.tf=cutree(h,k=N.cluster)==j
		rowColor[sub.tf]=colorL[j]
		if(grepl(colorL[j], "_")){#whether the color word has "_" or not
                    rowColor2[sub.tf] = unlist(strsplit(colorL[j], "_"))[1]#this is to remove _1, _2 to use the color words
                }else{
                    rowColor2[sub.tf] = rowColor[sub.tf]
                }
# 		clustername = paste(outdir, "Cluster_",str_replace(fname,".txt",""), tag, ".",colorL[j],".",j,".txt",sep="")
# 		print(clustername)
# 		write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
		#write.table(names(sub.tf[sub.tf==TRUE]),clustername,col.names=F,row.names=F)
		cluster[sub.tf]=j
	}
	
	#no need to specifically save the files
# 	write.table(rowColor, paste(outdir, "AllClusters_rowColor",str_replace(fname,".txt",""), tag, "_hclust.txt", sep = ""), quote = FALSE, sep = "\n", row.names = F, col.names = F)#save the clustering colors
	
	
# 	cluster <- as.matrix(cluster)
# 	rownames(cluster)<-rownames(my)
# 	imgDat = t(my[h$labels[h$order],])
# 
# 	if (flagArrangeColumn ==1) imgDat = imgDat[h2$labels[h2$order],]	# cluster column as well
# 	#P =      mlabel[h$labels[h$order],]
# 	imgDat[imgDat>maxval]=maxval
# 	imgDat[imgDat<minval]=minval
# 
# 	#color code
# 	#par(mar=c(10,33,5,0))
# 	mycol = colorL[1:N.cluster]
# 	cluster = cluster[h$labels[h$order],]


	#image(t(as.matrix(as.numeric(cluster))),col=mycol,  axes=FALSE)
	#box()
	#par(mar=c(10,1,5,1))
					
# 	#-------------------it has used the rowColor here in RowSideColors-----------------------------				
# 	heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
# 		 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, RowSideColors=rowColor2,
# 		 cexCol=1.5,margins=c(15,5) )#here RowSideColors=rowColor2 instead of rowColor

# 	dev.off()

# 	pngname2=paste(outdir,fname,".Dendro2.png",sep="")
# 	png(pngname2,width=1000,height=1200)
# 	    heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
# 		 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, #RowSideColors=F,
# 		 cexCol=1.5,margins=c(15,5) )
# 	dev.off()
	
	return(rowColor)
}

#' Calculate the adjusted rand index
#'
#' This function is to calculate the performance of the algorithm by 5 metrics, including the Rand index, HA, MA, FM and Jaccard. 
#'
#' @param w the ground-truth clusters
#' @param rowColor the predicted clusters
#' 
#' @examples
#' finalmetrics = ARI(gtc, finalrowColor)
#'
#' @export
ARI <- function(w, rowColor){#w: the ground-truth clusters; rowColor: the predicted clusters

	truelabel = w$cellType#get the label (categorical)

	truec = truelabel#make a copy
	levels(truec) = c(1:length(levels(truec)))#convert the categorical to numeric

	truecl = as.numeric(as.character(truec))#convert to numeric vector
	
	#convert a (categorical) vector to a numeric vector 
	t1 = factor(rowColor)#rowColor is a vector like "red red purple brown ..."
	levels(t1) = c(1:length(levels(t1)))#convert the categorical to numeric
	t1 = as.numeric(as.character(t1))
	
	metrics = adjustedRand(truecl, t1)#the ARI performance metrics
	
	return(metrics)
}
	

#' Run SHARP for single-cell RNA data clustering
#'
#' SHARP: \strong{S}ingle-cell RNA-Seq \strong{H}yper-fast and \strong{A}ccurate clustering via ensemble \strong{R}andom \strong{P}rojection. 
#'
#' @param E input single-cell expression matrix
#' @param gtc the ground-truth clusters
#' @param K number of applications of random projection
#' @param p the dimension to be reduced to
#' 
#' @examples
#' enresults = SHARP_test(E, gtc, K = 15)
#'
#' @export
# SHARP_test <- function(Files, isSparseFiles, K, p){
SHARP_test <- function(E, gtc, K, p){
# 	start_time0 <- Sys.time()#from scratch
# N.cluster=39; #change this one######## pre-determined in the later part by the unique number of clusters
# 	title="scRNA-Seq Clustering"
# 	flagArrangeColumn = 0

	colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan", "turquoise", "pink","khaki","magenta", "violet", "salmon", "goldenrod", "orchid", "seagreen", "slategray", "darkred", "darkblue", "darkcyan", "darkgreen", "darkgray", "darkkhaki", "darkorange", "darkmagenta", "darkviolet", "darkturquoise", "darksalmon", "darkgoldenrod", "darkorchid", "darkseagreen", "darkslategray", "deeppink", "lightcoral", "lightcyan")
	
# 	if(missing(Files)){
# # Files = c("Treutlein/GSE52583_FPKM_new.txt")
# # Files = c("Treutlein/GSE52583_TPM.txt")
# #             Files = c("Wang/Pancreas_Julia.published.dataTpms_selected.txt")
# # Files = c("Zeisel/GSE60361_C1-3005-Expression.tpms.txt")
# # Files = c("Zeisel/GSE60361_C1-3005-Expression_new.txt")
# # Files = c("Zeisel/Zeisel_fromSIMLR.txt")
# 
# # Files = c("Zeisel/Zeisel_expression_mat.txt")
# # Files = c("Zeisel/Zeisel_expression_mRNA_cpm.txt")
# 
# # Files = c("mECS/mECS_fromSIMLR.txt")
# # Files = c("Kolod/Kolod_fromSIMLR.txt")
# # Files = c("Usoskin/Usoskin_fromSIMLR.txt")
# # Files = c("Enge/tpms4MATLAB_selected.txt")#you need to make transpose of it because it is feature * samples
# # Files = c("GSE63473/Retina_all_exp.txt")#too large; may not be working due to memory problems
# # Files = c("GSE63473/Retina_all_exp_sparse.txt")
# # Files = c("GSE63473/Retina_900genes_exp_sparse.txt")
# # Files = c("GSE63473/Retina_1200genes_exp_sparse.txt")#less samples; more dense
# # Files = c("Baron_h1/Baron_human1_cpm.txt")#human1
# # Files = c("Baron_h2/Baron_human2_cpm.txt")#human2
# # Files = c("Baron_h3/Baron_human3_cpm.txt")#human3
# # Files = c("Baron_h4/Baron_human4_cpm.txt")#human4
# # Files = c("Baron_m1/Baron_mouse1_cpm.txt")#mouse1
# # Files = c("Baron_m2/Baron_mouse2_cpm.txt")#mouse2
# # Files = c("Goolam/Goolam_expression_log10.txt")
# #  Files = c("Macosko/retina_cpm_1800genes.txt")#human
# #  Files = c("Macosko/retina-tpms_1800genes_sparse.txt")
# # Files = c("Macosko/retina_cpm.txt")#human
# # Files = c("Darmanis/Darmanis_cpm_log10.txt")
# # Files = c("Deng/Deng_TPM.txt")#RPKM can not be comparable across samples
# # Files = c("Deng/Deng_TPM_normalized.txt")#RPKM can not be comparable across samples
# # Files = c("PBMC/PBMC_CPM_filtered0.txt")
# # Files = c("Yan/Yan_TPM.txt")
# # Files = c("Yan/Yan_RPKM.txt")
# # Files = c("Yan/Yan_TPM.txt")
# # Files = c("Klein/Klein_cpm.txt")
# 	 }
# 	 
# 	fname = Files
# 	print("Reading input expression matrix files...", quote = FALSE)
# 	if(missing(isSparseFiles) || (isSparseFiles == FALSE)){#default is no sparse files
# 	  E = read.delim(fname,row.names=1,header=T)#read the input file
# 	}else{
# 	  E = readMM(fname)#read the sparse matrix file
# 	}

	#timing
	start_time <- Sys.time()#we exclude the time for loading the input matrix
	
	ngenes = nrow(E)#number of genes
	ncells = ncol(E)#number of cells
	
# 	E = log2(E + 1)#log transform
	
	if(missing(K)){#default times of random projection
	  K = 15#K times of random projection
	}
# for (ff in 1:length(Files)){
# 	fname = Files
# 	E = read.delim(fname,row.names=1,header=T)#read the input file
# 	E = readMM(fname)#read the sparse matrix file

# 	E = read.delim(fname,header=T, check.names = F)
# 	E = read.delim(fname,row.names= NULL,header=T, check.names = F)
	
	if(missing(p)){#default dimensions to be reduced
	  p = ceiling(log2(ncol(E))/(0.2^2))#reduced 100 times of dimensions; about 200-dim
# 	  p = ceiling(log2(ncol(E))/(0.15^2))#reduced 100 times of dimensions; about 200-dim
	}
	
# 	dir1 = unlist(strsplit(Files, "/"))[1]#the folder
# 	f1 = unlist(strsplit(Files, "/"))[2]#the folder
	
	
# 	tcfile = "sampleID_sampleTitle_cell_type4MATLAB.txt"
# 	tcfile = "id2celltype.txt"
# 	tcfile = "id2celltype_900genes.txt"
# 	tcfile = "id2celltype_1200genes.txt"
# 	tcfile = "id2celltype_1500genes.txt"
# 	tcfile = "id2celltype_1800genes.txt"
# 	file1 = paste(dir1, "/", tcfile , sep = "")
# 	gtc = read.delim(file1, check.names = F)#read the file containing the ground-truth clusters
	
# 	#number of clusters
# 	if(missing(N.cluster)){#default number of clusters
# 	  N.cluster=length(unique(gtc$cellType))#number of clusters pre-determined
# 	}

# 	outdir = paste("allpara_Results/", dir1, "/", sep = "")
# 	if(!dir.exists(outdir)){#the folder to store the results
# 	  dir.create(outdir, recursive= T)
# 	}
# 	print(paste("For Dataset: ", dir1, sep = ""), quote = FALSE)
	

# 	E = log2(E + 1)#log transform
	
# 	dir1 = "Enge"
# 	dir1 = unlist(strsplit(Files, "/"))[1]#the folder
# # 	f1 = unlist(strsplit(Files, "/"))[2]#the folder
# 	
# 	outdir = paste("Results/", dir1, "/", sep = "")
# 	if(!dir.exists(outdir)){#the folder to store the results
# 	  dir.create(outdir)
# 	}
# 	
# # 	tcfile = "sampleID_sampleTitle_cell_type4MATLAB.txt"
# 	tcfile = "id2celltype.txt"
# # 	tcfile = "id2celltype_900genes.txt"
# # 	tcfile = "id2celltype_1200genes.txt"
# 	file1 = paste(dir1, "/", tcfile , sep = "")
# 	gtc = read.delim(file1, check.names = F)#read the file containing the ground-truth clusters
# 	
	
	
# 	#For different clusters (either small number of clusters or large number of clusters)
# 	# N.cluster = 39
# 	cL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")
# 	if(N.cluster <= length(cL)){
# 	  colorL = cL
# 	}else{
# 	  colorL = vector(mode = "character", length = N.cluster)
# 	  for (j in 1:N.cluster){
# 	    tmp = j%%length(cL)#residual
# 	    if (tmp != 0){
# 	      k = floor(j/length(cL)) + 1#how many rounds
# 	      colorL[j] = paste(cL[tmp], "_", k, sep = "")
# 	    }else{
# 	      k = floor(j/length(cL))#how many rounds
# 	      colorL[j] = paste(cL[tmp + length(cL)], "_", k, sep = "")
# 	    }
# 	  }
# 	}
# 	# print(colorL)
	
	
	print(paste("Number of cells: ", ncol(E), sep = ""), quote = FALSE)
	print(paste("Number of genes: ", nrow(E), sep = ""), quote = FALSE)
	print(paste("Ground-truth Number of clusters: ", length(unique(gtc$cellType)), sep = ""), quote = FALSE)
# 	p = ceiling(nrow(E)/200)#reduced 100 times of dimensions; about 200-dim
# 	p = ceiling(log(nrow(E))/(0.2^2))#reduced 100 times of dimensions; about 200-dim
# 	p = ceiling(log2(ncol(E))/(0.2^2))#reduced 100 times of dimensions; about 200-dim
	print(paste("The dimension has been reduced from ", nrow(E), " to ", p, sep=""), quote = FALSE)
	entag = paste("_enRP",p)
	allrpinfo <- vector("list", length = K)#declare a matrix of lists
	enresults = list()
# 	colnames(rpinfo) = c("E", "tag", "rowColor", "metrics")
	enrp<- matrix('0', ncol(E), K)#the ensemble results after several random projection;namely several rowColor's
	for (k in 1:K){
# 	  print("=========================================================================", quote = FALSE)
	  print(paste("The ", k, "-th time of random projection", sep=""), quote = FALSE)
# 	  inE = data.matrix(E)#conver from sparse matrix to normal matrix
	  newE = RPmat(E, p)#do the RP;the result is a list
	  E1 = t(newE$projmat)#for those which need tranpose
# 	  E1 = log10(oldE1 + 1)
	  
	  tag = paste("_RP", p, "_", k,  sep="")#k is the application times of random projection; p is the reduced dimension
# 	  rowColor= getrowColor(E1, tag, outdir, colorL)#hierarchical clustering
	  rowColor= getrowColor(E1, tag, colorL)#hierarchical clustering
	  metrics= ARI(gtc, rowColor)#performance evaluation
	  print(metrics)
	  
	  enrp[,k] = rowColor
	  
	  rpinfo = list()
	  #for different parameters, we do not need to save individual randome matrices
# 	  rpinfo$rpmat = E1#the after-random-projected matrix
# 	  rpinfo$R = newE$R#the random matrix
	  rpinfo$tag = tag#tag
	  rpinfo$rowColor = rowColor#the resulting clusters
	  rpinfo$metrics = metrics#the performance for each individual RPs
	  rpinfo$N.cluster = length(unique(rowColor))
	  
	  rpname = paste("RP_", k, sep = "")
	  allrpinfo[[rpname]] = rpinfo
# 	  rpinfo[[k]][[E]] = list(E)
# 	  rpinfo[[k]][[tag]] = list(tag)
# 	  rpinfo[[k]][[rowColor]] = list(rowColor)
# 	  rpinfo[[k]][[metrics]] = list(metrics)
	 }
	 
# 	 finalrowColor = wMetaC(enrp, N.cluster)
	 finalrowColor = wMetaC(enrp)
	 finalmetrics = ARI(gtc, finalrowColor)#performance evaluation
	 print("The ensemble performance metrics are:", quote = FALSE)
	 print(finalmetrics)
	 
	 
# 	save(enrp, file=paste(outdir,"enrp_", K, "times.RData", sep = ""))
# 	save(finalrowColor, file = paste(outdir,"finalrowColor_", K, "times.RData", sep = ""))
# 	save(finalmetrics, file = paste(outdir,"finalmetrics_", K, "times.RData", sep = ""))
	####################################
	end_time <- Sys.time()

# 	t = end_time - start_time#only the algorithm running time; excluding the loading time
# 	t0 = end_time - start_time0#all time
	
	t <- difftime(end_time, start_time, units='mins')#difference time in minutes
# 	t0 <- difftime(end_time, start_time0, units='mins')
	
	print(t, quote = FALSE)
# 	print(t0, quote = FALSE)
	####################################
# 	enresults$enrp = enrp#we have already saved the rowColor for each individual RP, thus we do not need to save enrp
# 	enresults$dataname = dir1
	enresults$finalrowColor = finalrowColor
	enresults$finalmetrics = finalmetrics
	enresults$N.pred_cluster = length(unique(finalrowColor))
	enresults$reduced_dim = p
	enresults$K = K
	enresults$allrpinfo = allrpinfo
	enresults$N.cluster = length(unique(gtc$cellType))
	enresults$N.cells = ncells
	enresults$N.genes = ngenes
	enresults$time = t
# 	enresults$cputime = t0
# 	save(enresults, file=paste(outdir,"enresults_", p, "dim_", K, "times.RData", sep = ""))
	return(enresults)
}

# 	enrowColor=vector(mode="character",length=nrow(my))
# 	cluster=vector(mode="character",length=nrow(my))
# 	unlink("Results/Cluster*.txt")
# 	for(j in 1:N.cluster){
# 		sub.tf=finalrowColor==j
# 		enrowColor[sub.tf]=colorL[j]
# 		clustername = paste("Results/Cluster",str_replace(fname,".txt",""), entag, ".",colorL[j],".",j,".txt",sep="")
# 		print(clustername)
# 		write.table(names(sub.tf[sub.tf==TRUE]),clustername,quote=F,col.names=F,row.names=F)
# 		#write.table(names(sub.tf[sub.tf==TRUE]),clustername,col.names=F,row.names=F)
# 		cluster[sub.tf]=j
# 	}
# 	
# 	cluster <- as.matrix(cluster)
# 	rownames(cluster)<-rownames(my)
# # 	imgDat = t(my[h$labels[h$order],])
# # 	#P <- mlabel[h$labels[h$order],]
# # 	#if (cl==2)	
# # 	if (flagArrangeColumn ==1) imgDat = imgDat[h2$labels[h2$order],]	# cluster column as well
# # 	#P =      mlabel[h$labels[h$order],]
# # 	imgDat[imgDat>maxval]=maxval
# # 	imgDat[imgDat<minval]=minval
# 
# 	#color code
# 	#par(mar=c(10,33,5,0))
# # 	mycol = colorL[1:N.cluster]
# # 	cluster = cluster[h$labels[h$order],]
# 	#image(t(as.matrix(as.numeric(cluster))),col=mycol,  axes=FALSE)
# 	#box()
# 	#par(mar=c(10,1,5,1))
# 							    
# 	#image(imgDat, col=bluered(length(bk)-1),axes=FALSE )
# 	#axis(1, at=seq(0,1,1/(ncol(mExp)-1)), labels=colnames(mExp), las=2, tick=FALSE,cex.axis=2)
# 	#axis(1, at=seq(0,1,1/(ncol(mExp)-1)), labels=rownames(imgDat), las=2, tick=FALSE,cex.axis=2)
# 	#box()
# 	#par(mar=c(10,1,5,14))
# 	    heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
# 		 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, RowSideColors=enrowColor,
# 		 cexCol=1.5,margins=c(15,5) )
# 
# 	#P[ (P>-3)&(P< -1) ] = 0
# 	#P[ (P<3)&(P> 1) ] = 0
# 	#P[P< -1]<- -1										    
# 	#P[P> 1]<- 1										    
#     #image(as.matrix(t(P)), col=colorRampPalette(c("blue","white","red"))(100), axes=FALSE )
#     #axis(1, at=seq(0,1,1/(ncol(P)-1)), labels=colnames(P), las=2, tick=FALSE,cex.axis=2)
#     #box()
# 	dev.off()
# 
# 	pngname2=paste(outdir, fname, entag, ".Dendro2.png",sep="")
# 	png(pngname2,width=1000,height=1200)
# 	    heatmap.2(as.matrix(my), Rowv=dend, Colv=F, scale="row",trace="none",dendrogram="row",labRow=F,
# 		 col=colorpanel(length(bk)-1,"blue","white","red"), key=T,keysize=0.7, #RowSideColors=F,
# 		 cexCol=1.5,margins=c(15,5) )
# 	dev.off()
	
# }
# stop("bawk")

# end_time <- Sys.time()
# 
# t = end_time - start_time
# print(t)

A = 100#number of times running SHARP
# K = c(3:9,seq(10, 60, by=5))#ensemble size
# K = c(seq(45, 60, by=5))#Enge2
# K = seq(5, 50, by = 5)
# K = seq(35, 50, by = 5)
K = 15
#  D = 20862#Enge
# D = 20490#Wang
#  D = 13473#Kolod
# D = 19487#Zeisel
# D = 23364#Treutlein
# D = 8990#
# Dim = ceiling(D/c(seq(20, 200, by=20)))
# Dim = ceiling(D/100)
#  Dim = ceiling(log2(2282)/(0.2^2))#default is the log2(N)/(0.2^2)#Enge
# Dim = ceiling(log2(479)/(0.2^2))#default is the log2(N)/(0.2^2)#Wang
#  Dim = ceiling(log2(704)/(0.2^2))#default is the log2(N)/(0.2^2)#Kolod
# Dim = ceiling(log2(3005)/(0.2^2))#default is the log2(N)/(0.2^2)#Zeisel
# Dim = ceiling(log2(182)/(0.2^2))#default is the log2(N)/(0.2^2)#Zeisel
# Dim = ceiling(log2(1937)/(0.2^2))#default is the log2(N)/(0.2^2)#Baron human1
# Dim = ceiling(log2(2860)/(0.2^2))#default is the log2(N)/(0.2^2)#Macosko 1800genes
# Dim = ceiling(log2(4385)/(0.2^2))#default is the log2(N)/(0.2^2)#Macosko 1800genes
#  Dim = ceiling(log2(44808)/(0.2^2))#default is the log2(N)/(0.2^2)#Macosko 1800genes

##########Input matrix files###########

# Files = c("Wang/Pancreas_Julia.published.dataTpms_selected.txt")
# Files = c("Zeisel/GSE60361_C1-3005-Expression.tpms.txt")
# Files = c("Zeisel/GSE60361_C1-3005-Expression_new.txt")
# Files = c("Zeisel/Zeisel_fromSIMLR.txt")

# Files = c("Zeisel/Zeisel_expression_mat.txt")
# Files = c("Zeisel/Zeisel_expression_mRNA_cpm.txt")

# Files = c("mECS/mECS_fromSIMLR.txt")
# Files = c("Kolod/Kolod_fromSIMLR.txt")
# Files = c("Usoskin/Usoskin_fromSIMLR.txt")
# Files = c("Enge/tpms4MATLAB_selected.txt")#you need to make transpose of it because it is feature * samples
# Files = c("GSE63473/Retina_all_exp.txt")#too large; may not be working due to memory problems
# Files = c("GSE63473/Retina_all_exp_sparse.txt")
# Files = c("GSE63473/Retina_900genes_exp_sparse.txt")
# Files = c("GSE63473/Retina_1200genes_exp_sparse.txt")#less samples; more dense
# Files = c("Baron_h1/Baron_human1_cpm.txt")#human1
# Files = c("Baron_h2/Baron_human2_cpm.txt")#human2
# Files = c("Baron_h3/Baron_human3_cpm.txt")#human3
# Files = c("Baron_h4/Baron_human4_cpm.txt")#human4
# Files = c("Baron_m1/Baron_mouse1_cpm.txt")#mouse1
# Files = c("Baron_m2/Baron_mouse2_cpm.txt")#mouse2
# Files = c("Goolam/Goolam_expression_log10.txt")
#  Files = c("Macosko/retina_cpm_1800genes.txt")#human
#  Files = c("Macosko/retina-tpms_1800genes_sparse.txt")
# Files = c("Macosko/retina_cpm.txt")#human
# Files = c("Darmanis/Darmanis_cpm_log10.txt")
# Files = c("Deng/Deng_TPM.txt")#RPKM can not be comparable across samples
# Files = c("Deng/Deng_TPM_normalized.txt")#RPKM can not be comparable across samples
# Files = c("PBMC/PBMC_CPM_filtered0.txt")
# Files = c("Yan/Yan_TPM.txt")
# Files = c("Yan/Yan_RPKM.txt")
# Files = c("Yan/Yan_TPM.txt")
# Files = c("Klein/Klein_cpm.txt")
#########################################
# dset = "Enge"
# dset = "Wang"
# dset = "Kolod"
# dset = "Zeisel"
# dset = "Treutlein"
# dset = "mECS"
# dset = "Baron_h1"
# dset = "Macosko"

# # eps = c(0.3, 0.25, 0.2, 0.15, 0.1, 0.05)
# allresults <- vector("list", length =  length(K))#declare a matrix of lists
# for(k in K){
# #   allresults = list()
# #   for(i in 1:length(eps)){
#     info = vector("list", length = A)
#     kname = paste("enSize_", k, sep = "")
# #     dim = ceiling(log2(2282)/(eps[i]^2))
# #     kname = paste("dim_", dim, sep = "")
#     for(j in 1:A){#times of using SHARP
#       print("=========================================================================", quote = FALSE)
#       print(paste("SHARP: Ensemble Size-", k, ", Run-", j, sep = ""), quote = FALSE)
# #       enresults = SHARP_test(Files, isSparseFiles, K, p, N.cluster)
# # 	enresults = SHARP_test(Files, K = K[k], p = dim)
# 	enresults = SHARP_test(Files, K = k)#using the default Dim
# #       enresults = SHARP_test(K = k, p = i, isSparseFiles = T)
# #       enresults = SHARP_test(K = k, p = i, isSparseFiles = T, N.cluster = 10)
# #       outname = paste("enresults_ensize", k, "_dim", enresults$reduced_dim, "_run",j, sep = "")
#       jname = paste("Run_", j, sep = "")
#       info[[jname]] = enresults
#     }
#     allresults[[kname]] = info
# #   }
# # #   dset = enresults$dataname#just use the last one to retrieve the dataset name
# #   
# }
# dset = enresults$dataname#just use the last one to retrieve the dataset name
# # save(allresults, file=paste("allpara_Results/",dset,"/SHARP_",dset,"_differentDim.RData", sep = ""))
# saveRDS(allresults, file=paste("allpara_Results/SHARP_", dset, ".rds", sep = ""))
# 
# print("Done!", quote = FALSE)
