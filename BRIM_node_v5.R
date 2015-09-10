library(igraph)
library(Matrix)
#This function calculates the modularity contribution for each node
Qscore = function(X,brim.output){
  #Qscore is designed to calculate the fraction of the modularity contributed by each node 
  #to its community's modularity 
  ##X is either the upper right or lower left quadrant of the 
  #block off-diagonal adjacency matrix representing a bipartite network as
  #described as A_tilde in the Barber paper. Alternatively, X can be a two column edge list
  #brim.output is the list returned by BRIM or fast.brim   
  bo <- brim.output
  bo$blue.memb <- bo$blue.memb[order(bo$blue.memb[,"blue.names"]),]
  bo$red.memb <- bo$red.memb[order(bo$red.memb[,"red.names"]),]
  bo$Qcoms <- bo$Qcoms[order(bo$Qcoms[,"community"]),]
  brim.object <- bo
  
  mdim <- min(dim(X))    
  if(mdim==2 & max(dim(X)) > 2){
      #Convert the edgelist to a sparseMatrix object
      esub <- X
      reds = as.integer(factor(esub[,1]))
      blues = as.integer(factor(esub[,2]))
      #Sparese matrix with the upper right block of the true Adjacency matrix. notices the dimension is reds x blues
      A = sparseMatrix(i=reds,j=blues,x=1,dims=c(length(unique(reds)),length(unique(blues))),index1=TRUE);
      if(nrow(A) < ncol(A)){A <- t(A)}
      #if(max(blues) > max(reds)){blues <- reds;}
      }
    #convert input matrix to sparse matrix, if it's not already
      if(mdim > 2){A = Matrix(as.matrix(X),sparse=TRUE)}  
      if(nrow(A) < ncol(A)){A <- t(A)}
    #reds = as.vector(esub[,1],mode="integer")
    p = nrow(A)
    q = ncol(A)
    N = p+q

    #Sparese matrix with the upper right block of the true Adjacency matrix. notices the dimension is reds x blues
    ki = rowSums(A)
    dj = colSums(A)
    m = sum(ki) # m = sum(dj) too
    R1 = brim.object$red.memb
    T1 = brim.object$blue.memb
    r1 = cbind(as.numeric(factor(R1[,1])),R1[,2])
    t1 = cbind(as.numeric(factor(T1[,1])),T1[,2]) 
    Rtrans = sparseMatrix(i=r1[,2],j=r1[,1],x=1,dims=c(max(r1[,2]),length(unique(r1[,1]))),index1=TRUE);    
    T2 = sparseMatrix(i=t1[,1],j=t1[,2],x=1,dims=c(max(t1[,1]),max(t1[,2])),index1=TRUE);
    Qcoms <- brim.object$Qcoms
    Qjk = vector(length=q)
      for(j in 1:max(t1[,1])){
        if(j %% 1000 == 0){print(paste(j,t1[j,]))}
        Bj = A[,j] - (ki*dj[j])/m;
        Qjk[j] = ((Rtrans[t1[j,2],] %*% Bj)/(2*m))*(1/Qcoms[t1[j,2],1])
      }  
    Qik = vector(length=p)
      for(i in 1:max(r1[,1])){
        if(i %% 1000 == 0){print(i)}
        Bi = A[i,] - (ki[i]*dj)/m;
        Qik[i] = ((Bi %*% T2[,r1[i,2]])/(2*m))*(1/Qcoms[r1[i,2],1])  
      }    
dout = list(blue.qscore=data.frame(T1,Qjk),red.qscore=data.frame(R1,Qik))
return(dout)
}

