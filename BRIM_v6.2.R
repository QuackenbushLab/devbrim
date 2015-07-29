#Change from v6.1: Corrected error in the case where project=FALSE. T0 was not
#populated correctly in v6.1.
#Change from v6.1: vectorized index calls in lines 108 and 141
#Changes from v6.0: New function max.component.bipartite, which takes a directed bipartite
#network and returns the edgelist of giant connected component with the red and blue
#nodes in separate columns, as required by BRIM and fast.brim
#Changes from v6.0: Throws an error if there is more than one connected component
#Changes from v6.0: Community IDs are relabeled so that they are sequential
#with no skips (eg, seq(1,Number communities,by=1))
#Changes from BRIM2_v2: Updated Q_com out to handle empty communities.
#Changes from BRIM2_v3/4: changed if statement for making a new community:
# old if statement - if(btr[g] < 0 && length(unique(R[,2])))
#Used as.integer(factor()) to obtain numerical code for column/row ids, per
#best practice as stated in the examples section of factor

###
#This function is based on the bipartite modularity
#as defined in "Modularity and community detection in bipartite networks"
#by Michael J. Barber, Phys. Rev. E 76, 066102 (2007)
#This function uses a slightly different implementation from the paper. It does not
#use the "adaptive BRIM" method for identifying the number of modules. Rather,
#simply continues to iterate until the difference in modularity between iterations
#is less that 10^-4. Starting from a random initial condition, this could take some time.
#Use fast.brim for quicker runtimes, it intializes the blue node memberships by projecting
#the blue nodes into a unipartite "blue" network and then identify communities in that
#network using a standard unipartite community detection algorithm run on the projected network.
#This was suggested by M. Barber, as the bisection method discussed in section D of his paper is
#rather tricky. This also speeds up run times.
###

# Input two column vectors, reds and blues, which together form the
# edgelist with node names
# edgelist should have one node type in each column #
#### WARNING: This code requires nrows > ncols ####
library(igraph)
library(Matrix)
library(nnet)

#X is either the upper right or lower left quadrant of the
#block off-diagonal adjacency matrix representing a bipartite network as
#described as A_tilde in the Barber paper or an edgelist.
#T0 is the intial community assignment for each "blue" node, assuming there
# are more blues than reds
BRIM = function(X,edge.list=FALSE,T0=cbind(1:q,rep(1,q)),weights=1){

  if(edge.list==TRUE){
    #Convert the edgelist to a sparseMatrix object
    esub <- X
    #make sure there's only one connected component
    g.component.test <- graph.data.frame(esub,directed=FALSE)
    if(!is.connected(g.component.test)){
        stop("More than one connected component,
              method requires only one connected component")
    }
    reds <- as.integer(factor(esub[,1]))
    red.names <- levels(factor(esub[,1]))
    blues <- as.integer(factor(esub[,2]))
    blue.names <- levels(factor(esub[,2]))
    #ensure that nrows > ncols
    if(length(red.names) < length(blue.names)){
      stop("Adjacency matrix dimension error: This code requires nrows > ncols")
      }

    #Sparse matrix with the upper right block of the true Adjacency matrix. notices the dimension is reds x blues
    A = sparseMatrix(i=reds,j=blues,x=weights,dims=c(length(unique(reds)),length(unique(blues))),index1=TRUE);
    rownames(A) <- red.names
    colnames(A) <- blue.names
  }


  if(edge.list==FALSE){
    #convert input matrix to sparse matrix, if it's not already
    A = Matrix(as.matrix(X),sparse=TRUE)
    red.names <- rownames(A)
    blue.names <- colnames(A)
    if(nrow(A) < ncol(A)){
    stop("Adjacency matrix dimension error: This code requires nrows > ncols")
    }
  }

p = nrow(A)
q = ncol(A)
N = p+q

#Sparese matrix with the upper right block of the true Adjacency matrix. notices the dimension is reds x blues
ki = rowSums(A)
dj = colSums(A)
m = sum(ki) # m = sum(dj) too
#initialize community assignments for red and blue nodes.
Tmat <- T0
R = cbind(1:p,rep(0,length=p))
cs = sort(unique(Tmat[,2]))
#variables to track modularity changes after each loop
Qhist <- vector();
Qnow  <- 0
deltaQ <- 1
#______________________________________
while(round(deltaQ,digits=4) > 0){
  btr <- BTR <- bt <- BT <- vector();

#calculate T tilde
for(i in 1:p){
  if(i %% 2500 == 0){print(paste(i/p,"percent through iteration"))}
  #find the optimal community for node i
  bt <- rep(0,length(cs))
  for(k in cs){
    if(length(Tmat[Tmat[,2]==k,2]) != 0){
    ind <- Tmat[,2] == k
    #bt[k] = sum(A[i,ind] - (ki[i]*dj[ind])/m)
    bt[k] = sum((A[i,] - (ki[i]*dj)/m)[ind])
    }
    }
  #note that which.max returns the FIRST max if more than one max
  h = which.max(bt)
    if(bt[h] < 0){
     print("making new comm")
     R[i,2] <- max(Tmat[,2])+1
     cs <- sort(c(cs,R[i,2]))
    }
  #if(length(h) > 1){h <- sample(h,1)}
  if(bt[h] > 0){
  R[i,2] <- h # assign blue vertex i to comm k such that Q is maximized
  bt[-h] <- 0 # BTR is zero if i is not in k (see definition of Q)

  }

  #BT <- rbind(BT,bt)
}
#calculate R tilde, i.e., B_transpose * R
for(j in 1:q)
{
  #initialize jth row of (B_transpose) * R
  btr = rep(0,length(cs))
  #calculate Q for assigning j to different communities
  for(k in cs)
  {
    #if node j is in community k, else BTR[j,k] = 0
    if(length(R[R[,2]==k,2]) != 0)
    {

    ind <- R[,2] == k
    #btr[k] = sum(A[ind,j]-(ki[ind]*dj[j])/m)
    btr[k] = sum((A[,j]-(ki*dj[j])/m)[ind])
    }
  }

  g = which.max(btr)
    #if there is no comm. assignment to increase modularity, make
    #a new community with that node only
    if(btr[g] < 0)
    {
      print("making new comm")
      Tmat[j,2] <- max(R[,2])+1
      btr <- rep(0,length(cs)+1)
      cs <- c(cs,Tmat[j,2])
      btr[length(cs)] <- sum(A[,j]-(ki*dj[j])/m)
      print(btr)
    }

    if(btr[g] > 0)
    {
    Tmat[j,2] <- g
    btr[-g] <- 0
    }

  #add another column to BTR if a new community is added
  if( !is.vector(BTR) && dim(BTR)[2] < length(btr)){BTR <- cbind(BTR,0) }
  BTR <- rbind(BTR,btr)
}

Tt =  t(sparseMatrix(i=Tmat[,1],j=Tmat[,2],x=1,dims=c(q,length(cs)),index1=TRUE))
Qthen <- Qnow
Qcom <- diag(Tt %*% BTR)/m
Qnow <- sum(Qcom)
Qhist = c(Qhist,Qnow)

print(paste("Q =",Qnow,sep=" "))
  if(round(Qnow,digits=4) != 0 && round(Qnow,digits=4) != 0){
  deltaQ = Qnow - Qthen
  }
}
  qcom_temp <- cbind(Qcom,sort(unique(cs)))
  #drop empty communities
  qcom_out <- qcom_temp[qcom_temp[,1] > 0,]
  #if communities were dropped, relabel so community labels can function
  #as row/column indices in the modularity matrix, B_ij.
  if(nrow(qcom_out) < nrow(qcom_temp)){
    qcom_out[,2] <- as.integer(factor(qcom_out[,2]))
    R[,2] <- as.integer(factor(R[,2]))
    Tmat[,2] <- as.integer(factor(Tmat[,2]))
  }
  colnames(qcom_out) <- c("Qcom","community")
  out = list(Qcoms=qcom_out,modularity=Qhist,red.memb=
               data.frame(red.names=red.names[R[,1]],com=R[,2]),blue.memb=
               data.frame(blue.names=blue.names[Tmat[,1]],com=Tmat[,2]))
  return(out)
}

fast.brim <- function(elist,cs.method="LCS",project=TRUE,weights=1)
{
  #make sure there's only one connected component
  g.component.test <- graph.data.frame(elist,directed=FALSE)
    if(!is.connected(g.component.test)){
      stop("More than one connected component detected,
              method requires only one connected component")
    }
  #cs.method is a string to specify which community detection method should be used
    G <- graph.data.frame(elist,directed=FALSE)

    #Use unipartite comm. structure method for first pass
    #project network into gene space to obtain intial community assignment for genes
    if(project){
        #SNP row indices
        reds = as.integer(factor(elist[,1]))
        red.names = levels(factor(elist[,1]))
        #Gene column indices
        blues = as.integer(factor(elist[,2]))
        blue.names = levels(factor(elist[,2]))
        N = max(blues)

        #Sparese matrix with the upper right block of the true Adjacency matrix. notices the dimension is reds x blues
        sM = sparseMatrix(i=reds,j=blues,x=1,dims=c(length(unique(reds)),length(unique(blues))),index1=T);
        #Project into gene space, projected adjacency matrix is has dim = genes x genes
        gM = t(sM) %*% sM;
        rm(sM)
        gc()
        colnames(gM) <- blue.names
        rownames(gM) <- blue.names
        G1 = graph.adjacency(gM,mode="undirected",weight=TRUE,diag=FALSE);
        #if(clusters(G1)$no > 1){print("Warning more than one component! May cause indexing error")}
        #V(G1)$name <- sort(unique(as.vector(esub[,2])))
        #remove loops and multiple edges
        gcc.initialize = simplify(max.component(G1))
    }

    #option to treat the bipartite network as if it is unipartite
    #for community initialization only
    if(!project)
      {
        gcc.initialize <- G
        blue.names = levels(factor(elist[,2]))
        #blue.indx <- V(G)$name %in% blue.names
      }

    if(cs.method=="LCS"){cs0 = multilevel.community(gcc.initialize)}
    if(cs.method=="LEC"){cs0 = leading.eigenvector.community(gcc.initialize)}
    if(cs.method=="FG"){cs0 = fastgreedy.community(gcc.initialize)}
    print(paste("modularity of projected graph",max(cs0$modularity)))

    #initial condition for genes community membership
    if(project){ T0 <- data.frame(as.integer(factor(blue.names)),membership(cs0)) }
    if(!project)
      {
        blue.membs <- membership(cs0)[blue.names]
        T0 <- data.frame(as.integer(factor(blue.names)),blue.membs)
      }
    #run BRIM using intial assignments T0
    bout <- BRIM(elist,edge.list=TRUE,T0=T0, weights=weights)

    return(bout)
}

max.component = function(g)
{
  # return largest connected component of the iGraph graph object g
  g.clust = clusters(g);
  maxclust.id = which(g.clust$csize == max(g.clust$csize))[1];
  h = induced.subgraph(g, which(g.clust$membership == maxclust.id)); # 1-indexed here
  return(h);
}

max.component.bipartite = function(g)
#This function returns an edgelist with the red nodes and blue nodes
#in separate columns, as required by BRIM and fast.brim
{
  if(!is.directed(g)){
    stop("function requires a directed graph!")
  }
  # return largest connected component of the iGraph graph object g
  g.clust = clusters(g);
  maxclust.id = which(g.clust$csize == max(g.clust$csize))[1];
  h = induced.subgraph(g, which(g.clust$membership == maxclust.id)); # 1-indexed here
  edgelist.out <- get.edgelist(h)
  return(edgelist.out);
}

edge.to.graph <- function(edgelist,return.gcc=FALSE)
{
    #edgelist should be a two column data.frame with colnames 'red' and 'blue'
    #this function produces an igraph graph object with a 'color' attribute
    #based on the colnames of edgelist. This can be accessed via
    #V(g)$color, which returns a vector indicating red/blue. Use V(g)$name
    #with V(g)$color to identify red/blue node names
    #if return.gcc=TRUE, returns the giant connected component of g
    if(sum(colnames(edgelist) %in% c("red","blue")) != 2)
    {
        stop("edgelist colnames must be labeled 'red' and 'blue'")
    }
    if(sum(is.na(edgelist$red) + is.na(edgelist$blue)) > 0)
    {
        stop("NA's detected. Remove these from edgelist")
    }
    if(sum(edgelist$red == "NA") + sum(edgelist$blue == "NA") > 0)
    {
        stop("NA's detected. Remove these from edgelist")
    }
    if(sum(edgelist$red == "") + sum(edgelist$blue == "") > 0)
    {
        stop("Empty strings detected. Remove these from edgelist")
    }
    if(sum(edgelist$red %in% edgelist$blue) > 0)
    {
        stop("edgelist contains one or more nodes that appear in both red and blue columns.
        Check to make sure network is truly bipartite.")
    }

    g <- graph.data.frame(edgelist,directed=FALSE)
    blue.indx <- V(g)$name %in% unique(edgelist$blue)
    V(g)$color <- "red"
    V(g)$color[blue.indx] <- "blue"

    if(!return.gcc){ g.out <- g}
    if(return.gcc){ gcc <- max.component(g); g.out <- gcc }
    blue.names <- V(g.out)$name[V(g.out)$name %in% unique(edgelist$blue)]
    red.names <- V(g.out)$name[V(g.out)$name %in% unique(edgelist$red)]
    edges <- edgelist[edgelist$blue %in% blue.names,]

    return(list(blues=blue.names,reds=red.names,G=g.out,edges=edges))
}

