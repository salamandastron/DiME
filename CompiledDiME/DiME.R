# We first load the expression data as well as the library functions needed for
# jackknife resampling of correlation coefficients into the R workspace
source("jackknifeCoexpression.R")
exprs <- read.csv("G4_expression.csv",
                    sep = ",",
                    header = TRUE)
nfeat <- dim(exprs)[1]
nsamples <- dim(exprs)[2]
exprs <- t(exprs)

cor.data <- cor.jack(exprs)

upper.r <- as.vector(cor.data[upper.tri(cor.data)])
upper.ind <- which(upper.tri(cor.data))
sorted.r <- sort((upper.r), index.return = TRUE)
nedges <- floor(length(upper.r)*0.001)
nedges.upper <- floor(0.5*nedges)
nedges.lower <- nedges - nedges.upper
selected.ind <- c(upper.ind[sorted.r$ix[1:nedges.lower]],
upper.ind[sorted.r$ix[(length(upper.r)-nedges.upper+1):length(upper.r)]])

# Now connect these correlation pairs to form an adjacency matrix for the network,
# and delete isolated nodes (genes with no selected coexpression partners)
adj <- matrix(0, nfeat, nfeat)
rownames(adj) <- as.character(rownames(cor.data))
colnames(adj) <- as.character(colnames(cor.data))
adj[selected.ind] <- 1
ind <- lower.tri(adj)
adj[ind] <- t(adj)[ind]
toDel <- c()
for (i in 1:nfeat) {
    if(sum(adj[i,])==0) {
    toDel <- c(toDel, i)
    }
}
adj <- adj[-toDel,-toDel]

# Build an iGraph object to represent the resulting network
require(igraph)
g <- graph.adjacency(adj,mode = "undirected")
edges <- get.edgelist(g)
write.table(edges,"G4_network.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# Disease Module Extraction using DiME

dyn.load(dynlibfile) # Change .so to .dll if you are using Windows
communityExtraction <- function(dim, popsize, adj, verbose){
out <- .C("commextr",
    d = as.integer(dim),
    ps = as.integer(popsize),
    adjm = (as.vector(adj)),
    res = as.integer(as.vector(rep(0,dim))),
    f = as.integer(0),
    v = as.integer(verbose))
    return(out$res)
}

W.score <- function(adj_matrix, th) {
    selected_nodes <- which(th!=0)
    unselected_nodes <- which(th==0)
    S <- length(selected_nodes)
    if (S==0 | S==length(th)) {
    return(0)
    }
    S_C <- length(unselected_nodes)
    Aij_S <- adj_matrix[selected_nodes,selected_nodes]
    Aij_cS <- adj_matrix[selected_nodes,unselected_nodes]
    OS <- sum(sum(Aij_S))
    BS <- sum(sum(Aij_cS))
    W <- S*S_C*(OS/(S*S) - BS/(S*S_C))
    return(W)
}

# Next, read the network file we have saved in the previous step, and
# prepare a copy of the network with vertex names in numbers in order
# to allow B-score calculation later on

library(igraph)
edges <- read.table("G4_network.txt",
            sep = "\t",
            header = FALSE)
edges <- cbind(as.character(edges[,1]), as.character(edges[,2]))
g <- graph.edgelist(as.matrix(edges), directed = FALSE)
IDs <- V(g)$name
inds <- c(1:length(V(g)))
V(g)$name <- as.character(c(1:length(V(g))))
edgelist <- get.edgelist(g)
edgelist <- data.frame(as.character(edgelist[,1]), as.character(edgelist[,2]))
write.table(edgelist,
            "G4_edgelist.txt",
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

# Perform the DiME algorithm iteratively on the network,
# removing newly extracted modules in each iteration

nRuns <- 50 # Number of extraction trials to be done for each community
Bscores <- c()
Wscores <- c()
bestCommunities <- list()
level <- 0
while (level <= 1) {
    level <- level + 1
    currentBest <- rep(0,length(V(g)))
    # Running in verbose mode
    writeLines(paste("Now extracting community ",as.character(level),"...",sep=""))
    adj <- get.adjacency(g)
    for (i in 1:nRuns) {
        res <- communityExtraction(dim = length(V(g)),
                popsize = 50,
                adj = adj,
                verbose=0)
        if (W.score(adj,res) > W.score(adj,currentBest)) {
                    currentBest <- res
        }
        cat(".")
    }
    wbest <- W.score(adj, currentBest)
    cat("\n")
    if(sum(currentBest) > 2) {
        bestCommunities[[level]] <- as.numeric(V(g)$name[currentBest!=0])
        Wscores[level] <- wbest
    } else {
        cat("\n...No feasible community can be extracted any more. Extraction done.\n")
        break
    }
    g <- delete.vertices(g,V(g)[currentBest!=0])
}
    
write.table(Wscores,
            "G4_Wscores.txt",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)
            
# Write a list of modules extracted for B-score evaluation
temp <- lapply(bestCommunities, write, "G4_commlist.dat", append = TRUE, ncolumns = 1000, sep = "\t")

for (x in 1:length(Wscores)) {
    filename = paste("G4_community_", as.character(x), ".txt", sep = "")
    community = as.character(IDs[bestCommunities[[x]]])
    write.table(community,
    	        filename,
                sep = "\t",
                col.names = FALSE,
                row.names = FALSE,
                quote = FALSE)
}
