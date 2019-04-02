library(igraph)
    
G4_exprs <- as.matrix(read.table("G4_expression.csv", sep = ",", header = TRUE))
NT_exprs <- as.matrix(read.table("NT_expression.csv", sep = ",", header = TRUE))

# Read in only modules with B-scores below the 0.001 threshold

B.tol <- 0.001
    
graph.filename <- "G4_network.txt"
e <- read.table(graph.filename,
        sep = "\t",
        header = FALSE)
e <- cbind(as.character(e[,1]), as.character(e[,2]))
g <- graph.edgelist(as.matrix(e), directed = FALSE)
IDs <- V(g)$name
adj <- get.adjacency(g)
    
Bscores <- read.table("G4_edgelist.txt.table",
        sep = " ",
        header = FALSE)
Bscores <- Bscores[,c(1,3,7)]
goodB <- which(Bscores[,3]<B.tol)
listOfGenes <- list()
count <- 1
while (count <= length(goodB)) {
    filename <- paste("G4_community_",
                      as.character(Bscores[goodB[count],1]),
                      ".txt",
                      sep = "")
    geneList <- read.table(filename,
                           sep = "\t",
                           header = FALSE)
    geneList <- as.character(geneList[,1])
    listOfGenes[[count]] <- geneList
    count <- count + 1
}

# Build an adjacency matrix with same number of rows as the number of
# statistically significant modules, and fill in each entry (i,j) with
# the number of connections between module i and j in the entire
# coexpression network

adjMod <- matrix(0,length(listOfGenes), length(listOfGenes))
for (i in 1:length(listOfGenes)) {
    for (j in i:length(listOfGenes)) {
        adjMod[i,j] <- sum(adj[listOfGenes[[i]], listOfGenes[[j]]])
    }
}

# Set diagonals to 1 in order to prevent isolated modules (modules with no
# connections to other modules) from being lost in the edgelist (these self-loops
# can easily be removed in Cytoscape)
diag(adjMod) <- rep(1, dim(adjMod)[1])

gmod <- graph.adjacency(adjMod, mode="undirected", weighted = TRUE)
emod <- cbind(get.edgelist(gmod),log(E(gmod)$weight)) # Weights are logged to prevent disproportionally "thick" edges
write.table(emod,
            "G4_modnet_edgelist.txt",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

# Record module sizes and calculate log fold-change of gene expression values
# averaged for each module as node attributes
modFC <- function(x) {
    ntExprs <- mean(2^NT_exprs[x,])
    exprs <- mean(2^G4_exprs[x,])
    return(log2(exprs/ntExprs))
}
fc <- sapply(listOfGenes,modFC)
lengths <- sapply(listOfGenes,length)
    
nodeAttrib <- data.frame(V = rep(1:length(listOfGenes)),
                         FC = fc,
                         size = lengths)
write.table(nodeAttrib,
            "G4_modnet_nodeAttrib.txt",
            sep = "\t",
            col.names = FALSE,
            row.names = FALSE,
            quote = FALSE)

