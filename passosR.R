Data loading and converting

The export file can be read by R using the function readLines:
  
  D <- readLines("http://www.bibliometrix.org/datasets/savedrecs.bib")
D is a large character object.

It can be converted in a data frame using the function convert2df:
  
  library(bibliometrix)

M <- convert2df(D, dbsource = "isi", format = "bibtex")
Bibliometric Analysis

T

results <- biblioAnalysis(M, sep = ";")

Functions summary and plot



S=summary(object = results, k = 10, pause = FALSE)
plot(x = results, k = 10, pause = FALSE)


Analysis of Cited References

M$CR[1]
The figure shows the reference string of the first manuscript. In this case, the separator field is sep = ".  ".



CR <- citations(M, field = "article", sep = ".  ")
CR[1:10]

To obtain the most frequent cited first authors:
  
  CR <- citations(M, field = "author", sep = ".  ")
CR[1:10]


function localCitations generates the frequency table of the most local cited authors. Local citations measure how many times an author included in this collection have been cited by other authors also in the collection.

To obtain the most frequent local cited authors:
  
  CR <- localCitations(M, results, sep = ".  ")
CR[1:10]

Authors’ Dominance ranking

The function dominance calculates the authors’ dominance ranking as proposed by Kumar & Kumar, 2008.

Function arguments are: results (object of class bibliometrix) obtained by biblioAnalysis; and k (the number of authors to consider in the analysis).

DF <- dominance(results, k = 10)
DF

The Dominance Factor is a ratio indicating the fraction of multi authored articles in which a scholar appears as first author.



Authors’ h-index

The function Hindex calculates the authors’ H-index and its variants (g-index and m-index) in a bibliographic collection.

indices <- Hindex(M, authors="BORNMANN L", sep = ";")

# Bornmann's impact indices:

indices$H

##       Author h_index g_index m_index
## 1 BORNMANN L       4       7     0.8
# Bornmann's citations

indices$CitationList

To calculate the h-index of the first 10 most productive authors (in this collection):
  
  authors=gsub(","," ",names(results$Authors)[1:10])

indices <- Hindex(M, authors, sep = ";")

indices$H


Lotka’s Law coefficient estimation

The function lotka estimates Lotka’s law coefficients for scientific productivity (Lotka A.J., 1926).

Lotka’s law describes the frequency of publication by authors in any given field as an inverse square law, where the number of authors publishing a certain number of articles is a fixed ratio to the number of authors publishing a single article. This assumption implies that the theoretical beta coefficient of Lotka’s law is equal to 2.

Using lotka function is possible to estimate the Beta coefficient of our bibliographic collection and assess, through a statistical test, the similarity of this empirical distribution with the theoretical one.

L <- lotka(results)

# Author Productivity. Empirical Distribution
L$AuthorProd

# Beta coefficient estimate
L$Beta
## [1] 3.04525
# Constant
L$C
## [1] 0.6018257
# Goodness of fit
L$R2
## [1] 0.9353053
# P-value of K-S two sample test
L$p.value
## [1] 0.08786641
The table L$AuthorProd shows the observed distribution of scientific productivity in our example.

The estimated Beta coefficient is 3.05 with a goodness of fit equal to 0.94. Kolmogorov-Smirnoff two sample test provides a p-value 0.09 that means there is not a significant difference between the observed and the theoretical Lotka distributions.

You can compare the two distributions using plot function:
  
  # Observed distribution
  Observed=L$AuthorProd[,3]

# Theoretical distribution with Beta = 2
Theoretical=10^(log10(L$C)-2*log10(L$AuthorProd[,1]))

plot(L$AuthorProd[,1],Theoretical,type="l",col="red",ylim=c(0, 1), xlab="Articles",ylab="Freq. of Authors",main="Scientific Productivity")
lines(L$AuthorProd[,1],Observed,col="blue")
legend(x="topright",c("Theoretical (B=2)","Observed"),col=c("red","blue"),lty = c(1,1,1),cex=0.6,bty="n")


Bibliometric network matrices



Bipartite networks

cocMatrix is a general function to compute a bipartite network selecting one of the metadata attributes.

For example, to create a network Manuscript x Publication Source you have to use the field tag “SO”:
  
  A <- cocMatrix(M, Field = "SO", sep = ";")
A is a rectangular binary matrix, representing a bipartite network where rows and columns are manuscripts and sources respectively.

The generic element aijaij is 1 if the manuscript ii has been published in source jj, 0 otherwise.

The j−thj−th column sum ajaj is the number of manuscripts published in source jj.

Sorting, in decreasing order, the column sums of A, you can see the most relevant publication sources:
  
  sort(Matrix::colSums(A), decreasing = TRUE)[1:5]

Following this approach, you can compute several bipartite networks:
  
  # Citation network
  
  A <- cocMatrix(M, Field = "CR", sep = ".  ")

# Author network

A <- cocMatrix(M, Field = "AU", sep = ";")

#Country network

Authors’ Countries is not a standard attribute of the bibliographic data frame. You need to extract this information from affiliation attribute using the function metaTagExtraction.

M <- metaTagExtraction(M, Field = "AU_CO", sep = ";")


# A <- cocMatrix(M, Field = "AU_CO", sep = ";")
metaTagExtraction allows to extract the following additional field tags: Authors’ countries (Field = "AU_CO"); First author of each cited reference (Field = "CR_AU"); and Publication source of each cited reference (Field = "CR_SO").

# Author keyword network

A <- cocMatrix(M, Field = "DE", sep = ";")


#Keyword Plus network
A <- cocMatrix(M, Field = "ID", sep = ";")
Etc.


Bibliographic coupling

the most frequently used coupling networks: Authors, Sources, Keywords and Countries


biblioNetwork uses two arguments to define the network to compute:
  
  analysis argument can be “collaboration”, “coupling” or “co-citation”

network argument can be “authors”, “references”, “sources”, “countries”,“keywords” or “author_keywords”

The following code calculates a classical article coupling network:
  
  
  NetMatrix <- biblioNetwork(M, analysis = "coupling", network = "references", sep = ".  ")

Articles with only a few references, therefore, would tend to be more weakly bibliographically coupled, if coupling strength is measured simply according to the number of references articles contain in common.

This suggests that it might be more practicable to switch to a relative measure of bibliographic coupling.


couplingSimilarity function 

calculates Jaccard or Salton similarity coefficient among manuscripts of a coupling network.

NetMatrix <- biblioNetwork(M, analysis = "coupling", network = "sources", sep = ";")

# calculate jaccard similarity coefficient
S <- couplingSimilarity(NetMatrix, type="jaccard")

# plot journals' similarity (with min 3 manuscripts)
diag <- Matrix::diag
MapDegree <- 3
NETMAP <- S[diag(NetMatrix)>=MapDegree,diag(NetMatrix)>=MapDegree]
diag(NETMAP) <- 0

H <- heatmap(max(NETMAP)-as.matrix(NETMAP),symm=T, cexRow=0.3,cexCol=0.3)


#

#Bibliographic co-citation


Using the function biblioNetwork, you can calculate a classical reference co-citation network:
  
  # NetMatrix <- biblioNetwork(M, analysis = "co-citation", network = "references", sep = ".  ")
  
  
  #Bibliographic collaboration
  
  
  
  Using the function biblioNetwork, you can calculate an authors’ collaboration network:
  
  # NetMatrix <- biblioNetwork(M, analysis = "collaboration", network = "authors", sep = ";")
  
  
  or a country collaboration network:
  
  # NetMatrix <- biblioNetwork(M, analysis = "collaboration", network = "countries", sep = ";")
  
  
  ##Visualizing bibliographic networks
  
  
  Country Scientific Collaboration

# Create a country collaboration network

M <- metaTagExtraction(M, Field = "AU_CO", sep = ";")
NetMatrix <- biblioNetwork(M, analysis = "collaboration", network = "countries", sep = ";")

# define functions from package Matrix
diag <- Matrix::diag 
colSums <-Matrix::colSums

# delete not linked vertices
ind <- which(Matrix::colSums(NetMatrix)-Matrix::diag(NetMatrix)>0)
NET <- NetMatrix[ind,ind]

# Select number of vertices to plot
n <- 20    # n. of vertices
NetDegree <- sort(diag(NET),decreasing=TRUE)[n]
NET <- NET[diag(NET)>=NetDegree,diag(NET)>=NetDegree]

# delete diagonal elements (self-loops)
diag(NET) <- 0

# Create igraph object
bsk.network <- graph.adjacency(NET,mode="undirected")

# Compute node degrees (#links) and use that to set node size:
deg <- degree(bsk.network, mode="all")
V(bsk.network)$size <- deg*1.1

# Remove loops
bsk.network <- simplify(bsk.network, remove.multiple = F, remove.loops = T) 

# Choose Network layout
#l <- layout.fruchterman.reingold(bsk.network)
l <- layout.circle(bsk.network)
#l <- layout.sphere(bsk.network)
#l <- layout.mds(bsk.network)
#l <- layout.kamada.kawai(bsk.network)


## Plot the network
plot(bsk.network,layout = l, vertex.label.dist = 0.5, vertex.frame.color = 'blue', vertex.label.color = 'black', vertex.label.font = 1, vertex.label = V(bsk.network)$name, vertex.label.cex = 0.5, main="Country collaboration")


Co-Citation Network

# Create a co-citation network

NetMatrix <- biblioNetwork(M, analysis = "co-citation", network = "references", sep = ".  ")

# define functions from package Matrix
diag <- Matrix::diag 
colSums <-Matrix::colSums

# delete not linked vertices
ind=which(Matrix::colSums(NetMatrix)-Matrix::diag(NetMatrix)>0)
NET=NetMatrix[ind,ind]

# Select number of vertices to plot
n <- 10    # n. of vertices
NetDegree <- sort(diag(NET),decreasing=TRUE)[n]
NET <- NET[diag(NET)>=NetDegree,diag(NET)>=NetDegree]

# delete diagonal elements (self-loops)
diag(NET) <- 0

# Create igraph object
bsk.network <- graph.adjacency(NET,mode="undirected")

# Remove loops
bsk.network <- simplify(bsk.network, remove.multiple = F, remove.loops = T) 

# Choose Network layout
l = layout.fruchterman.reingold(bsk.network)

## Plot
plot(bsk.network,layout = l, vertex.label.dist = 0.5, vertex.frame.color = 'blue', vertex.label.color = 'black', vertex.label.font = 1, vertex.label = V(bsk.network)$name, vertex.label.cex = 0.5, main="Co-citation network")


Keyword Coupling

# Create a co-citation network

NetMatrix <- biblioNetwork(M, analysis = "coupling", network = "keywords", sep = ";")

# define functions from package Matrix
diag <- Matrix::diag 
colSums <-Matrix::colSums

# delete not linked vertices
ind=which(Matrix::colSums(NetMatrix)-Matrix::diag(NetMatrix)>0)
NET=NetMatrix[ind,ind]

# Select number of vertices to plot
n <- 10    # n. of vertices
NetDegree <- sort(diag(NET),decreasing=TRUE)[n]
NET <- NET[diag(NET)>=NetDegree,diag(NET)>=NetDegree]

# delete diagonal elements (self-loops)
diag(NET) <- 0

# Plot Keywords' Heatmap (most frequent 30 words)
n=30
NETMAP=NetMatrix[ind,ind]
MapDegree <- sort(diag(NETMAP),decreasing=TRUE)[n]
NETMAP <- NETMAP[diag(NETMAP)>=MapDegree,diag(NETMAP)>=MapDegree]
diag(NETMAP) <- 0

H <- heatmap(max(NETMAP)-as.matrix(NETMAP),symm=T, cexRow=0.3,cexCol=0.3)


# Create igraph object
bsk.network <- graph.adjacency(NET,mode="undirected")

# Remove loops
bsk.network <- simplify(bsk.network, remove.multiple = T, remove.loops = T) 

# Choose Network layout
l = layout.fruchterman.reingold(bsk.network)


## Plot
plot(bsk.network,layout = l, vertex.label.dist = 0.5, vertex.frame.color = 'black', vertex.label.color = 'black', vertex.label.font = 1, vertex.label = V(bsk.network)$name, vertex.label.cex = 0.5, main="Keyword coupling")
