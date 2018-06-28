######### All credits given to Dr. Gabriel Moreno-Hagelsieb for writing this script
######### This script will produce every file used for di-nucleotide
######### signatures, including transforming the distance file from
######### a "genome1, genome2, distance" format into a matrix proper,
######### both clusters (agglomerative and divisive) in Newick format
######### to view as a tree if so the user wants, asn files containing
######### groups of redundant genomnes at different thresholds of
######### dissimilarity.

library(cluster)
library(MCMCpack)
library(ape)
library(reshape2)

#### this is a function to provide names in each cluster at given threshold.
#### Input of the function are dendrogram (tree) and threshold of interest.
group.fun<-function(dendrogram,threshold) {
    g<-cutree(dendrogram,h=threshold)
    groupnumber<-length(table(g))
    namesingroup<-NULL
    b<-NULL
    for (i in 1:groupnumber) {
        namesingroup[i]<-list(rownames(signaturetbl)[g==i])
        b[i]<-noquote(paste("Group", i))
    }
    data.frame(noquote(I(unlist(lapply(namesingroup,paste,collapse=",")))),row.names=noquote(b))
}
### done with the function

############## import the distance matrix (Here your file Monica!)
dx<-read.table("DATA/neighbourhoods.tbl.bz2",header=T)

########## creating dendrograms
## create a dendrogram using the divisive (top down) approach
## and the dx dissimilarity matrix
divTree<-diana(dx,diss=T)
divDend<-as.hclust(divTree)

## Output the tree in the Newick format
phy <- as.phylo(as.hclust(divTree))
write.tree(phy, file="Clusters/neighbourhoods.dend")

## Produce groups of similar/redundant genomes, divisive, at
## several thresholds
## The code will give five seperate files that contains all the groups
## for the corresponding threshold.
## looping
for (i in 1:5){
    thrd<-i*0.1
    thrname = formatC(thrd,digits=2,format="f")
    groups<-group.fun(divDend,thrd)
    write.table(groups,file=paste(paste("Clusters/neighbourhoods",thrname,sep="-"),".tbl",sep=""),row.names=T,col.names=F,quote=F,append=F)
}

## calculate the divisive coefficients
coeff.div<-coef(divTree)
coeff.div
