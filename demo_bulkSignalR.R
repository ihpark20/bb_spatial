
library(BulkSignalR)
library(igraph)
library(dplyr)
library(doParallel) # for parallel processing

n.proc <- 4 # number of cores to use
cl <- makeCluster(n.proc)
registerDoParallel(cl)
#stopCluster(cl)

data(sdc)
head(sdc)

bsrdm <- prepareDataset(counts=sdc)
bsrdm <- learnParameters(bsrdm, plot.folder = "demo/bulkSignalR/figures", filename="sdc")

#?su
#q?learnParameters
bsrinf <- initialInference(bsrdm)
LRinter.dataframe <- LRinter(bsrinf)

head(LRinter.dataframe)

LRinter.dataframe.selected <- LRinter.dataframe[LRinter.dataframe$qval<= 0.001,]
LRinter.dataframe.selected <- LRinter.dataframe.selected[order(LRinter.dataframe.selected$qval),]

bsrinf.redBP <- reduceToBestPathway(bsrinf)
bsrinf.L    <- reduceToLigand(bsrinf)
bsrinf.R    <- reduceToReceptor(bsrinf)

bsrinf.redP  <- reduceToPathway(bsrinf)  
bsrinf.redPBP <- reduceToBestPathway(bsrinf.redP) 

bsrsig.redBP <- getLRGeneSignatures(bsrinf.redBP, qval.thres=0.001)

scoresLR <- scoreLRGeneSignatures(bsrdm, bsrsig.redBP, name.by.pathway=FALSE)

simpleHeatmap(scoresLR[1:20,], path="demo/bulkSignalR/figures/", filename="sdc_scoresLR", column.names=TRUE, height=5, width=9, pointsize=10, hcl.palette="Cividis")

bsrsig.redPBP <- getLRGeneSignatures(bsrinf.redPBP, qval.thres=0.001)
     
scoresPathway <- scoreLRGeneSignatures(bsrdm, bsrsig.redPBP, name.by.pathway=TRUE)
     
simpleHeatmap(scoresPathway[1:10,], 
                   path = "demo/bulkSignalR/figures/",
                   filename = "sdc_scoresPathway",
                   column.names = TRUE, 
                   width = 9, 
                   height = 5, 
                   pointsize = 12,
                   hcl.palette = "Blue-Red 2"
                   )

pathway1 <- "Elastic fibre formation"
signatureHeatmaps(
        pathway = pathway1,
        bsrdm = bsrdm,
        bsrsig = bsrsig.redPBP,
        path = "demo/bulkSignalR/figures/",
        filename = "sdc_signatureheatmap",
        h.width  = 15,
        h.height = 10 ,
        show_column_names = TRUE)

        
alluvialPlot(bsrinf,
              keywords = c("COL4A1"),
              type = "L",
              qval.thres = 0.001,
              path = "demo/bulkSignalR/figures/",
              filename = "sdc_alluvial", 
              width  = 16, 
              height = 12
              )

pathway1 <- "PD-1 signaling"
pathway2 <- "Interferon gamma signaling"
pathways <- c(pathway1,pathway2)

bubblePlotPathwaysLR(bsrinf,
    pathways = pathways, 
    qval.thres  = 1,
    path = "demo/bulkSignalR/figures/",
    color = "red",
    filename  = "sdc_bubble", 
    width  = 16, 
    height = 7,
    pointsize = 8
    #filter.L = c("ADAM12")
    #filter.R = c("ITGA3")
    ) 

chordDiagramLR (bsrinf,
                  path = "demo/bulkSignalR/figures/",
                  filename = "sdc_chord",
                  pw.id.filter = "R-HSA-202733",
                  limit = 20,
                  width = 5, 
                  height = 4.5
    )


gLR <- getLRNetwork(bsrinf.redBP, qval.thres=1e-8)
# As an alternative to Cytoscape, you can play with igraph package functions.
plot(gLR,
     edge.arrow.size=0.1,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.1)
dev.off()


u.gLR <- as.undirected(gLR) # most algorithms work for undirected graphs only
comm <- cluster_edge_betweenness(u.gLR)
cb <- cohesive_blocks(u.gLR)
plot(cb,u.gLR,
     vertex.label.color="black",
     vertex.label.family="Helvetica",
     vertex.label.cex=0.1,
     edge.color="black")
dev.off()
