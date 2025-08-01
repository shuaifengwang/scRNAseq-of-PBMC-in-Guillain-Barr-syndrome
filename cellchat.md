# cellchat

library(CellChat)
data.input <- GetAssayData(sce2, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(sce2)
meta <- data.frame(group = labels, row.names = names(labels))

cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")

CellChatDB <- CellChatDB.human
 
showDatabaseCategory(CellChatDB)
dir.create("../5-cellchat");setwd("../5-cellchat")
ggsave("CellChatDB.pdf", width = 8, height = 6)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use # set the used database in the object
unique(CellChatDB$interaction$annotation)

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)  # 默认 raw.use = TRUE

cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
head(df.netp)
write.csv(df.netp, "net_pathway.csv")

cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

 netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
 netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


pathways.show <- "MIF" #Cxcl Tnf
# netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)   
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)


levels(cellchat@idents)
netVisual_bubble(cellchat, sources.use = CD8+ T cell, targets.use = NULL, remove.isolate = FALSE)
ggsave("L-R pairs.pdf", width = 6, height = 6)


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 15, height = 6, font.size = 10)
# # save as TIL/SNA_CXCL_signalingRole.pdf

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 5)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 5)
ht1 + ht2
ggsave("cellchat_SignalingPattern.pdf", width = 8, height = 5)
# save as TIL/SNA_SignalingPattern.pdf

plotGeneExpression(cellchat, signaling = "MIF")
ggsave("GeneExpression_MIF.pdf", width = 6, height = 4)
plotGeneExpression(cellchat, signaling = "Tnf")
ggsave("GeneExpression_Tnf.pdf", width = 6, height = 4)
plotGeneExpression(cellchat, signaling = "Cxcl")
ggsave("GeneExpression_Cxcl.pdf", width = 6, height = 4)
