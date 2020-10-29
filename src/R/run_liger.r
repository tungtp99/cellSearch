.libPaths(c('/home/tung/.BioTBDataDev/App/dependencies', '/home/tung/.BioTBDataDev/App/analytics', '/home/tung/R/x86_64-pc-linux-gnu-library/3.6'))
library(R.utils)
library(Matrix)
library(liger)


list.file <- dir()
list.not.study <- grep('[.]', list.file, value=TRUE)
list.study <- setdiff(list.file, list.not.study)

#c('GSE135922_all', 'GSE123904_human', 'GSE140228_10X', 'GSE133549_human_hamony', 'GSE135922_EC', 'GSE133549_human_reference', 'Sharma_14HCC', 'E-MTAB-6149_all', 'GSE131907', 'GSE132465', 'SDY997_landscape', 'GSE145281_RN', 'madissoon2020_esophagus', 'GSE146771_10x', 'GSE144469_CD3', 'GSE146771_smartseq', 'madissoon2020_lung', 'GSE124310', 'GSE150430', 'GSE138266', 'GSE123814_BCC', 'young2018', 'GSE114725', 'GSE128223', 'james2020', 'GSE139324', 'GSE144735', 'GSE140819_Melanoma', 'GSE139555', 'GSE99254', 'stewart2019_mature', 'GSE143363', 'GSE126030', 'E-MTAB-6149_t', 'roider2019', 'GSE145281_BT', 'GSE115978', 'AN1801', 'stewart2019_fetal', 'GSE146811_human', 'GSE110686')
#list.study <- c('GSE135922_all', 'GSE123904_human', 'GSE140228_10X')#, 'GSE133549_human_hamony', 'GSE135922_EC', 'GSE133549_human_reference', 'Sharma_14HCC', 'E-MTAB-6149_all', 'GSE131907', 'GSE132465', 'SDY997_landscape', 'GSE145281_RN', 'madissoon2020_esophagus', 'GSE146771_10x', 'GSE144469_CD3', 'GSE146771_smartseq', 'madissoon2020_lung', 'GSE124310', 'GSE150430', 'GSE138266', 'GSE123814_BCC', 'young2018', 'GSE114725', 'GSE128223', 'james2020', 'GSE139324', 'GSE144735', 'GSE140819_Melanoma', 'GSE139555', 'GSE99254', 'stewart2019_mature', 'GSE143363', 'GSE126030', 'E-MTAB-6149_t', 'roider2019', 'GSE145281_BT', 'GSE115978', 'AN1801', 'stewart2019_fetal', 'GSE146811_human', 'GSE110686')
list.matrix <- list()


for (i in 1:length(list.study)) {
  tmp <- Signac::ReadSpMt(file.path('.', list.study[i], 'matrix.hdf5'), '/')
  barcodes <- colnames(tmp)
  barcodes.cd4 <- grep('00220', barcodes, value=TRUE)
  pos.cd4 <- match(barcodes.cd4, barcodes)
  #tmp <- tmp[,pos.cd4]
  #if (length(pos.cd4)  > 0) {
  list.matrix[list.study[i]] <- tmp
  #}
}

#list.matrix <- list.matrix[c(TRUE, TRUE, FALSE, TRUE, TRUE)]


ligerex = createLiger(list.matrix) #Can also pass in more than 2 datasets
ligerex = normalize(ligerex)
ligerex = selectGenes(ligerex, var.thresh = 0.1)
ligerex = scaleNotCenter(ligerex)
ligerex = optimizeALS(ligerex, k = 20, num.cores = 4) 
ligerex = quantileAlignSNF(ligerex) 
ligerex = runTSNE(ligerex)
ligerex = runUMAP(ligerex)
plotByDatasetAndCluster(ligerex) #Can also pass in different set of cluster labels to plot
plotFeature(ligerex, "nUMI")
pdf("word_clouds.pdf")
plotWordClouds(ligerex)
dev.off()
pdf("gene_loadings.pdf")
plotGeneLoadings(ligerex)
dev.off()




for (i in 1:length(list.matrix)) {
  Signac::WriteSpMtAsSpMat(file=file.path('../220',paste(names(list.matrix)[i],'.hdf5', sep ='')), group='/', mat=list.matrix[[i]])
}


##['GSE145281', 'GSE135922', 'GSE143363', 'roider2019', 'GSE128223', 'GSE140228', 'GSE139324']
##['stewart2019', 'GSE126030', 'young2018', 'GSE144735', 'madissoon2020', 'roider2019', 'E-MTAB-6149', 'GSE128223', 'GSE140228', 'GSE139324']