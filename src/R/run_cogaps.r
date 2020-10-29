
params <- new("CogapsParams")
params <- setParam(params, "nPatterns", 100)

list.study <- c('SDY997_landscape', 'GSE148837', '56535e196ac6-5c3e-80e7-f408e40e246a')
list.matrix <- lapply(list.study, function(study) {
  Signac::ReadSpMt(file.path('/home/tung/RepresentData', study, 'matrix.hdf5'), '/')
})
names(list.matrix) <- list.study

list.genes <- rownames(matrix)
for (matrix in list.matrix) {
  list.genes <- intersect(list.genes, rownames(matrix))
}


for (i in 1:length(list.matrix)) {
  matrix <- list.matrix[i]
  pos.rows <- match(list.genes, rownames(matrix))
  list.matrix[i] <- matrix[pos.row,]
}

P <- NaN
if (length(list.genes) < 100) {
  print("List genes smaller than 100 genes")
} else {
  if (P != NaN) { 
    cogap <- CoGAPS(data = matrix, distributed="single-cell")
  } else {
    cogap <- CoGAPS(data = matrix, distributed="single-cell")
  }
}
