run_PCA_Seurat <- function(data.iso,res.cluster = 0.5, scale.factor = 10000,
                           normalization.method = "LogNormalize",
                           selection.method = "mvp",
                           min.cells=25, min.features=5986,
                           nPCs=50,project="TCGA-cancer",res_isomas=NULL){
  library(Seurat)

  if(!is.null(res_isomas)){
    res.cluster = res_isomas$res.cluster
    scale.factor = res_isomas$scale.factor
    normalization.method = res_isomas$normalization.method
    selection.method = res_isomas$selection.method
    min.cells = res_isomas$min.cells
    min.features = res_isomas$min.features
    nPCs = res_isomas$nPCs
    project = res_isomas$project
  }

  iso <- CreateSeuratObject(counts = as.matrix(data.iso), project = project,
                            min.cells = min.cells, min.features = min.features)

  cat(paste0("dim of iso: ", paste(dim(iso),collapse = ", "),"\n"))
  # dim of iso: 53963, 517

  iso <- NormalizeData(iso, normalization.method = normalization.method, scale.factor = scale.factor)
  iso <- FindVariableFeatures(iso, selection.method = selection.method)
  varFeatures = VariableFeatures(iso) #3315
  cat(paste0("selection.method of mean.var.plot (mvp) for FindVariableFeatures was used,
             length of VariableFeatures(iso): ",length(varFeatures),"\n"))
  iso <- ScaleData(iso, features = rownames(iso))

  ##--Run PCA
  if(ncol(iso) <= nPCs){
    cat(paste0("ncol(iso)=",ncol(iso),", nPCs=",nPCs,"\n"))
    cat("#columns of iso is smaller than #components, downsizing nPCs to ncol(iso)/3\n")
    nPCs = max(floor(ncol(iso)/3),2)
  }
  iso <- RunPCA(iso, features = VariableFeatures(object = iso),npcs = nPCs,nfeatures.print = 6)
  iso <- FindNeighbors(iso, dims = seq(nPCs))
  iso <- FindClusters(iso, resolution = res.cluster)
  return(iso)
}

