plot_isomas_heatmap <- function(iso, dims=2, nfeat.pos=12, nfeat.neg=12,
                               anno_row=T, anno_col=T, anno_legend=F,
                               anno_col_keys = c("TUMOR_STAGE","seurat_clusters","Mutation"),
                               title=NULL, get.data=F, lgd = T,trimLarge=T,
                               show_rwnms=T, show_clnms=F,show_gname=F,
                               iso_to_gene = NULL,fs_row=16,fs_col=16,
                               fs=16,plot=T){

  if(is.numeric(dims)){
    PCs = paste0("PC_",dims)
  }else if(!grepl("^PC_",dims)){
    stop("a character dims should start with 'PC_'\n")
  }else{
    PCs = dims
  }

  scale.data = as.data.frame(iso@assays$RNA@scale.data)
  feature.loadings = as.data.frame(iso@reductions$pca@feature.loadings)

  ##--
  feature.pos = data.frame(feature.loadings[feature.loadings[,PCs]>0,PCs,drop=F])
  feature.pos = feature.pos[order(feature.pos[,PCs],decreasing = T),,drop=F]
  feature.neg = data.frame(feature.loadings[feature.loadings[,PCs]<0,PCs,drop=F])
  feature.neg = feature.neg[order(feature.neg[,PCs],decreasing = F),,drop=F]

  if(!is.null(nfeat.pos)){
    nfeat.pos = min(nrow(feature.pos),nfeat.pos) #3/21/2023 @4C08
  }else{
    nfeat.pos = nrow(feature.pos)
  }

  if(!is.null(nfeat.neg)){
    nfeat.neg = min(nrow(feature.neg),nfeat.neg)
  }else{
    nfeat.neg = nrow(feature.neg)
  }

  feature.select = c(rownames(feature.pos)[1:nfeat.pos],rownames(feature.neg)[1:nfeat.neg])

  if(anno_row){
    annotation_row = data.frame(PC_direction=c(rep("Positive",nfeat.pos),
                                               rep("Negative",nfeat.neg)))
    annotation_row$PC_direction = factor(annotation_row$PC_direction,
                                         levels = c("Positive","Negative"),
                                         ordered = T)
    rownames(annotation_row) = feature.select
  }else{
    annotation_row = NA
  }
  ##--

  cell.embeddings = as.data.frame(iso@reductions$pca@cell.embeddings)
  cell.PCs = cell.embeddings[,PCs,drop=F]
  cell.PCs = cell.PCs[order(cell.PCs[,PCs]),PCs,drop=F]

  scale.data.select = scale.data[feature.select,rownames(cell.PCs)]
  # scale.data.select[abs(scale.data.select)>2] = 2 #this was wrong
  scale.data.select[scale.data.select > 2] = 2
  scale.data.select[scale.data.select < -2] = -2 #by HT on 5/3/22
  #
  if(anno_col){
    annotation_col = iso@meta.data[colnames(scale.data.select), anno_col_keys,drop=F]
    # colnames(annotation_col)[2] = "Clusters"
    annotation_col[,PCs] = cell.PCs[rownames(annotation_col),PCs]
  }else{
    annotation_col = NA
  }

  if (show_gname && !is.null(iso_to_gene)){ #2/18/22 @13547
    rownames(scale.data.select) = paste0(rownames(scale.data.select),"(",
                                         iso_to_gene[rownames(scale.data.select)],")")
  }else{
    warning("if gene names not shown, check iso_to_gene availability")
  }

  if(plot){
    pheatmap(scale.data.select,cluster_rows = F,cluster_cols = F,
             show_rownames = show_rwnms, show_colnames = show_clnms,scale = "row",
             main = title,
             annotation_col = annotation_col,
             annotation_row = annotation_row,
             annotation_names_row = F,
             annotation_legend = anno_legend,fontsize = fs,
             legend = lgd, fontsize_row = fs_row,fontsize_col = fs_col,
             legend_breaks = seq(-2,2), legend_labels = c("<-2","-1","0","1",">2"),
             drop_levels = T,color = colorRampPalette(c("magenta","black","yellow"))(10))
  }

  if(get.data){
    return(list(scale.data.select=scale.data.select,
                annotation_row=annotation_row,
                annotation_col=annotation_col,
                feature.select=feature.select))
  }

}
