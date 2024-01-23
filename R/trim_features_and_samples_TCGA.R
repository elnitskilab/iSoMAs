trim_features_and_samples_TCGA <- function(data.iso,gene_to_iso=NULL,iso_to_gene=NULL){
  ##--trim features: keep multiple-isoform genes only
  if(!is.null(gene_to_iso) & !is.null(iso_to_gene)){
    genes_list = names(gene_to_iso) #29,181
    genes_list = sort(genes_list[unlist(lapply(genes_list,function(x){length(gene_to_iso[[x]])>1}))]) #15,448
    idx_kp = iso_to_gene[rownames(data.iso)] %in% genes_list #59,866 out of 73,599
  }else{
    idx_kp = seq(nrow(data.iso))
    cat("gene_to_iso or iso_to_gene is NULL, all features kept\n")
  }

  ##--trim samples: remove normal samples based on TCGA barcode
  sampTypeCode = as.numeric(substr(colnames(data.iso),14,15))
  data.iso = data.iso[idx_kp,sampTypeCode<10,drop=F] #59866, 517

  # data.iso is a matrix-like object with unnormalized data with cells as columns and features as rows
  cat(paste0("dim of trimmed data.iso: ", paste(dim(data.iso),collapse = ", "),"\n"))
  return(data.iso)
}

