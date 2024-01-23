get_pvals_sig <- function(pvals_all,myGene='TP53',jPC1=4,pval.thres=1e-3){
  idx = which(rownames(pvals_all)==myGene)
  PCs_sig = which(pvals_all[idx,jPC1:ncol(pvals_all)] < pval.thres)
  PCs = paste0("PC_",PCs_sig)
  pvals_sig = pvals_all[myGene,PCs,drop=F]
  print(pvals_sig)
  return(pvals_sig)
}
