plot_isomas_pvals <- function(pvals_all,sorted = T,pval.thres = 1e-3,
                             showNsamp = T,idx.plot = 50, title = NULL,
                             p.log10 = T, lgd =T,jPC1=4,
                             show_rwnms=T, show_clnms=T, plot.sig=F,
                             get.data=F,scale.by="none",fs=12){
  # revised by Hua Tan on 11/23/2022 @4C08

  library(pheatmap)

  if(sorted){
    minP = apply(pvals_all[,jPC1:ncol(pvals_all)],1,min)
    pvals_all = pvals_all[order(minP),]
  }
  if(!is.null(pval.thres)){
    idx_sig = apply(pvals_all[,jPC1:ncol(pvals_all)],1,function(x){min(x) < pval.thres})
    pvals_all = pvals_all[idx_sig,,drop=F]
    cat('significant p-values: ',nrow(pvals_all),'\n')
  }

  minP = apply(pvals_all[,jPC1:ncol(pvals_all)],1,min)
  minP.idx = apply(pvals_all[,jPC1:ncol(pvals_all)],1,which.min)
  # print(table(minP.idx))
  print(head(minP))

  if(showNsamp){
    rownames(pvals_all) = paste0(rownames(pvals_all),
                                 "(",pvals_all$Mutant,"/",
                                 pvals_all$WildType,", ",
                                 minP,")")
  }

  if(length(idx.plot)==1){
    idx.plot = seq(idx.plot)
  }

  pvals_all.plot = as.matrix(pvals_all[idx.plot,jPC1:ncol(pvals_all)])

  mat.plot = pvals_all.plot

  # print(dim(pvals_all))
  nc = ncol(mat.plot)
  print(head(mat.plot[,c(1:2,(nc-1):nc)]))

  if(plot.sig){
    cat("plot significant only\n")
    idx_r = apply(mat.plot,1,function(x){!any(is.na(x))})
    mat.plot = mat.plot[idx_r,,drop=F]
    idx_r2 = apply(mat.plot,1,function(x){min(x,na.rm = T)<pval.thres})
    mat.plot = mat.plot[idx_r2,,drop=F]
    idx_c = apply(mat.plot,2,function(x){min(x,na.rm = T)<pval.thres})
    mat.plot = mat.plot[,idx_c,drop=F]

    nc = ncol(mat.plot)
    print(head(mat.plot[,c(1:2,(nc-1):nc)]))

    if(sum(idx_c)<2){
      cat('significant PCs are located in less than two PC axes, plotting significant with pheatmap is impossible,
          switch back to pvals_all.plot\n')
      mat.plot = pvals_all.plot
      print(dim(mat.plot))
    }
  }

  if(p.log10){
    cat("plot -log10(p-values)\n")
    pval.thres.log10 = ceiling(-log10(pval.thres))
    mat.plot = -log10(mat.plot)
    mat.plot[mat.plot>pval.thres.log10] = pval.thres.log10
    pheatmap(mat.plot,cluster_rows = F, cluster_cols = F,scale = scale.by,main = title,
             show_rownames = show_rwnms, show_colnames = show_clnms,
             fontsize = fs,legend =lgd,
             legend_breaks = seq(0,pval.thres.log10),
             legend_labels = c(as.character(seq(0,pval.thres.log10-1)),paste0(">=",pval.thres.log10)))
  }else{
    cat('plot original p-values\n')
    pheatmap(mat.plot,cluster_rows = F, cluster_cols = F,scale = scale.by,main = title)
  }
  if(get.data){
    return(mat.plot)
  }
}
