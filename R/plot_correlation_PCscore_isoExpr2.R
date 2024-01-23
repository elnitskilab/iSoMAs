plot_correlation_PCscore_isoExpr2 <-
  function(res_isomas,iso,nfeat.pos = NULL, nfeat.neg = NULL,
           nfeat.pos.mark=50,nfeat.neg.mark=50,
           PCs="PC_2",fs=24){
  ###-----correlation between pc score and isoform expression
  # extend to all isoforms used in the PCA, but still only mark the top positive and negative isoforms
  # revised by Hua Tan on 5/9/2023 @4C08

  library(ggplot2)

  N = nrow(iso@reductions$pca@feature.loadings) #3315 for LUAD
  if(is.null(nfeat.pos)){
    nfeat.pos = N
  }
  if(is.null(nfeat.neg)){
    nfeat.neg = N
  }

  data.heat = plot_isomas_heatmap(iso,dims = PCs, nfeat.pos = nfeat.pos, nfeat.neg = nfeat.neg,
                                 anno_row = T, anno_col = T, anno_col_keys = c("Mutation"),
                                 anno_legend = T, get.data = T,fs = fs,
                                 show_rwnms = T,
                                 iso_to_gene = NULL, show_gname = F, #remember not to show gene name
                                 title=NA,plot=F)


  pcscore = iso@reductions$pca@cell.embeddings[,PCs]
  isoExpr = data.heat$scale.data.select
  pcscore = pcscore[colnames(isoExpr)]
  # print(all(names(pcscore)==colnames(isoExpr))) #TRUE
  isoExpr = data.frame(t(isoExpr)) #rows: samples; columns:isoforms
  res.cor.list = lapply(seq(ncol(isoExpr)),function(j){
    temp = cor.test(pcscore,isoExpr[,j],method = 'spearman')
    return(c(temp$estimate,temp$p.value))
  })
  names(res.cor.list) = colnames(isoExpr)
  res.cor.df = do.call(rbind,res.cor.list)
  colnames(res.cor.df) = c('Spearman.rho','P.value')
  res.cor.df = data.frame(res.cor.df)
  res.cor.df$Log10Pvalue = -log10(res.cor.df$P.value)
  res.cor.df$PC_direction = data.heat$annotation_row[rownames(res.cor.df),'PC_direction']
  idx.pos = which(res.cor.df$PC_direction == 'Positive')
  idx.neg = which(res.cor.df$PC_direction == 'Negative')
  res.cor.df$Color = 'gray50'
  res.cor.df$Color[idx.pos[1:nfeat.pos.mark]] = "#FF0000"
  res.cor.df$Color[idx.neg[1:nfeat.neg.mark]] = "#32CD32"
  res.cor.df$Shape = 1
  res.cor.df$Shape[idx.pos[1:nfeat.pos.mark]] = 24
  res.cor.df$Shape[idx.neg[1:nfeat.neg.mark]] = 25
  res.cor.df$Size = 1
  res.cor.df$Size[idx.pos[1:nfeat.pos.mark]] = 3
  res.cor.df$Size[idx.neg[1:nfeat.neg.mark]] = 3

  df = res.cor.df
  ggplot(df,aes(x=Spearman.rho,y=Log10Pvalue))+
    geom_point(color=df$Color,size=df$Size,shape=df$Shape)+
    labs(x=expression(paste('Spearman correlation (',rho,')')),
         y='-Log10 (P-value)')+
    theme(
      axis.text = element_text(size = fs),
      axis.title = element_text(size = fs),
      legend.position = 'none',
    )
}
