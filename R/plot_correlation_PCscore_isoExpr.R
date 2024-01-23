plot_correlation_PCscore_isoExpr <- function(res_isomas,iso,nfeat.pos = 50, nfeat.neg = 50,PCs="PC_2"){
  ###-----correlation between pc score and isoform expression
  library(ggplot2)
  data.heat = plot_isomas_heatmap(iso,dims = PCs, nfeat.pos = nfeat.pos, nfeat.neg = nfeat.neg,
                                 anno_row = F, anno_col = T, anno_col_keys = c("Mutation"),
                                 anno_legend = T, get.data = T,fs = 20,
                                 show_rwnms = T,
                                 iso_to_gene = NULL, show_gname = F, #remember not to show gene name
                                 title=NA)


  pcscore = data.heat$annotation_col[,PCs] #LUAD: Named num [1:517]
  names(pcscore) = rownames(data.heat$annotation_col)
  isoExpr = data.heat$scale.data.select #100 X 517

  if(all(colnames(isoExpr) %in% names(pcscore))){
    pcscore = pcscore[colnames(isoExpr)]
  }else{
    stop("names of PCscores and isoExpr not consistent.")
  }

  isoExpr = data.frame(t(isoExpr)) #rows: samples; columns:isoforms, 517 X 100
  res.cor.list = lapply(seq(ncol(isoExpr)),function(j){
    temp = cor.test(pcscore,isoExpr[,j],method = 'spearman')
    return(c(temp$estimate,temp$p.value))
  })
  names(res.cor.list) = colnames(isoExpr)
  res.cor.df = do.call(rbind,res.cor.list)
  colnames(res.cor.df) = c('Spearman.rho','P.value')
  res.cor.df = data.frame(res.cor.df) #100 X 2
  res.cor.df$Log10Pvalue = -log10(res.cor.df$P.value)
  res.cor.df$Color = 'black'
  res.cor.df$Color[res.cor.df$Spearman.rho>0] = "#FF0000"
  res.cor.df$Color[res.cor.df$Spearman.rho<0] = "#32CD32"

  df = res.cor.df

  fs = 24
  ggplot(df,aes(x=Spearman.rho,y=Log10Pvalue))+
    geom_point(color=df$Color,size=2)+
    labs(x=expression(paste('Spearman correlation (',rho,')')),
         y='-Log10 (P-value)')+
    theme(
      axis.text = element_text(size = fs),
      axis.title = element_text(size = fs),
      legend.position = 'none',
    )
}
