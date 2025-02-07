plot_feature_loadings_along_PCs <- function(res_isomas,PCs="PC_2",PC.ord=NULL,fs=16){
  library(ggplot2)

  if(is.null(PC.ord)){
    PC.ord = PCs[1]
    cat(paste0('PC.ord not designated, set to the first PC-axis in PCs: ',PC.ord,'\n'))
  }
  feature.loadings = res_isomas$pca@feature.loadings[,PCs,drop=F]
  feature.loadings.ord = data.frame(feature.loadings[order(feature.loadings[,PC.ord],decreasing = T),,drop=F])
  feature.loadings.ord$order = seq(nrow(feature.loadings.ord))
  data.long = reshape2::melt(feature.loadings.ord,
                             id.vars = "order",
                             variable.name="PC",
                             value.name="Loading")

  ggplot(data.long,aes(x=order,y=Loading,color=PC))+
    geom_point(shape=21,size=2)+
    geom_smooth(size=2)+
    labs(x="Feature Order",y="PC Loading")+
    theme(
      axis.text = element_text(size = fs),
      axis.title = element_text(size = fs+2),
      legend.position = c(0.5,0.05),
      legend.text = element_text(size = fs),
      legend.background = element_blank(),
      legend.direction = 'horizontal',
      legend.title = element_blank()
    )
}
