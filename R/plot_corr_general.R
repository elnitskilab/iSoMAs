plot_corr_general <- function(vec1,vec2,legend.pos=c(0.1,0.9),
                              Log2 = T,title=NULL, smoothPoint=T, addPval=T,
                              fs = 16,sz=6,xlab=NULL,ylab=NULL,xylim=T,
                              rm0=F){
  # plot correlation between two vectors with same length
  library(dplyr)
  library(ggplot2)
  library(ggpmisc)
  library(ggrepel)

  data = data.frame(g1=vec1,g2=vec2)

  if(rm0){
    idx_kp = data$g1!=0 & data$g2!=0
    data = data[idx_kp,,drop=F]
    if(nrow(data)<3){
      cat("samples less than 3 after removing 0's\n")
      return(NULL)
    }
  }

  if(Log2){
    data = log2(1+data)
  }
  # print(head(data))
  # fs = 14
  p <- ggplot(data,aes(x=g1,y=g2))+
    geom_point(color = "blue",size = 2, alpha = 0.8)+
    # labs(x=gene1,y=gene2)+
    labs(x=xlab,y=ylab)+
    ggtitle(title)+
    theme(
      axis.text = element_text(size = fs),
      axis.title = element_text(size = fs+2),
      plot.title = element_text(size = fs+4),
      legend.position = legend.pos,
      panel.background = element_rect(fill = "white",
                                      colour = "black",
                                      size = 0.5,
                                      linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'dashed',
                                      colour = "gray75"),
      panel.grid.minor = element_line(size = 0.25, linetype = 'dashed',
                                      colour = "gray75")
    )
  if(xylim){
    xyMin = min(data,na.rm = T)
    xyMax = max(data,na.rm = T)
    p=p+lims(x=c(xyMin,xyMax),y=c(xyMin,xyMax))
  }
  if(smoothPoint){
    p = p+geom_smooth(data=data,method = lm,se = T,formula = y ~ x,color='black')
  }

  if(addPval){
    stats = cor.test(data$g1,data$g2,method = "kendall")
    R = round(stats$estimate,2)
    P = signif(stats$p.value,2)
    n = nrow(data)
    anno = paste0("R=",R, "\nP=",P,"\nn=",n)
    # xma = max(data$g1,na.rm = T)
    # yma = max(data$g2,na.rm = T)
    # p = p+geom_text(label=anno,x=xyMax,y=xyMin,hjust=1,vjust=0,color="red",size=6)
    if(xylim){
      p = p+geom_text(label=anno,
                      x=xyMin,
                      y=xyMax,
                      hjust=0,vjust=1,color="red",size=sz)
    }else{
      p = p+geom_text(label=anno,
                      x=min(data$g1,na.rm = T),
                      y=max(data$g2,na.rm = T),
                      hjust=0,vjust=1,color="black",size=sz)
    }
  }
  return(p)
}
