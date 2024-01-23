plot_PCscore_along_PCs <- function(res_isomas,data.maf,myGene.mut="TP53",PCs=1:5, method="wilcox",
                                   label.comp="p.format",fs=16,xlab="Principal Component",
                                   legend.direction = "horizontal",sz=6,
                                   legend.position = c(0.75,0.1),
                                   vars.comp = c("WildType","Mutant"),
                                   showNsamp = T,ynpc=0.95,
                                   title.legend=NULL,title.plot=NULL){
  # remove non-SNP samples, show samle size (Nsamp)
  # revised by Hua Tan on 6/27/2022 @13547

  library(dplyr)
  library(ggpubr)

  cell.embeddings = data.frame(res_isomas$pca@cell.embeddings)
  cell.embeddings$Mutation = get_maf_of_interest_TCGA(
    data.maf,
    myGene=myGene.mut,
    mySamples = rownames(cell.embeddings)
  )$Mutation
  # Mutant    Other WildType
  # 213       37      267

  if(is.numeric(PCs)){
    PCs = paste0("PC_",PCs)
  }

  if(showNsamp){
    Nsamp = table(cell.embeddings$Mutation)
    cell.embeddings$Mutation = paste0(cell.embeddings$Mutation,"(",Nsamp[cell.embeddings$Mutation],")")
    vars.comp = paste0(vars.comp,"(",Nsamp[vars.comp],")")
  }

  data.long = reshape2::melt(cell.embeddings,variable.name = "Variable", value.name="Value")
  # data.long = filter(data.long, Variable %in% paste0("PC_",c(1,2,4,6,7,11,14))) #for TP53
  # data.long = filter(data.long, Variable %in% paste0("PC_",c(1,4,6,11,21))) #for KEAP1
  data.long = filter(data.long, Variable %in% PCs,Mutation %in% vars.comp) #6/27/2022

  data.long$Mutation = factor(data.long$Mutation,
                              levels = vars.comp,
                              ordered = T)

  # data.long$Variable = as.character(data.long$Variable)

  # fs = 16
  ggboxplot(data.long, x = "Variable", y = "Value",
            color = "Mutation", palette = "jco",
            short.panel.labs = F, font.label = 20,
            add = "jitter",add.params = list(size=0.5))+
    stat_compare_means(aes(group = Mutation),
                       label = label.comp,
                       method = method,
                       label.x.npc = 0,
                       label.y.npc = ynpc,
                       size=sz)+
    labs(x=xlab, y="PC score",color=title.legend)+
    # labs(x=xlab, y="PCA-expression",color=title.legend)+
    grids(linetype = "dashed",color = "grey92")+
    ggtitle(title.plot)+
    theme(
      plot.title=element_text(size = fs+2,face='bold',hjust = 0),
      axis.text = element_text(size = fs),
      axis.title = element_text(size = fs+2),
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.text = element_text(size = fs),
      legend.title = element_text(size = fs+2)
    )
}
