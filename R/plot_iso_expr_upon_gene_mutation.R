plot_iso_expr_upon_gene_mutation <- function(data.iso,Mutations,gene_to_iso, myGene.mut="TP53",
                                             myGene.iso="TPX2", isoIDs=NULL,Log2=T,normData=T,
                                             stat.method = "wilcox",label = "p.signif",
                                             plot.percent = T,get.data=F, title=NULL,
                                             legend.position=c(0.75,0.1),
                                             vars.comp = c("WildType","Mutant"),
                                             legend.direction = "horizontal",label.y.npc=0.975,
                                             plot.gene=T,rmNorm=F,angle.p=90,vjust.p=0.75,
                                             showNsamp = T,fs = 12,sz=8){

  library(ggpubr)
  library(dplyr)

  # method = {t.test,wilcox.test,anova,kruskal.test}
  if(normData){
    cat("normalize total counts to 1,000,000 for each sample...\n")
    colsum = colSums(data.iso)
    for (j in seq(ncol(data.iso))){
      data.iso[,j] = data.iso[,j]/colsum[j]*1e6
    }
  }

  if(is.null(isoIDs)){  #4/22/22
    isoIDs = gene_to_iso[[myGene.iso]]
  }

  data.iso = as.data.frame(t(data.iso[isoIDs,])) #LUAD: 576 * 13

  if(rmNorm){
    samples.all = row.names(data.iso)
    sample.type.code = as.numeric(substr(samples.all,14,15))
    data.iso = data.iso[sample.type.code<10,,drop=F]
    cat(paste0('after removing normals: ',paste(dim(data.iso),collapse =', '),'\n'))
  }

  if(plot.percent){
    # percent should be calculated on count not log2count
    data.iso = round(100*data.iso/(rowSums(data.iso,na.rm = T)+1e-6),2)
    value.name = paste0(myGene.iso," (%)")
  }else if(Log2){
    data.iso = log2(1+data.iso)
    valueName = "Log2FPKM"
  }else{
    valueName = 'FPKM'
  }

  # ymax = (max(data.iso) - min(data.iso))*1.1

  ##--attach mutant information
  if(showNsamp){
    Mutations = Mutations[rownames(data.iso),,drop=F] #7/12/22
    Nsamp = table(Mutations$Mutation)
    Mutations$Mutation = paste0(Mutations$Mutation,"(",Nsamp[Mutations$Mutation],")")
    vars.comp = paste0(vars.comp,"(",Nsamp[vars.comp],")")
  }

  if(plot.gene){
    data.iso[,myGene.iso] = rowSums(data.iso,na.rm = T,dims = 1)
  }

  data.iso$Mutation = Mutations[rownames(data.iso),'Mutation']

  data.iso.long = reshape2::melt(data.iso,
                                 id.vars = "Mutation",
                                 variable.name = "Isoform",
                                 value.name=valueName)
  data.iso.long = filter(data.iso.long, Mutation %in% vars.comp)

  data.iso.long$Mutation = factor(data.iso.long$Mutation,
                                  levels = vars.comp,
                                  ordered = T)

  print(str(data.iso.long))

  p = ggboxplot(data.iso.long, x = "Isoform", y = valueName,
                color = "Mutation", palette = "jco",
                short.panel.labs = F, font.label = 20,
                add = "jitter",add.params = list(size=0.5))+
    stat_compare_means(aes(group = Mutation),
                       label = label,
                       method = stat.method,
                       angle=angle.p,
                       vjust=vjust.p,
                       label.x.npc = 0,
                       label.y.npc = label.y.npc,
                       size=sz)+
    grids(linetype = "dashed",color = "grey92")+
    labs(x=myGene.iso,color=myGene.mut)+
    # lims(y=c(0,ymax))+
    ggtitle(title)+
    theme(
      plot.title=element_text(size = fs+2,face='bold',hjust = 0),
      axis.text.x = element_text(size = fs),
      # axis.text.x = element_text(size = fs,angle = 45,hjust = 1),
      # axis.text.x = element_text(size = fs,angle = 90,vjust = 0.5),
      axis.text.y = element_text(size = fs),
      axis.title = element_text(size = fs+2),
      legend.position = legend.position,
      legend.direction = legend.direction,
      legend.text = element_text(size = fs),
      legend.title = element_text(size = fs+2)
    )

  print(p)

  if(get.data){
    return(data.iso.long)
  }
}
