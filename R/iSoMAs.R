#' iSoMAs: iSoform expression and somatic Mutation Association
#'
#' \code{iSoMAs} is an efficient computational pipeline based on principal component analysis (PCA) technique
#' for exploring the role of somatic mutations in shaping the landscape of gene isoform expression at the
#' transcriptome level.
#'
#' iSoMAs integrates sample-matched RNA-seq and DNA-seq data to study the association between gene somatic
#' mutation and gene isoform expression via both cis- and trans-acting mechanisms.
#' The iSoMAs workflow consists of two steps:
#' In the first step, the high-dimensional gene isoform expression matrix (d=59,866 derived from the 15,448
#' multi-isoform genes) for each cancer type is trimmed by mean.var.plot method (built in the Seurat toolkit)
#' into a more informative expression matrix, which keeps only the most variable isoforms but is still
#' high-dimensional (d~3,500). After that, the PCA is performed to further reduce the dimension of the
#' informative expression matrix into a much lower-dimensional PC score matrix (d<=50) by calculating a PC
#' loading matrix. Each column of the PC loading matrix performs a particular linear combination of the top
#' variable isoforms into a meta-isoform, with the combination coefficients stored in the corresponding column
#' of the PC loading matrix. All the meta-isoforms comprise the coordinates of the new low-dimensional space.
#' In the second step, a differential PC score analysis is conducted along each of the 50 PC-coordinates based
#' on the mutation status of the studied gene by Wilcoxon rank-sum test. Following the differential PC score analysis,
#' the significant genes (candidate iSoMAs genes) are identified if the minimum of the 50 p-values [defined as
#' `minP = min{P1,P2,...,P50}`] is smaller than a predefined threshold (i.e.,`minP<1e-3`). The iSoMAs genes are
#' eventually determined after a two-step multiple testing correction procedure on all candidate iSoMAs genes.
#'
#' @param data.iso Gene expression matrix at the isoform-level. Each row represents a transcript (isoform) and
#' each column represents a sample.
#' @param data.maf Gene somatic mutation tab in mutation annotation format (maf)
#' @param gene_to_iso A list of list that maps each gene to its isoforms (transcripts)
#' @param iso_to_gene A list mapping each isoform to its gene
#' @param genes_test A list of genes to test. If NULL, iSoMAs will automatically determines a qualified list of genes
#' for test based on the input data.maf, mut.num.thres and mut.percent.thres.
#' @param stat.method Statistical method used for the differential PC score analysis, choose one of:
#'   \itemize{
#'   \item wilcox for Wilcoxon's rank-sum test.
#'   \item t.test for Student's t-test.
#'   }
#' @param mut.num.thres Threshold number of samples for determining genes_test.
#' @param mut.percent.thres Threshold percent of samples for determing genes_test. iSoMAs will take the maximum of
#' mut.num.thres and mut.percent.thres to determine the final mut.samp.thres.
#' @param varType Designate the specific Variant_Type of somatic mutation in data.maf.
#' @param varClassification Designate the specific Variant_Classification of somatic mutation in data.maf.
#' @param res.cluster Argument used in \code{Seurat::FindClusters}.
#' @param scale.factor Argument used in \code{Seurat::NormalizeData}.
#' @param normalization.method Argument used in \code{Seurat::NormalizeData}.
#' @param selection.method Argument used in \code{Seurat::FindVariableFeatures}.
#' @param nPCs Number of PCs in PCA.
#' @param return.Seurat Return the Seurat object.
#' @minP.sorted Sort genes based on minP across the nPCs coordinates.
#' @param project Name of the project.
#' @param filename Name of file for saving the iSoMAs results.
#' @ntop.gene Number of top genes from genes_test for test.
#'
#' @return This function returns a list with the following 20 elements:
#' pca, pvals_all, pvals_all_sorted, genes_test, min.cells, min.features, dim_data.iso, dim_iso, project,
#' mut.samp.thres, varType, varClassification, group1, group2, nPCs, res.cluster, scale.factor,
#' normalization.method, selection.method, iso, time_stamp
#'
#' @import Seurat PCA process
#' @export

iSoMAs <- function(data.iso,data.maf,gene_to_iso=NULL, iso_to_gene=NULL,
                  genes_test=NULL,
                  stat.method="wilcox",
                  mut.num.thres = 5,mut.percent.thres = 0.02,
                  varType = "SNP",varClassification = NULL,
                  res.cluster = 0.5,
                  scale.factor = 10000,
                  normalization.method = "LogNormalize",
                  selection.method = "mvp",
                  nPCs = 50,return.Seurat=F,minP.sorted=F,
                  project = "TCGA-LUAD",
                  filename=NULL,
                  ntop.gene=100){


  # iSoMAs.R: main function for iSoMAs pipeline
  # copyright (c) Hua Tan, warm.tan@gmail.com
  # 11/2/2022 @4C08

  packages.iSoMAs = c("Matrix", "Seurat", "pbapply", "dplyr", "ggplot2", "ggpmisc", "ggrepel", "ggpubr", "pheatmap")
  tryCatch(install.packages(setdiff(packages.iSoMAs, rownames(installed.packages()))), error=function(e) e,
           finally = cat(paste0('\niSoMAs requires the following packages: ',paste(packages.iSoMAs,collapse = ', '),'\n')))

  library(Matrix)
  library(Seurat)
  library(pbapply)

  ###-----global variables
  jPC1=4   #the PC-coordinates starts with column 4, the first 3 columns are sample number of {WildType,Other,Mutant}
  group1 = "WildType" #test differential PC score between WildType and Mutant groups
  group2 = "Mutant"

  data.iso = trim_features_and_samples_TCGA(data.iso,gene_to_iso,iso_to_gene)

  ###-----numbers in comments are based on LUAD cancer-----###
  min.features = max(1,floor(nrow(data.iso)*0.1)) #5986
  min.cells = max(10,floor(ncol(data.iso)*0.05)) #25
  cat(paste0('min.cells=',min.cells,"; min.features=",min.features,"\n"))

  iso = run_PCA_Seurat(data.iso,res.cluster = res.cluster, scale.factor = scale.factor,
                       normalization.method = normalization.method,
                       selection.method = selection.method,
                       min.cells=min.cells, min.features=min.features,
                       nPCs=nPCs,project = project)

  pca = iso@reductions$pca
  myCoord = data.frame(pca@cell.embeddings)
  PCI = ncol(myCoord)

  if(any(rownames(myCoord) != colnames(data.iso))){
    stop("sample names after PCA are not consistent with data.iso")
  }else{
    mySamples = colnames(data.iso)
  }

  cat("Conducting MAF_PCA...\n")
  ##get genes_test from maf data

  if(is.null(genes_test)){
    cat("genes_test is NULL, generating genes_test from data.maf")
    mut.samp.thres = max(mut.num.thres,ceiling(ncol(data.iso)*mut.percent.thres)) #11
    cat(paste0('mut.samp.thres = ',mut.samp.thres,'\n'))

    data.maf.tmp = get_maf_of_interest_TCGA(data.maf,varType = varType,
                                            varClassification = varClassification)#199580, 120

    genes_maf = sort(table(data.maf.tmp$Hugo_Symbol),decreasing = T) #18888
    genes_maf = genes_maf[genes_maf>=mut.samp.thres] #5725

    cat('top 10 mutated genes:\n')
    print(head(genes_maf,10))
    cat('bottom 10 mutated genes:\n')
    print(tail(genes_maf,10))

    genes_test = names(genes_maf)
    rm(data.maf.tmp)
  }
  if(!is.null(ntop.gene)){
    cat(paste0('truncating to top ',ntop.gene,' genes\n'))
    genes_test = genes_test[1:min(length(genes_test),ntop.gene)]
  }
  cat(paste0("genes_test: ",length(genes_test),"---",paste(head(genes_test,20),collapse = ", "),"...\n"))

  if(length(genes_test)>0){
    pvals_all = pblapply(genes_test,function(myGene.mut){
      myCoord$Mutation = get_maf_of_interest_TCGA(data.maf,myGene=myGene.mut,
                                                  mySamples=mySamples,
                                                  varType=varType,
                                                  varClassification=varClassification,
                                                  verbose = F)$Mutation

      tabMut = table(myCoord$Mutation)
      # Mutant    Other WildType
      # 213       37      267
      if(any(!c("WildType","Mutant") %in% names(tabMut)) || tabMut[group1]<2 || tabMut[group2]<2){
        #by HT on 6/24/2022 @13547
        cat(paste0(Sys.time(),"---no sufficient samples for differential analysis for ",myGene.mut,"***\n"))
        print(tabMut)
        return(NULL)
      }

      ###-----in case there is no 'Other' samples
      tabMut0 = c(Mutant=0,Other=0,WildType=0)
      tabMut0[names(tabMut)] = tabMut
      tabMut = tabMut0
      # print(tabMut)

      if(tolower(stat.method)=="wilcox"){
        pvals = unlist(lapply(seq(PCI),function(j){
          wilcox.test(myCoord[myCoord$Mutation==group1,j],
                      myCoord[myCoord$Mutation==group2,j],
                      alternative="two.sided")$p.value
        }))
      }else if(tolower(stat.method)=="t.test"){
        pvals = unlist(lapply(seq(PCI),function(j){
          t.test(myCoord[myCoord$Mutation==group1,j],
                 myCoord[myCoord$Mutation==group2,j],
                 alternative="two.sided")$p.value
        }))
      }else{
        cat("stat.method must be {wilcox, t.test}\n")
      }
      return(c(tabMut,signif(pvals,2)))
    })
    names(pvals_all) = genes_test
    pvals_all = pvals_all[!unlist(lapply(pvals_all,is.null))]
    pvals_all = do.call("rbind",pvals_all)
    pvals_all = as.data.frame(pvals_all)
    colnames(pvals_all)[jPC1:ncol(pvals_all)] = paste0("PC_",seq(ncol(pvals_all)-jPC1+1))

    # cat("head of pvals_all:\n")
    # print(head(pvals_all[,]))
    # cat("tail of pvals_all:\n")
    # print(tail(pvals_all))
  }else{
    cat("since genes_test is empty, skipping MAF_PCA + Wilcoxon\n")
    pvals_all = NULL
  }

  if(return.Seurat){
    iso.rt = iso
  }else{
    iso.rt = NULL
  }

  if(minP.sorted){
    a = pvals_all[,jPC1:ncol(pvals_all)]
    a.min = apply(a,1,min)
    a.sorted = sort(a.min)
    pvals_all_sorted = pvals_all[names(a.sorted),,drop=F]
    pvals_all_sorted$minP = a.sorted
  }else{
    pvals_all_sorted = NULL
  }

  res_isomas = list(pca=pca, pvals_all=pvals_all, pvals_all_sorted=pvals_all_sorted,
                   genes_test=genes_test,
                   min.cells=min.cells,min.features=min.features,
                   dim_data.iso=dim(data.iso), dim_iso=dim(iso),
                   project=project, mut.samp.thres=mut.samp.thres,
                   varType=varType,varClassification=varClassification,
                   group1=group1,group2=group2,
                   nPCs=nPCs,res.cluster=res.cluster,scale.factor=scale.factor,
                   normalization.method = normalization.method,
                   selection.method=selection.method,
                   iso=iso.rt,
                   time_stamp=Sys.time())

  if(is.null(filename)){
    filename = paste0("res_SoMAS_",project,".RData")
  }

  save(res_isomas, file = filename)
  cat("OK\n")
  return(res_isomas)
}

