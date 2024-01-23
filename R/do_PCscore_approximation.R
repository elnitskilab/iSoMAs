do_PCscore_approximation <- function(res_isomas,iso,add.method = 'accumulatively',
                                     sign.PC='both',verbose = 100){
  # perform PC score approximation with the top features added accumulatively (
  # or individually)

  pca = res_isomas$pca
  pvals_all=res_isomas$pvals_all
  pvals_sig = get_pvals_sig(pvals_all,myGene = myGene.mut)
  PC = names(pvals_sig)[1]

  scale.data = iso@assays$RNA@scale.data #[1] 53963   517
  scale.data.loadings = t(scale.data[rownames(pca@feature.loadings),]) #[1] 517 3315
  feature.loadings = pca@feature.loadings #3315,50

  cat(paste0("scale.data: ",paste(dim(scale.data),collapse = ", "),"\n"))
  cat(paste0("scale.data.loadings: ",paste(dim(scale.data.loadings),collapse = ", "),"\n"))
  cat(paste0("feature.loadings: ",paste(dim(feature.loadings),collapse = ", "),"\n"))

  feature.loadings.PC = pca@feature.loadings[,PC,drop=T]

  # by Hua Tan on 7/27/2022 @4C08
  if(tolower(sign.PC) == "both"){
    ord = order(abs(feature.loadings.PC),decreasing = T)
  }else if(tolower(sign.PC) == "pos"){
    ord = order(feature.loadings.PC,decreasing = T)
  }else if(tolower(sign.PC) == "neg"){
    ord = order(feature.loadings.PC,decreasing = F)
  }else{
    cat('sign.PC must be {both, pos, neg}\n')
  }

  scale.data.loadings.ord = scale.data.loadings[,ord] #517 3315
  feature.loadings.ord = feature.loadings[ord,] #3315   50
  feature.loadings.PC.ord = feature.loadings.PC[ord] #3315

  Nsamp = nrow(pca@cell.embeddings) #517
  Nfeat = nrow(pca@feature.loadings) #3315
  if(tolower(add.method) == 'individually'){
    idx_kp_list = seq(Nfeat) #individually
    # xLab = 'Top n isoforms (individually)'
    xLab = 'Index (indiv.)'
  }else if (tolower(add.method) == 'accumulatively'){
    idx_kp_list = lapply(seq(Nfeat),function(i){seq(1,i)}) #accumulatively
    # xLab = 'Top n isoforms (accumulatively)'
    xLab = 'Index (accum.)'
  }else{
    cat('add.method must be {individually,accumulatively}\n')
  }

  cell.embeddings.pool = matrix(NA,nrow = Nsamp,ncol = length(idx_kp_list)) #517 3315
  for(i in seq(length(idx_kp_list))){
    idx_kp = idx_kp_list[[i]]
    if(verbose & i%%verbose==0){
      cat(paste0(paste(idx_kp,collapse = ', '),'\n\n'))
    }
    cell.embeddings.pool[,i] = scale.data.loadings.ord[,idx_kp,drop=F] %*%
      feature.loadings.ord[idx_kp,PC,drop=F] #[1] 517  1
  }
  cell.embeddings.pool = data.frame(cell.embeddings.pool)
  rownames(cell.embeddings.pool) = rownames(scale.data.loadings.ord)
  PCscore=as.numeric(res_isomas$pca@cell.embeddings[,PC])
  names(PCscore) = PC
  return(list(PCscore=PCscore,
              PCscore.approx=cell.embeddings.pool))
}
