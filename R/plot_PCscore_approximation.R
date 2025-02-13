plot_PCscore_approximation <- function(res_PCapprox,plt.index=c(1,10,100,1000,2000,3000,3315),fs=16){
  # plot PC score approximation
  # res_PCapprox is generated by do_PCscore_approximation
  PCscore = res_PCapprox$PCscore
  ylab = names(PCscore)
  plt.clnms = paste0("X",plt.index)
  plist = lapply(plt.clnms,function(j){
    plot_corr_general(res_PCapprox$PCscore.approx[,j],
                      PCscore,
                      xlab = j,ylab = ylab,legend.pos = c(0.1,0.9),
                      fs = fs,xylim = F,Log2 = F)
  })
  names(plist) = plt.clnms
  # p = Seurat::CombinePlots(plist[1:min(6,length(plist))],ncol = 3)
  # print(p)
  return(plist)
}














