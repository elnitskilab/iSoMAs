get_maf_of_interest_TCGA <- function(data.maf,myGene=NULL,mySamples=NULL,varType="SNP",
                                     varClassification=NULL,aa.mut=NULL,verbose=T){

  library(dplyr)

  data.maf$sample = substr(data.maf$Tumor_Sample_Barcode,1,16)

  if(!is.null(myGene)){
    data.maf = filter(data.maf,Hugo_Symbol %in% myGene) #LUAD, TP53: 294*120
    if(verbose){
      cat(paste0(myGene," has ", nrow(data.maf), " variants in total\n"))
    }

    if(!is.null(varType)){
      sample.other = unique(filter(data.maf,Variant_Type != varType)$sample)
      # cat(paste0("samples of Other = ",length(sample.other),"\n"))
    }else{
      sample.other = NULL
    }
  }

  if(!is.null(varType)){
    data.maf = filter(data.maf, Variant_Type %in% varType) #255, 121
    if(verbose){
      cat(paste0(myGene," has ", nrow(data.maf), " ", paste(varType,collapse = "/")," variants with varType designated\n"))
    }
  }

  if(!is.null(varClassification)){
    data.maf = filter(data.maf, Variant_Classification %in% varClassification)
    cat(paste0(myGene," has ", nrow(data.maf), " ", paste(varClassification,collapse = "/")," variants with varClassification designated\n"))
  }

  if(!is.null(aa.mut)){
    sample.other = union(sample.other,unique(filter(data.maf, !HGVSp_Short %in% aa.mut)$sample))
  }

  if(is.null(mySamples)){
    return(data.maf)
  }else{
    mut.samp = unique(data.maf$Tumor_Sample_Barcode) #243
    mut.samp = substr(mut.samp,1,16)
    Mutation = rep("Mutant",length(mut.samp))
    names(Mutation) = mut.samp
    samples = data.frame(barcode=mySamples,
                         sample=substr(mySamples,1,16))
    samples$Mutation = Mutation[samples$sample]
    samples$Mutation[is.na(samples$Mutation)] = "WildType"

    if(length(sample.other)>0){
      samples$Mutation[samples$sample %in% sample.other] = "Other"
    }
    rownames(samples) = samples$barcode
    if(verbose){
      print(table(samples$Mutation))
    }
    return(samples)
  }
}
