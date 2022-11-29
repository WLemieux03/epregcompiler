#' Load the eplet dictionary
#'
#' This function loads a precalculated dictinary for allele to eplet conversion.
#'
#' @param pth A string of the path of the .xlsx file containing the dictionary.
#' @return A list of dataframes of the dictionary.
#' @export
load.dict <- function(pth){
  sn <- readxl::excel_sheets(pth)
  dict <- list()
  for (n in sn){
    dict[[n]] <- as.data.frame(readxl::read_excel(pth, n))
    rownames(dict[[n]]) <- dict[[n]]$allele
    dict[[n]] <- dict[[n]][,-1]
  }
  dict
}

#' Assign epitopes to a genotype
#'
#' This function calculates the epitype from a genotype.
#'
#' @param ind A named vector of strings with a genotype. The names of \code{ind} should be of the format "loci_d", Where "loci" is the loci of the allele and d is either "1" or "2".
#' @return A list of dataframes of the dictionary.
#' @export
assign.epitope <- function(ind, dict, loci=c("A", "B", "C", "DRB1", "DQB1", "DQA1", "DPB1", "DPA1"), binary=T){
  loc <- sapply(ind, function(i){unlist(strsplit(i, "\\*"))[1]})
  epn <- unique(unlist(sapply(loci, function(l){colnames(dict[[l]])})))
  ep <- rep(0, length(epn)); names(ep) <- epn

  ind <- sub(".+N$", NA, ind)
  ind <- sub("[QL]$", "", ind)
  ind <- stringr::str_extract(ind, "[ABCDRQP]{1,3}1?\\*\\d{2,3}:\\d{2,3}")

  for (i in 1:length(ind)){
    if(is.na(ind[i])){
      next
    }
    tmp <- dict[[loc[i]]][ind[i],]
    ep[names(tmp)] <- unlist(ep[names(tmp)]+tmp)
  }
  if(binary){
    ep <- (!(ep==0))*1
  }
  ep
}

#' Assign epitopes to a dataframe of genotypes
#'
#' This function calculates the epitype from a dataframe containing individuals as rows and genotypes in specifically named columns. Supplementary information and genotypes are kept.
#'
#' @param df A dataframe of the genotypes. The names of columns containing genotypes should be of the format "loci_d", Where "loci" is the loci of the allele and d is either "1" or "2".
#' @param dict A list of dataframes of the dictionary containing the conversion between alleles and epitype.
#' @param loci A string or vector of strings of the loci to be considered in the epitype.
#' @param binary A logical indicating whether the sum or a 0/1 value should be returned for each eplet in the epitype.
#' @param mc A logical indicating whether parallel processing should be used.
#' @return A dataframe containing the information of \code{df} with the epytipe appended as supplementary columns.
#' @export
multi.assign <- function(df, dict, loci=c("A", "B", "C", "DRB1", "DQB1", "DQA1", "DPB1", "DPA1"), binary=T, mc=F){
  coln <- apply(expand.grid(c(1,2), loci), 1, function(x){paste(x[2], x[1], sep="_")})
  if(mc){
    return(cbind(df, t(apply(df[,intersect(coln, colnames(df))], 1, assign.epitope, dict=dict, loci=loci, binary=binary))))
  } else {
    plan(cluster, workers=parallel::makeCluster(Sys.getenv("SLURM_CPUS_PER_TASK")))
    return(cbind(df, t(future.apply::future_apply(df[,intersect(coln, colnames(df))], 1, assign.epitope, dict=dict, loci=loci, binary=binary))))
  }
}


