
loci <- c("A", "B", "C", "DRB1", "DQB1", "DQA1", "DPB1", "DPA1", "DRB345")

dimal <- list(A=setdiff(seq(-24, 341), 0),
              B=setdiff(seq(-24, 338), 0),
              C=setdiff(seq(-24, 342), 0),
              DRB1=setdiff(seq(-29, 237), 0),
              DQB1=setdiff(seq(-32, 229), 0),
              DQA1=setdiff(seq(-23, 232), 0),
              DPB1=setdiff(seq(-29, 229), 0),
              DPA1=setdiff(seq(-31, 229), 0),
              DRB345=setdiff(seq(-29, 237), 0))

#' Get IMGT/HLA version
#'
#' This function takes a date and returns the IMGT/HLA version current at that
#' date. The function should work until 2024-01-01.
#'
#' @param x A string date in YYYY-MM-DD format
#' @return A string with the IMGT/HLA version without field separators
#' @export
IMGTHLA.version <- function(x){
  cuts <- as.Date(c("2010-04-01","2010-07-01","2010-10-01",
            "2011-01-01", "2011-04-01","2011-07-01","2011-10-01",
            "2012-01-01", "2012-04-01","2012-07-01","2012-10-01",
            "2013-01-01", "2013-04-01","2013-07-01","2013-10-01",
            "2014-01-01", "2014-04-01","2014-07-01","2014-10-01",
            "2015-01-01", "2015-04-01","2015-07-01","2015-10-01",
            "2016-01-01", "2016-04-01","2016-07-01","2016-10-01",
            "2017-01-01", "2017-04-01","2017-07-01","2017-10-01",
            "2018-01-01", "2018-04-01","2018-07-01","2018-10-01",
            "2019-01-01", "2019-04-01","2019-07-01","2019-10-01",
            "2020-01-01", "2020-04-01","2020-07-01","2020-10-01",
            "2021-01-01", "2021-04-01","2021-07-01","2021-10-01",
            "2022-01-01", "2022-04-01","2022-07-01","2022-10-01",
            "2023-01-01", "2023-04-01","2023-07-01","2023-10-01"))
  labs <- c("300","310","320","330",
            "340","350","360","370",
            "380","390","3100","3110",
            "3120","3130","3140","3150",
            "3160","3170","3180","3190",
            "3200","3210","3220","3230",
            "3240","3250","3260","3270",
            "3280","3290","3300","3310",
            "3320","3330","3340","3350",
            "3360","3370","3380","3390",
            "3400","3410","3420","3430",
            "3440","3450","3460","3470",
            "3480","3490","3500","3510",
            "3520","3530")
  as.character(cut.Date(as.Date(x), cuts, labs))
}

#' Get alignment files
#'
#' This function recovers the alignment files from the IMGT/HLA github account
#' (https://raw.githubusercontent.com/ANHIG/IMGTHLA/). The latest version is
#' downloaded by default.
#'
#' @param loci A string or vector of strings of the loci to be downloaded
#' @param dest A string with the destination folder for the downloaded files
#' @param IMGTHLA A string of the IMGT/HLA version to download the alignments for
#' @return NULL
#' @export
get.aln <- function(loci, dest="", IMGTHLA="Latest"){
  dir.create(dest)
  if(IMGTHLA=="Latest"){
    lnk <- "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/msf/"
  } else {
    lnk <- paste0("https://raw.githubusercontent.com/ANHIG/IMGTHLA/tree/", IMGTHLA, "/msf/")
  }
  for (l in loci){
    utils::download.file(paste0(lnk,l,"_prot.msf"), paste0(dest,l,"_prot.msf"))
  }
  NULL
}

#' Load the allele alignments
#'
#' This function loads downloaded allele alignments in .msf format.
#'
#' @param pth A string of the path to the folder containing the alignment files
#' @return A list of the alignments and complementary informations
#' @export
load.aln <- function(loci, dest=""){
  df <- list()
  msk <- list()
  num <- list()

  for (l in loci){
    df[[l]] <- bios2mds::import.msf(paste0(dest,l,"_prot.msf"), T, T)

    msk[[l]] <- df[[l]][[1]]!="-"

    num[[l]] <- rep(0, length(msk[[l]]))
    num[[l]][msk[[l]]] <- dimal[[l]]
  }
  return(list(aln=df, msk=msk, num=num))
}


id.missing <- function(a, msk){
  beg1 <- min(which(a!="-"))-1
  beg1 <- c(rep(T, beg1), rep(F, length(a)-beg1))
  beg2 <- min(which(msk))-1
  beg2 <- c(rep(T, beg2), rep(F, length(msk)-beg2))
  if(any(beg1 & !beg2)){
    a[beg1] <- "*"
  }

  end1 <- max(which(a!="-"))
  end1 <- c(rep(F, end1), rep(T, length(a)-end1))
  end2 <- max(which(msk))
  end2 <- c(rep(F, end2), rep(T, length(msk)-end2))
  if(any(end1 & !end2)){
    a[end1] <- "*"
  }

  a
}


closest.allele <- function(aln, wh, msk, a=NULL, imp=c("full", "step")){
  imp <- match.arg(imp)
  imp <- switch(imp, "full"=closest.allele.full, "step"=closest.allele.step)
  imp(aln, wh, msk, a)
}


closest.allele.step <- function(aln, wh, msk, a=NULL){
  if(is.null(a)){
    a <- id.missing(aln[[wh]], msk)
    mit <- sapply(aln, function(x){sum(id.missing(x, msk)=="*")})
    mth <- sum(id.missing(aln[[wh]], msk)=="*")
  } else {
    a <- id.missing(a, msk)
    mit <- sapply(aln, function(x){sum(id.missing(x, msk)=="*")})
    mth <- sum(id.missing(a, msk)=="*")
  }

  tmp <- sapply(aln, dif, seq2=a, gap=F)

  m <- min(tmp[names(tmp[mit<mth])], na.rm=T)
  tmp <- tmp[mit<mth][tmp[mit<mth]==m & !is.na(tmp[mit<mth])]
  tmp
}


closest.allele.full <- function(aln, wh, msk, a=NULL){
  if(is.null(a)){
    a <- id.missing(aln[[wh]], msk)
    mit <- sapply(aln, function(x){sum(id.missing(x, msk)=="*")})
  } else {
    a <- id.missing(a, msk)
    mit <- sapply(aln, function(x){sum(id.missing(x, msk)=="*")})
  }

  tmp <- sapply(aln, dif, seq2=a, gap=F)

  m <- min(tmp[names(tmp[mit==0])], na.rm=T)
  tmp <- tmp[mit==0][tmp[mit==0]==m & !is.na(tmp[mit==0])]
  tmp
}


impute.allele <- function(wh, aln, msk, imp=c("full", "step")){
  if("*" %in% id.missing(aln[[wh]], msk)){
    imp <- match.arg(imp)
    imp <- switch(imp, "full"=impute.allele.full, "step"=impute.allele.step)
    imp(wh, aln, msk)
  } else {
    aln[[wh]]
  }

}


impute.allele.full <- function(wh, aln, msk){
  a <- aln[[wh]]
  b <- aln[[names(closest.allele.full(aln, wh, msk))[1]]]

  mi <- id.missing(a, msk)=="*"

  if(any(mi)){
    a[mi] <- b[mi]
  }

  a
}


impute.allele.step <- function(wh, aln, msk){
  a <- aln[[wh]]
  b <- aln[[names(closest.allele.step(aln, wh, msk))[1]]]

  mi <- id.missing(a, msk)=="*"
  while(any(mi)){
    if(any(mi)){
      a[mi] <- b[mi]
    }

    mi <- id.missing(a, msk)=="*"
    if(any(mi)){
      b <- aln[[names(closest.allele.step(aln, wh, msk, a))[1]]]
    }
  }

  a
}


impute.whole <- function(db, impute=c("step", "full"), verbose=T){
  tmt <- tm <- Sys.time()
  impute=match.arg(impute)
  aln <- db$aln; num <- db$num; msk <- db$msk

  dbi <- list()
  for(locus in names(aln)){
    alni <- lapply(names(aln[[locus]]), function(al){
      a <- impute.allele(al, aln[[locus]], msk[[locus]], imp=impute)
    })
    names(alni) <- names(aln[[locus]])
    dbi[[locus]] <- alni
    if(verbose){print(paste("Locus", locus, "done (", round(difftime(Sys.time(), tmt, units = "mins"),2), "/", round(difftime(Sys.time(), tm, units = "mins"),2), "minutes )"))}
    tmt <- Sys.time()
  }
  list("aln"=dbi, "num"=num, "msk"=msk)
}


allele.diff <- function(a1, a2, msk, num, aln=NULL, locus=NULL){
  if(!is.null(locus)){
    msk <- msk[[locus]]
    num <- num[[locus]]
    aln <- aln[[locus]]
  }

  if(!is.null(aln)){
    a1 <- aln[[a1]]
    a2 <- aln[[a2]]
  }

  ma1 <- a1!="-"
  ma2 <- a2!="-"

  su1 <- ma1 &! msk
  su2 <- ma2 &! msk

  mi1 <- msk &! ma1
  mi2 <- msk &! ma2
  if(any(mi1|mi2)){
    mi <- xor(mi1, mi2)
    ni <- mapply(function(x,y,z){c(0,which(c(x,y)))[z+1]}, mi1, mi2, mi)

    for(i in 1:length(ni)){
      if(ni[i]==0){next()}
      j <- ni[i]
      r <- switch(j, 2, 1)
      tmp <- get(paste0("a", j))
      ref <- get(paste0("a", r))[i]
      tmp[i] <- ref
      assign(paste0("a", j), tmp)
    }
  }

  i <- mapply(identical, a1, a2)
  MM <- NULL
  for(j in which(!i)){
    MM <- c(MM, paste0(num[j], a1[j], ">", a2[j]))
  }
  MM
}


multi.compare <- function(aln, wh, msk, num, locus=NULL){
  if(!is.null(locus)){
    msk <- msk[[locus]]
    num <- num[[locus]]
    aln <- aln[[locus]]
  }

  A <- closest.allele(aln, wh, msk)
  res <- list()
  for (a in names(A)){
    res[[paste(a,wh, sep="//")]] <- allele.diff(a, wh, msk, num, aln)
  }
  res
}


view.align <- function(allele1, allele2, num, aln=NULL){
  if(!is.null(aln)){
    locus <- unlist(strsplit(allele1, "\\*"))[1]
    if(! allele1 %in% names(aln[[locus]])){
      allele1 <- allele.list.4f2l[[locus]][[allele1]]
    }
    if(! allele2 %in% names(aln[[locus]])){
      allele2 <- allele.list.4f2l[[locus]][[allele2]]
    }
    tmp <- rbind(aln[[locus]][[allele1]], aln[[locus]][[allele2]])
    rownames(tmp) <- c(allele1, allele2)
    colnames(tmp) <- num[[locus]]
  } else {
    tmp <- rbind(allele1, allele2)
    rownames(tmp) <- c("allele1", "allele2")
    colnames(tmp) <- num
  }

  tmp
}


view.closest <- function(allele, aln, msk, num, imp=c("full", "step")){
  imp <- match.arg(imp)
  imp <- switch(imp, "full"=closest.allele.full, "step"=closest.allele.step)
  locus <- unlist(strsplit(allele, "\\*"))[1]
  if(! allele %in% names(aln[[locus]])){
    allele <- allele.list.4f2l[[locus]][[allele]]
  }
  allele2 <- names(imp(aln[[locus]], allele, msk[[locus]]))[1]
  tmp <- rbind(aln[[locus]][[allele2]], aln[[locus]][[allele]])
  rownames(tmp) <- c(allele2, allele)
  colnames(tmp) <- num[[locus]]
  tmp
}
