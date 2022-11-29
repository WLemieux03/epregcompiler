
load.epreg <- function(path){
  sn <- readxl::excel_sheets(path)
  ep <- list()
  for (n in sn){
    ep[[n]] <- readxl::read_excel(path, n)
  }
  ep
}


clean.or.pre <- function(ep){
  sn <- names(ep)
  sw <- function(x) switch(x, "DR11"="DRB", "DP"="DPB")
  adj <- NULL
  for (n in sn){
    ide <- grepl("or", ep[[n]]$Polymorphic)
    if(any(ide)){
      tmp <- ep[[n]][ide,]

      nam <- grep("or", ep[[n]]$Polymorphic, value=T)
      ep[[n]]$Polymorphic[ide] <- sapply(nam, function(x){unlist(strsplit(x, " or "))[1]})

      nam <- sapply(nam, function(x){unlist(strsplit(x, " or "))[2]})
      lo <- sapply(nam, function(x){sw(gsub("[\\s\u00A0]","",unlist(strsplit(x, " on "))[2]))})
      tmp$Polymorphic <- sapply(nam, function(x){unlist(strsplit(x, " on "))[1]})

      for (i in 1:dim(tmp)[1]){
        ep[[lo[i]]] <- rbind(ep[[lo[i]]], tmp[i,])
      }
      adj <- rbind(adj, data.frame("locus"=rep(n, dim(tmp)[1]), "toCorr"=lo, "Eplet"=tmp$Eplet, "abv"=is.abv(tmp$Antibody)))
    }
  }
  list("ep"=ep, "adj"=adj)
}


clean.or.post <- function(epmat, adj){
  if(is.null(adj)){
    return(epmat)
  }
  for (i in 1:dim(adj)[1]){
    rownames(epmat[[adj$toCorr[i]]]) <- sub(paste(adj$toCorr[i], adj$Eplet[i], adj$abv[i], sep="\\."), paste(adj$locus[i], adj$Eplet[i], adj$abv[i], sep="\\."), rownames(epmat[[adj$toCorr[i]]]))
  }
  epmat
}


is.abv <- function(X){
  c("oth", "abv")[((X=="Yes" | X=="Confirmed" | X=="Provisional") & !is.na(X))+1]
}


ep.matrix <- function(ep, strict=T, supp=T){
  sn <- setdiff(names(ep), "DRDQDP")
  AA <- list()
  for (n in sn){
    aa <- stringr::str_extract_all(ep[[n]]$Polymorphic, "\\d+[[[:alpha:]]/]+")
    aa0 <- stringr::str_extract_all(stringr::str_extract_all(ep[[n]]$Polymorphic, "\\(.+\\)"), "\\d+[[[:alpha:]]/]+")
    aa1 <- mapply(setdiff, aa, aa0)
    if(strict){aa <- aa1}
    an <- lapply(stringr::str_extract_all(aa, "\\d+"), as.integer)
    ai <- lapply(aa, function(x){unlist(stringr::str_extract_all(x, "\\D+"))})
    epnam <- paste(n,ep[[n]]$Eplet, is.abv(ep[[n]]$Antibody), sep=".")
    if(supp & (n %in% c("DRB", "DQB", "DPB"))){
      ptrn <- switch(n, "DRB"="^[qp]*r", "DQB"="^[rp]*q", "DPB"="^[rq]*p")
      sep <- ep$DRDQDP[grep(ptrn, ep$DRDQDP$Eplet),]
      saa <- stringr::str_extract_all(sep$Polymorphic, "\\d+[[[:alpha:]]/]+")
      saa0 <- stringr::str_extract_all(stringr::str_extract_all(sep$Polymorphic, "\\(.+\\)"), "\\d+[[[:alpha:]]/]+")
      saa1 <- mapply(setdiff, saa, saa0)
      if(strict){saa <- saa1}
      san <- lapply(stringr::str_extract_all(saa, "\\d+"), as.integer)
      if(n=="DPB"){
        cuts <- c(-Inf, 25, 156, Inf)
        labs <- c(0, 2, 3)
        san <- lapply(san, function(x){x-as.integer(as.character(cut(x, breaks=cuts, labels=labs)))})
      }
      sai <- lapply(saa, function(x){unlist(stringr::str_extract_all(x, "\\D+"))})
      sepnam <- paste("DRDQDP",sep$Eplet, c("oth", "abv")[((sep$Antibody=="Yes" | sep$Antibody=="Confirmed" | sep$Antibody=="Provisional") & !is.na(sep$Antibody))+1], sep=".")
      an <- c(an, san)
      ai <- c(ai, sai)
      epnam <- c(epnam, sepnam)

    }
    aaref <- sort(unique(unlist(an)))

    ai <- lapply(1:length(ai), function(x){ai[[x]][order(an[[x]])]})
    an <- lapply(1:length(ai), function(x){an[[x]][order(an[[x]])]})

    aamat <- t(mapply(function(x,y){
      res <- rep(NA, length(aaref))
      res[aaref %in% x] <- y
      res
    }, an, ai))
    rownames(aamat) <- epnam
    colnames(aamat) <- aaref

    AA[[n]] <- aamat
  }
  AA
}


ep.calc <- function(a, num, aamat){
  sapply(rownames(aamat), function(x){
    ptrn <- sub("/", "|", aamat[x,!is.na(aamat[x,])])
    aam <- a[num %in% as.integer(colnames(aamat)[!is.na(aamat[x,])])]
    all(mapply(function(y,z){grepl(y,z)}, ptrn, aam))
  })*1
}


ind.ep.calc <- function(ind, num, aamat, aln, msk, restrict=NULL, as.list=F, impute=T){
  EPR <- NULL
  res <- lapply(colnames(ind), function(x){
    epr <- switch(locus, "A"="ABC", "B"="ABC", "C"="ABC", "DRB1"="DRB",
                  "DQB1"="DQB", "DQA1"="DQA", "DPB1"="DPB", "DPA1"="DPA",
                  "DRB3"="DRB", "DRB4"="DRB", "DRB5"="DRB", "DRB345"="DRB")
    EPR <<- c(EPR, epr)
    emm <- sapply(ind[[x]], function(y){
      if(is.na(y) | grepl("N$", y)){
        if(length(epr)==1){
          rep <- rownames(aamat[[epr]])
          tmp <- rep(0, length(rep))
          names(tmp) <- rep
          return(tmp)
        } else{
          rep <- unlist(lapply(epr, function(e){rownames(aamat[[e]])}))
          tmp <- rep(0, length(rep))
          names(tmp) <- rep
          return(tmp)
        }
      }
      if(impute){
        a <- impute.allele(y, aln[[x]], msk[[x]])
      } else {
        a <- aln[[x]][[y]]
      }

      if(length(epr)==1){
        return(ep.calc(a, num[[x]], aamat[[epr]]))
      } else{
        foreach(z=epr, .combine=c) %do% {ep.calc(a, num[[x]], aamat[[z]])}
      }

    })
    apply(emm, 1, max)
  })
  names(res) <- colnames(ind)

  if(as.list){return(res)}

  if(is.null(restrict)){
    EPR <- names(aamat)[names(aamat) %in% EPR]
  } else {
    EPR <- restrict
  }

  rep <- unlist(lapply(EPR, function(e){rownames(aamat[[e]])}))
  tmp <- sapply(rep, function(x){
    sapply(res, function(y){if(x %in% names(y)){y[[x]]}else{0}})
  })
  apply(tmp, 2, max)
}


pop.ep.calc <- function(pop, num, aamat, aln, msk, loci=c("A","B","C","DRB1","DQB1"), restrict=NULL, as.list=F, impute=T){
  res <- apply(pop, 1, function(p){
    ind <- list()
    for (locus in loci){
      ind[[locus]] <- p[paste0(locus, c("_1", "_2"))]
      ind[[locus]] <- sapply(ind[[locus]], function(a){grep(sub("\\*","\\\\*",a), names(df[[locus]]), value=T)[1]})
    }
    ind <- data.frame(ind)
    ind.ep.calc(ind, num, aamat, aln, msk, restrict=NULL, as.list=F, impute=T)
  })
  t(res)
}


ep.dict <- function(aln, msk, num, aamat, restrict=NULL, as.list=F, impute=NULL, verbose=T){
  tmt <- tm <- Sys.time()
  dict <- list()
  for(locus in names(aln)){
    epr <- switch(locus, "A"="ABC", "B"="ABC", "C"="ABC", "DRB1"="DRB",
                  "DQB1"="DQB", "DQA1"="DQA", "DPB1"="DPB", "DPA1"="DPA",
                  "DRB3"="DRB", "DRB4"="DRB", "DRB5"="DRB", "DRB345"="DRB")
    emm <- sapply(names(aln[[locus]]), function(al){
      if(grepl("N$", al)){
        if(length(epr)==1){
          rep <- rownames(aamat[[epr]])
          tmp <- rep(0, length(rep))
          names(tmp) <- rep
          return(tmp)
        } else{
          rep <- unlist(lapply(epr, function(e){rownames(aamat[[e]])}))
          tmp <- rep(0, length(rep))
          names(tmp) <- rep
          return(tmp)
        }
      }
      if(!is.null(impute)){
        a <- impute.allele(al, aln[[locus]], msk[[locus]], imp=impute)
      } else {
        a <- aln[[locus]][[al]]
      }

      if(length(epr)==1){
        return(ep.calc(a, num[[locus]], aamat[[epr]]))
      } else{
        foreach(z=epr, .combine=c) %do% {ep.calc(a, num[[locus]], aamat[[z]])}
      }
    })
    dict[[locus]] <- t(emm)
    if(verbose){print(paste("Locus", locus, "done (", round(difftime(Sys.time(), tmt, units = "mins"),2), "/", round(difftime(Sys.time(), tm, units = "mins"),2), "minutes )"))}
    tmt <- Sys.time()
  }
  dict
}


which.ep <- function(aamat){
  lapply(aamat, rownames)
}


ind.ep.match <- function(ind, dict, eplist, restrict=NULL){
  EPR <- NULL
  res <- lapply(colnames(ind), function(x){
    epr <- switch(locus, "A"="ABC", "B"="ABC", "C"="ABC", "DRB1"="DRB",
                  "DQB1"="DQB", "DQA1"="DQA", "DPB1"="DPB", "DPA1"="DPA",
                  "DRB3"="DRB", "DRB4"="DRB", "DRB5"="DRB")
    EPR <<- c(EPR, epr)
    emm <- sapply(ind[[x]], function(y){
      if(is.na(y)){
        tmp <- rep(0, dim(dict[[x]])[2])
        names(tmp) <- colnames(dict[[x]])
        return(tmp)
      }

      return(dict[[x]][y,])
    })
    apply(emm, 1, max)
  })
  names(res) <- colnames(ind)


  if(is.null(restrict)){
    EPR <- names(eplist)[names(eplist) %in% EPR]
  } else {
    EPR <- restrict
  }

  rep <- unlist(lapply(EPR, function(e){eplist[[e]]}))
  tmp <- sapply(rep, function(x){
    sapply(res, function(y){if(x %in% names(y)){y[[x]]}else{0}})
  })
  apply(tmp, 2, max)
}


pop.ep.match <- function(pop, dict, eplist, loci=c("A","B","C","DRB1","DQB1"), restrict=NULL){
  res <- apply(pop, 1, function(p){
    ind <- list()
    for (locus in loci){
      ind[[locus]] <- p[paste0(locus, c("_1", "_2"))]
      ind[[locus]] <- sapply(ind[[locus]], function(a){grep(sub("\\*","\\\\*",a), rownames(dict[[locus]]), value=T)[1]})
    }
    ind <- data.frame(ind)
    ind.ep.match(ind, dict, eplist, restrict=NULL)
  })
  t(res)
}

#ind <- data.frame(A=c("A*02:01:01:01", "A*03:02:01:01"), B=c("B*18:07:01","B*07:223"), DRB1=c("DRB1*04:105:02","DRB1*07:90"), DQB1=c("DQB1*05:01:01:01", "DQB1*06:01:01:01"))
#ind <- data.frame(A=c("A*01:01:01:01", "A*01:01:01:02N"), B=c("B*07:02:01:01","B*07:02:01:01"), DRB1=c("DRB1*01:01:01:01","DRB1*01:01:01:01"))


tst.allele <- function(allele, epdict, conversion, dest=""){
  dir.create(dest, F)
  locus <- unlist(strsplit(allele, "\\*"))[1]
  if(! allele %in% rownames(epdict[[locus]])){
    allele <- conversion[[locus]][[allele]]
  }
  write(colnames(epdict[[locus]])[epdict[[locus]][allele,]==1], paste0(dest,gsub(":","-",gsub("\\*", "_", allele)), ".txt"))
}


to.xlsx <- function(epdict, out="epdict.xlsx"){
  d <- lapply(epdict, as.data.frame)
  d <- lapply(d, function(x){cbind(allele=rownames(x), x)})
  writexl::write_xlsx(d, out)
  return(NULL)
}


to.AAlist <- function(epmat, out="aalist.xlsx"){
  epsp <- lapply(epmat, function(x){
    data.frame("AA"=apply(x, 1, function(y){
      paste(paste0(names(y), y)[!is.na(y)], collapse=" ")
    }))
  })
  to.xlsx(epsp, out)
}


ep.all.alleles <- function(ep, luminex=F){
  sn <- names(ep)
  a <- NULL
  cn <- c("All Alleles", "Luminex Alleles")[luminex+1]
  for (n in sn){
    a <- c(a, unique(unlist(strsplit(unlist(ep[[n]][,cn]), ", "))))
  }
  loc <- sort(unique(sapply(a, function(x){strsplit(x, "\\*")[[1]][1]})))
  L <- list()
  for(l in loc){
    L[[l]] <- unique(sort(grep(paste0("^",l, "\\*"), a, value=T)))
  }
  L
}


ep.alleles <- function(ep, luminex=F){
  sn <- names(ep)
  EP <- list()
  cn <- c("All Alleles", "Luminex Alleles")[luminex+1]
  for (n in sn){
    EP[[n]] <- strsplit(unlist(ep[[n]][,cn]), ", ")
    names(EP[[n]]) <- unlist(ep[[n]][,"Eplet"])
  }
  EP
}


dict.from.epreg <- function(ep, supp=T, luminex=F){
  alleles <- ep.all.alleles(ep, luminex)
  eplets <- ep.alleles(ep, luminex)
  AA <- list()
  for (locus in names(alleles)){
    epr <- switch(locus, "A"="ABC", "B"="ABC", "C"="ABC", "DRB1"="DRB",
                  "DQB1"="DQB", "DQA1"="DQA", "DPB1"="DPB", "DPA1"="DPA",
                  "DRB3"="DRB", "DRB4"="DRB", "DRB5"="DRB")
    aamat <- t(sapply(alleles[[locus]], function(x){sapply(eplets[[epr]], function(y){x %in% y})}))
    epnam <- paste(epr,names(eplets[[epr]]), c("oth", "abv")[((ep[[epr]]$Antibody=="Yes" | ep[[epr]]$Antibody=="Confirmed" | ep[[epr]]$Antibody=="Provisional") & !is.na(ep[[epr]]$Antibody))+1], sep=".")

    if(supp & (epr %in% c("DRB", "DQB", "DPB"))){
      ptrn <- switch(epr, "DRB"="^[qp]*r", "DQB"="^[rp]*q", "DPB"="^[rq]*p")
      sep <- grepl(ptrn, ep$DRDQDP$Eplet)

      saamat <- t(sapply(alleles[[locus]], function(x){sapply(eplets$DRDQDP[sep], function(y){x %in% y})}))
      sepnam <- paste("DRDQDP",names(eplets$DRDQDP[sep]), c("oth", "abv")[((ep$DRDQDP[sep,]$Antibody=="Yes" | ep$DRDQDP[sep,]$Antibody=="Confirmed" | ep$DRDQDP[sep,]$Antibody=="Provisional") & !is.na(ep$DRDQDP[sep,]$Antibody))+1], sep=".")
      aamat <- cbind(aamat, saamat)
      epnam <- c(epnam, sepnam)

    }
    colnames(aamat) <- epnam

    AA[[locus]] <- aamat*1
  }
  AA
}


ep.from.dict <- function(allele, dict){
  locus <- unlist(strsplit(allele, "\\*"))[1]
  if(!allele %in% rownames(dict[[locus]])){
    ch <- grep(sub("\\*", "\\\\\\*", allele), rownames(dict[[locus]]), value=T)
    if(length(ch)>16){stop("Too much choices!")}
    al <- readline(paste(c("Which specificity?", ch), collapse="\n"))
  } else {
    al <- allele
  }
  ep <- dict[[locus]][al,]
  names(ep)[ep==1]
}


compound.dicts <- function(dict1, dict2){
  loci <- union(names(dict1), names(dict2))
  dict <- list()
  for(locus in loci){
    dict[[locus]] <- (dict1[[locus]]+dict2[[locus]])/2
  }
  dict
}


which.has <- function(dict, locus, ep){
  if(!ep %in% colnames(dict[[locus]])){
    ep <- grep(ep, colnames(dict[[locus]]), value=T)
    print(ep)
  }
  invest <- dict[[locus]][,ep]
  names(invest)[invest==1]
}


twoField <- function(dict){
  sn <- names(dict)
  for (n in sn){
    al.2f <- unique(unlist(str_extract_all(rownames(dict[[n]]), "^\\D{1,3}\\d?\\*\\d{2,4}:\\d{2,4}")))
    al.4f <- sapply(al.2f, function(a){grep(sub("\\*","\\\\*",a), rownames(dict[[n]]), value=T)[1]})
    dict[[n]] <- dict[[n]][al.4f,]
    rownames(dict[[n]]) <- al.2f
  }
  dict
}
