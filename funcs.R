blya <- function(all_data, cellcount, vcellcount, pop_ls, srch_v, kids_names, pop, verbose){
#Checks if parent population is really above in hierarchy
  out <- NULL
  ndCol <- 0
  ndKidCol <- NULL
  node <- NULL
  out_nodes <- NULL
  count1 <- 0
  count2 <- 0
  forDf <- NULL
  forName <- NULL
  for (i in pop_ls){
    #print(i)
    count1 <- count1+1
    for (j in srch_v){
     # print(j)
      if (length(grep(j, i)) != 0 & length(grep(pop, i)) != 0){
        if (grep(j, i) - grep(pop, i) == 1 & i[grep(j, i)] == tail(i, 1)){
          out <- c(out, count1)
          if (verbose == T) {print("Nodes checking")}
        }
        if (grep(j, i) - grep(pop, i) > 1 & i[grep(j, i)] == tail(i, 1)){
              node <- i[grep(j, i)-1]
              if (verbose == T) {print(c("For",j ,"Node found", node))}
          a <- strsplit(nicename(vcellcount[grep(node, nicename(vcellcount))]), ",")
          na <- grep(node, nicename(vcellcount))
          count2 <- 0
          for (z in a){
            count2 <- count2 + 1
            if (length(grep(pop, z)) != 0){
              if (z[grep(node, z)] == tail(z, 1)){
                out_nodes <- all_data[names(all_data) == cellcount[na[count2]]]
                ndCol <- out_nodes
                #forName <- paste(nicename(cellcount[na[count2]]),"node")
                #assign(tempdf, forDf, envir = .GlobalEnv)
                ndKidCol <- all_data[gsub(" ", "", nicename(names(all_data))) == kids_names[count1]]
                forDf <- c(forDf, (ndKidCol/ndCol)*100)
                forName <- c(forName, paste(j, node, pop, sep="_"))
              }
            }
          }
        }
      }
    }
    #cat(paste("Internal nodes were found. Following populations will be added:\n", forName))
  }
  forDf <- as.data.frame(forDf)
  names(forDf) <- forName
  #assign(real_kids_cols, out, envir = .GlobalEnv)
  return(list(out, forDf))
}
nicename <- function(i){
  name <- gsub("[+]","pos", i)
  name <- gsub("[-]", "neg", name)
  name <- gsub(" [|]", "", name)
  name <- gsub(" Count", "", name)
  return(name)
}
filenameparser <- function(df, flname = "exp_ID date tissue genotype treatment sample_ID"){
    if (flname == casefold("no")){
        return(df)
        }else{
            map <- strsplit(flname, " ")[[1]]
            col <- strsplit(df$filename, "_")
            col <- col[col != c("Mean", "SD")]
            col <- sapply(col, function(x) x[-length(x)])
            v <- grep("ICs|Unstained|Isotype|unst|iso|IC", col, invert = TRUE, ignore.case = T)
            newdf <- as.data.frame(lapply(seq_along(map), function(i) sapply(col[v], "[[", i)))
            names(newdf) <- map
            newdf$ICs <- rep(0, nrow(newdf))
            newdf$Unstained <- rep(0, nrow(newdf))
            newdf1 <- as.data.frame(lapply(seq_along(col[-v][[1]]), function(i) sapply(col[-v], "[[", i)))
            #Process ICs
            ic <- sapply(newdf1, function(x) grep("ICs|Isotype|iso|IC", x, ignore.case = T))
            ic.logic <- sapply(ic[ic>0], function(x) is.numeric(x)[[1]])
            ic <- ic[ic.logic][[1]]
            ic.col <- sapply(seq_along(ic.logic), function(x){ if (ic.logic[x] == TRUE){return(x)}})
            ic.col <- ic.col[sapply(ic.col, is.integer)][[1]]
            newdf.IC <- newdf1[ic,1:ic.col-1]
            newdf.IC[ic.col:length(map)] <- NA
            names(newdf.IC) <- map
            newdf.IC$ICs <- rep(1, nrow(newdf.IC))
            newdf.IC$Unstained <- rep(0, nrow(newdf.IC))
            newdf.IC <- cbind(df[ic,], newdf.IC)
            sensedf <- cbind(df[v,], newdf)
            sensedf <- rbind(newdf.IC,sensedf)
            #Process Unstained
            unst <- sapply(newdf1, function(x) grep("Unstained|unst", x, ignore.case = T))
            unst.logic <- sapply(unst[unst>0], function(x) is.numeric(x)[[1]])
            unst <- unst[unst.logic][[1]]
            unst.col <- sapply(seq_along(unst.logic), function(x){ if (unst.logic[x] == TRUE){return(x)}})
            unst.col <- unst.col[sapply(unst.col, is.integer)][[1]]
            newdf.UN <- newdf1[ic,1:ic.col-1]
            newdf.UN[unst.col:length(map)] <- NA
            names(newdf.UN) <- map
            newdf.UN$ICs <- rep(0, nrow(newdf.IC))
            newdf.UN$Unstained <- rep(1, nrow(newdf.IC))
            newdf.UN <- cbind(df[unst,], newdf.UN)
            sensedf <- rbind(newdf.UN,sensedf)
            return(sensedf)
            }
}
popConstruct <- function(popul_file, all_data, verbose = FALSE){
  cellcount <- names(all_data)[grep("Count", x = names(all_data))]
  ncellcount <- gsub(c(" Count"), "", cellcount)
  ncellcount <- gsub(" [|]", "", ncellcount)
  ncellcount <- gsub(" ", "", ncellcount, perl = TRUE)
  vcellcount <- strsplit(ncellcount, "/")
  inpops <- file(popul_file, open = "r")
  read <- readLines(inpops)
  close(inpops)
  param <- sapply(1:length(read), filePrep, read = read)
  param <- param[!sapply(param, is.null)]
  #param <- sapply(1:length(param), function(i) param[i] <- sapply(param[i], function(x) x[[1]][!is.na(x[[1]])]))
  finaldf <- NULL
  resname <- NULL
  forDF <- NULL
  #while(TRUE){
  populations <- param[sapply(1:length(param),function(i) names(param[i][[1]]) == "populations")]
  populations <- populations[[1]][[1]][!is.na(populations[[1]][[1]])]
  everything <- sapply(populations, function(l){
    if (verbose==T) cat("========= Start of line ============\n")
    #print(l)
    rd <- strsplit(l, "_", fixed = TRUE)
    if (verbose==T) print(c("Population list", rd))
    pop <- rd[[1]][1]
    if (verbose==T) print(c("Main population name", pop))
    nrd <- gsub(" ", "", strsplit(rd[[1]][2],"-")[[1]])
    #print(nrd)
    par <- strsplit(gsub(" ", "", nrd[1]), ",")
    if (verbose==T) print(c("Parent populations name", par))
    #Current population
    kids <- NULL
    kids_cols <- NULL
    kids_ls <- NULL
    kids_names <- NULL
    parents_cols <- NULL
    parents_ls <- NULL
    parents_names <- NULL
    pcol <- NULL
    plab <- NULL
    tcol <- NULL
    tlab <- NULL
    v <- NULL
    vect <- NULL
    if (length(nrd) > 1){
      kids <- gsub(" ", "", strsplit(nrd[2], ",")[[1]])
      if (verbose==T) print(c(length(kids), "Inherited populations", kids))
    } else{
      kids <- NULL
      kids_cols <- NULL
      kids_ls <- NULL
      kids_names <- NULL
    }
    #Parent population to compute frequency
    #Inherited populations to compute a percentage of pop
    #Put all data needed to a temporary dataframe to organize it
    tempdf <- NULL
    name <- NULL
    for (i in  gsub(" ", "", nicename(cellcount[grep(nicename(pop), gsub(" ", "", nicename(cellcount)))]))){
      tempdf <- c(tempdf, all_data[gsub(" ", "", nicename(names(all_data))) == i])
      if (length(grep(nicename(gsub(" ","",pop)), tail(strsplit(i, "/")[[1]], 1))) > 0){
        name <- c(name, paste(nicename(i),"target"))
      }else{
      name <- c(name, nicename(i))
      }
    }
    vect <- NULL
    for (i in par[[1]]){
      v <- NULL
      for (j in vcellcount){
        v <- c(v, i %in% tail(j,1))
      }
      vect <- c(vect, v)
      tempdf <- c(tempdf, all_data[names(all_data) == cellcount[v]])
      name <- c(name, paste(nicename(cellcount[v]),"parent"))
    }
    tempdf <- as.data.frame(tempdf)
    names(tempdf) <- name
    parents_names <- names(tempdf)[grep("parent", names(tempdf))]
    parents_cols <- grep("parent", names(tempdf))
    kids <- nicename(kids)
    pop <- nicename(pop)
    par <- nicename(par)
    kids_names <- names(tempdf)[grep(paste(c("parent", "target"), collapse = "|"), names(tempdf), invert = TRUE)]
    kids_cols <- grep(paste(c("parent", "target"), collapse = "|"), names(tempdf), invert = TRUE)
    kids_ls <- strsplit(kids_names, "/")
    parents_ls <- strsplit(parents_names, "/")
    pcol <- grep("parent", names(tempdf))
    tcol <- grep("target", names(tempdf))
    resdf <- NULL
    for (i in pcol){
      resdf <- c(resdf,tempdf[,tcol]/tempdf[,i]*100)
      tlab <- strsplit(tail(strsplit(names(tempdf)[tcol],"/")[[1]],1)," ")[[1]][1]
      plab <- strsplit(tail(strsplit(names(tempdf)[pcol],"/")[[1]],1)," ")[[1]][1]
      resname <- c(resname, paste(tlab,plab, sep="_"))
      if (verbose==T) {print(paste("resname, parents cycle: ", resname))}
    }
    bLYA <- blya(all_data,cellcount,vcellcount,kids_ls,kids,kids_names,pop, verbose)
    real_kids_cols <- bLYA[[1]]
    if (verbose==T) {cat(paste("**********\n Real kids cols value\n", real_kids_cols, "**********"))}
    for (i in real_kids_cols){
        if (verbose==T) {print("kids cycle entered")}
      resdf <- list(resdf,(tempdf[names(tempdf) == kids_names[i]])[[1]]/tempdf[,tcol]*100)
      tlab <- tail(strsplit(kids_names[i],"/")[[1]],1)
      plab <- strsplit(tail(strsplit(names(tempdf)[tcol],"/")[[1]],1)," ")[[1]][1]
      resname <- c(resname, paste(tlab,plab, sep="_"))
      if (verbose==T) {print(paste("resname, kids cycle: ", resname))}
    }
    resdf <- as.data.frame(resdf)
    #print(resdf)
    if (verbose==T) {cat("\nCreated temporary Data frame\n")}
    if ( length(line) == 0 ) {
      break
    }
    if (verbose==T) cat("========= End of line ============\n\n")

    finaldf <- resdf
    forDF <- bLYA[[2]]
  finaldf <- as.data.frame(finaldf)
  names(finaldf) <- resname
  forDF <- as.data.frame(forDF)
  finaldf <- c(finaldf, forDF)
  finaldf <- as.data.frame(finaldf)
  return(finaldf)
  })
  everything.df <- NULL
  for (x in 1:length(everything)){
      everything.df <- c(everything.df, everything[[x]])}
  everything.df <- as.data.frame(everything.df)
  everything.df$filename <- all_data$filename
  flname <- param[sapply(1:length(param), function(i) names(param[i][[1]]) == "flname")][[1]][[1]]
  everything.df <- everything.df[grep("Mean|SD", everything.df$filename, invert = T),]
  everything.df <- filenameparser(everything.df, flname)
  #Calculations with ICs
  ics <- param[sapply(1:length(param), function(i) names(param[i][[1]]) == "isotypes")][[1]][[1]]
  ics <- gsub(" ", "", nicename(strsplit(ics, ",")[[1]]))
  icsDF <- everything.df[everything.df$ICs == 1,]
  unstDF <- everything.df[everything.df$Unstained == 1,]
  #Find tissues or exps of dates to substract ICs correctly
  map <- strsplit(casefold(flname), " ")[[1]]
  t <- grep("tissue|organ|org", map)
  e <- grep("exp|manip", map)
  d <- grep("date|day", map)
  n <- names(everything.df)
  l <- 1:length(n)
  lst <- c(l[casefold(n) == map[t]],l[casefold(n) == map[e]],l[casefold(n) == map[d]])
  #Substract ICs
  icsDF.calc <- sapply(1:length(everything.df[grep(paste(ics, collapse = "|"), names(everything.df))]), function(x) {
      print(c(1,x))
      tv <- everything.df[grep(paste(ics, collapse = "|"), names(everything.df))][x]
      print(c("It's tv", tv))
      sapply(1:nrow(tv), function(i){
          #print(c("It's i", i))
          #print(c(2,i))
          #print(c(3,icsDF[apply(icsDF, 1, function(d) all(d[lst]== everything.df[i, lst]))]))
          #print(c(4,icsDF[apply(icsDF, 1, function(d) all(d[lst] == everything.df[i, lst]))][colnames(icsDF) == colnames(tv)]))
          if (is.numeric(icsDF[apply(icsDF, 1, function(d) all(d[lst] == everything.df[i, lst]))][colnames(icsDF) == colnames(tv)])){
          res <- tv[i,]-icsDF[apply(icsDF, 1, function(d) all(d[lst] == everything.df[i, lst]))][colnames(icsDF) == colnames(tv)]
          } else{res <- tv[i,]-0}
          res <- unlist(res)
          unname(res)
          #print(res)
          })
      })
  icsDF.calc <- as.data.frame(icsDF.calc)
  names(icsDF.calc) <- names(everything.df[grep(paste(ics, collapse = "|"), names(everything.df))])
  everything.df[grep(paste(ics, collapse = "|"), names(everything.df))] <- icsDF.calc
  everything.df[everything.df$ICs == 1,] <- icsDF
  everything.df[everything.df$Unstained == 1,] <- unstDF

  #Means processing
  meandf <- all_data[grep("Mean|SD", all_data$filename, invert = T),]
  meandf <- meandf[grep("Mean", names(meandf))]
  numLevels <- 2
  meanNames <- lapply(names(meandf), function(x) tail(strsplit(x, "/")[[1]], numLevels))
  meanNames <- lapply(meanNames, function(x) strsplit(x, "[|]"))
  channelNames <- lapply(meanNames, function(x) tail(tail(x, 1)[[1]], 1))
  channelNames <- unname(channelNames)
  channelNames <- lapply(channelNames, function(x) {
     x <- gsub("Mean| |-A)$|Comp-", "", x)
     x <- gsub("[(]","",x)
     x <- gsub("[-]", ".",x)
     })
  channelNames <- unlist(channelNames)
  meanNames <- lapply(meanNames, function(x) lapply(x, function(i) i[grep("Mean",i,invert = T)]))
  meanNames.full <- sapply(1:length(meanNames), function(x){
      a <- sapply(meanNames[[x]], nicename)
      paste(c(a[1:length(a)], channelNames[x], "MFI"), collapse = " ")
  })
  names(meandf) <- meanNames.full
  everything.df <- cbind(everything.df, meandf)
  return(everything.df)
  }
filePrep <- function(lnNum, read){
    if (read[lnNum] == "filename"){
        return(list("flname" = read[lnNum+1]))
    }
    if (read[lnNum] == "isotypes"){
        return(list("isotypes" = read[lnNum+1]))
    }
    if (read[lnNum] == "populations"){
        return(list("populations" = read[lnNum+1:length(read)]))
    }
}

