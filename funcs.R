blya <- function(all_data, cellcount, vcellcount, pop_ls, srch_v, kids_names, pop){
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
          #Nodes checking
        }
        if (grep(j, i) - grep(pop, i) > 1 & i[grep(j, i)] == tail(i, 1)){
              node <- i[grep(j, i)-1]
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
    if (flname == 0){
        return(NULL)
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
            }
}
popConstruct <- function(popul_file, all_data){
  cellcount <- names(all_data)[grep("Count", x = names(all_data))]
  ncellcount <- gsub(c(" Count"), "", cellcount)
  ncellcount <- gsub(" [|]", "", ncellcount)
  ncellcount <- gsub(" ", "", ncellcount, perl = TRUE)
  vcellcount <- strsplit(ncellcount, "/")
  inpops <- file("pops.txt", open = "r")
  finaldf <- NULL
  resname <- NULL
  forDF <- NULL
  while(TRUE){
    #print(resname)
    #l <- "Treg_CD45+ - RORgt+, Tbet+, ST2hi, ST2neg, Gata3+, CD25+"
    l <- readLines(inpops, n=1)
    if (l == "end"){
      cat("\nEnd of file")
      break
    }
    cat("========= Start of line ============\n")
    #print(l)
    rd <- strsplit(l, "_", fixed = TRUE)
    print(c("Population list", rd))
    pop <- rd[[1]][1]
    print(c("Main population name", pop))
    nrd <- gsub(" ", "", strsplit(rd[[1]][2],"-")[[1]])
    #print(nrd)
    par <- strsplit(gsub(" ", "", nrd[1]), ",")
    print(c("Parent populations name", par))
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
      print(c(length(kids), "Inherited populations", kids))
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
      #print(gsub(" ", "", nicename(names(all_data))) == i)
      if (length(grep(nicename(gsub(" ","",pop)), tail(strsplit(i, "/")[[1]], 1))) > 0){
        #print(length(grep(nicename(gsub(" ","",pop)), tail(strsplit(i, "/")[[1]], 1))) > 0)
        name <- c(name, paste(nicename(i),"target"))
        #print(paste(nicename(i),"target"))
      }else{
      name <- c(name, nicename(i))
      #print(nicename(i))
      }
    }
    #View(tempdf)
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
   # if (is.null(chknode(all_data, vcellcount, kids, kids_ls, real_kids_cols, pop)) == FALSE){
    #  nodesdf <- as.data.frame(out_nodes)
     # tempdf <- cbind(tempdf, nodesdf)
    #  names(tempdf) <- name
    #}
    #Sorry
   # wrong_parents <- blya(parents_ls, par, pop)
    #wrong_kids <- blya(kids_ls, par, pop)
    #if (is.null(wrong_kids) == FALSE | is.null(wrong_parents) == FALSE){
    #  cat(paste("Probably, something is wrong. Check the following data: \n", wrong_kids,"\n", wrong_parents))
    #}

    pcol <- grep("parent", names(tempdf))
    tcol <- grep("target", names(tempdf))
    resdf <- NULL
    for (i in pcol){
      resdf <- c(resdf,tempdf[,tcol]/tempdf[,i]*100)
     # print(paste("typeof(tempdf): ",typeof(tempdf)))
      #print(paste("typeof(resdf): ",typeof(resdf)))
      #cat(paste(resdf, "(resdf)\n"))
      tlab <- strsplit(tail(strsplit(names(tempdf)[tcol],"/")[[1]],1)," ")[[1]][1]
      plab <- strsplit(tail(strsplit(names(tempdf)[pcol],"/")[[1]],1)," ")[[1]][1]
      resname <- c(resname, paste(tlab,plab, sep="_"))
      print(paste("resname, parents cycle: ", resname))
    }
    bLYA <- blya(all_data,cellcount,vcellcount,kids_ls,kids,kids_names,pop)
    real_kids_cols <- bLYA[[1]]
    cat(paste("**********\n Real kids cols value\n", real_kids_cols, "**********"))
    print(typeof(resdf))
    for (i in real_kids_cols){
      print("kids cycle entered")
      resdf <- list(resdf,(tempdf[names(tempdf) == kids_names[i]])[[1]]/tempdf[,tcol]*100)
      print(typeof(resdf))
      tlab <- tail(strsplit(kids_names[i],"/")[[1]],1)
      plab <- strsplit(tail(strsplit(names(tempdf)[tcol],"/")[[1]],1)," ")[[1]][1]
      resname <- c(resname, paste(tlab,plab, sep="_"))
      print(paste("resname, kids cycle: ", resname))
    }
    resdf <- as.data.frame(resdf)
    print(typeof(resdf))
    #print(resdf)
    cat("\nCreated temporary Data frame\n")
    if ( length(line) == 0 ) {
      break
    }
    cat("========= End of line ============\n\n")

    finaldf <- c(finaldf, resdf)
    forDF <- c(forDF, bLYA[[2]])
  }
  finaldf <- as.data.frame(finaldf)
  names(finaldf) <- resname
  forDF <- as.data.frame(forDF)
  finaldf <- c(finaldf, forDF)
  finaldf <- as.data.frame(finaldf)
  finaldf$filename <- all_data$filename
  close(inpops)
  return(finaldf)
}
