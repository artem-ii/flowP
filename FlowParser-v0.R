
setwd() ###### SET THIS!

library(readxl)
library(ggplot2)
kidsNodes <- function(all_data, cellcount, vcellcount, pop_ls, srch_v, kids_names, pop, verbose){
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
filenameparser <- function(df, flname = "exp_ID date tissue genotype treatment sample_ID", verbose = F){
  if (flname == casefold("no")){
    return(df)
  }else{
    map <- strsplit(flname, " ")[[1]]
    col <- strsplit(df$filename, "_")
    col <- col[col != c("Mean", "SD")]
    col <- col[grep("Comp|Contr", col, ignore.case = T, invert = T)]
    col <- sapply(col, function(x) x[-length(x)])
    if (verbose == T) print(c("filenameparser, col ", col))
    v <- grep("ICs|Unstained|Isotype|unst|iso|IC", col, invert = TRUE, ignore.case = T)
    newdf <- as.data.frame(lapply(seq_along(map), function(i) sapply(col[v], "[[", i)))
    names(newdf) <- map
    newdf$ICs <- rep(0, nrow(newdf))
    newdf$Unstained <- rep(0, nrow(newdf))
    newdf1 <- as.data.frame(lapply(seq_along(col[-v][[1]]), function(i) sapply(col[-v], "[[", i)))
    if (verbose == T) View(newdf1)
    #Process ICs
    ic <- sapply(newdf1, function(x) grep("ICs|Isotype|iso|IC", x, ignore.case = T))
    if (verbose == T) print(c("filenameparser, ic ", ic))
    ic.logic <- sapply(ic[ic>0], function(x) is.numeric(x)[[1]])
    if (verbose == T) print(c("filenameparser, ic.logic ", ic.logic))
    ic <- ic[ic.logic][[1]]
    ic.col <- sapply(seq_along(ic.logic), function(x){ if (ic.logic[x] == TRUE){return(x)}})
    ic.col <- ic.col[sapply(ic.col, is.integer)][[1]]
    newdf.IC <- newdf1[ic,1:ic.col-1]
    newdf.IC[ic.col:length(map)] <- NA
    names(newdf.IC) <- map
    newdf.IC$ICs <- rep(1, nrow(newdf.IC))
    newdf.IC$Unstained <- rep(0, nrow(newdf.IC))
    newdf.IC <- cbind(newdf.IC, df[ic,])
    sensedf <- cbind(newdf, df[v,])
    sensedf <- rbind(newdf.IC, sensedf)
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
    newdf.UN <- cbind(newdf.UN, df[unst,])
    sensedf <- rbind(newdf.UN,sensedf)
    sensedf <- sensedf[,]
    return(sensedf)
  }
}
popConstruct <- function(popul_file, all_data, refSample = 0, refSampleRow = 0, verbose = FALSE){
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
    nrd <- gsub(" ", "", strsplit(rd[[1]][2],">")[[1]])
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
    kids.nodes <- kidsNodes(all_data,cellcount,vcellcount,kids_ls,kids,kids_names,pop, verbose)
    real_kids_cols <- kids.nodes[[1]]
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
    forDF <- kids.nodes[[2]]
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
  everything.df <- everything.df[grep("Mean|SD|Comp", everything.df$filename, invert = T),]
  everything.df <- filenameparser(everything.df, flname, verbose)
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
    if(verbose == T) {print(c(1,x))}
    tv <- everything.df[grep(paste(ics, collapse = "|"), names(everything.df))][x]
    if(verbose == T) {print(c("It's tv", tv))}
    sapply(1:nrow(tv), function(i){
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
  meandf <- all_data[grep("Mean|SD|Comp", all_data$filename, invert = T),]
  meandf <- meandf[grep("Mean", names(meandf))]
  numLevels <- 2
  meanNames <- lapply(names(meandf), function(x) tail(strsplit(x, "/")[[1]], numLevels))
  meanNames <- lapply(meanNames, function(x) strsplit(x, "[|]"))
  channelNames <- lapply(meanNames, function(x) tail(tail(x, 1)[[1]], 1))
  channelNames <- unname(channelNames)
  channelNames <- lapply(channelNames, function(x) {
    x <- gsub("Mean| |-A)$|Comp-", "", x)
    x <- gsub("[(]","",x)
    x <- gsub("[-]", "_",x)
  })
  channelNames <- unlist(channelNames)
  meanNames <- lapply(meanNames, function(x) lapply(x, function(i) i[grep("Mean",i,invert = T)]))
  meanNames.full <- sapply(1:length(meanNames), function(x){
    a <- sapply(meanNames[[x]], nicename)
    gsub(" ", "", paste(c(a[1:length(a)], channelNames[x], "MFI"), collapse = "_"))
  })
  channels <- param[sapply(1:length(param),function(i) names(param[i][[1]]) == "channels")]
  channels <- channels[[1]][[1]][!is.na(channels[[1]][[1]])]
  channels <- sapply(channels, strsplit, ",")
  channels <- sapply(channels, function(x) gsub(" ", "", x))
  channels.c <- sapply(channels, function(x) gsub("[-]","_", strsplit(x, "=")[[1]][1]))
  channels.a <- sapply(channels, function(x) strsplit(x, "=")[[1]][2])
  meanlst <- sapply(1:length(channels.c), function(y){
    a <- meanNames.full[grep(channels.c[[y]], meanNames.full)]
    a <- gsub(channels.c[[y]], channels.a[[y]], unlist(a))
    meanNames.full[grep(channels.c[[y]], meanNames.full)] <- unlist(a)
    return(meanNames.full[grep(channels.a[[y]], meanNames.full)])
    })
  seq.1 <- sort(sapply(1:length(channels.c), function(x) nchar(channels.c[x])),decreasing = T)
  for (i in 1:length(seq.1)){
    p <- meanlst[names(channels.c) == names(seq.1[i])][[1]]
    meanNames.full[grep(channels.c[names(channels.c) == names(seq.1[i])], meanNames.full)] <- p[grep(paste(c(channels.a[names(channels.c) == names(seq.1[i])],"MFI"),collapse = "_"), p)]
  }
  names(meandf) <- meanNames.full
  c.a <- everything.df$ICs == 0
  c.b <- everything.df$Unstained == 0
  meandf.s <- meandf[c.a & c.b,]
  pos <- sapply(meandf.s, function(x){
    i <- x>0
    i[is.na(i)] <- TRUE
    return(i)})
  meandf[everything.df$ICs == 0 & everything.df$Unstained == 0,] <- as.data.frame(sapply(1:ncol(meandf.s), function(x){
    if(!all(pos[,x])){
      m <- min(meandf.s[x])
      meandf[everything.df$ICs == 0 & everything.df$Unstained == 0, x][[1]] <- ((meandf[everything.df$ICs == 0 & everything.df$Unstained == 0, x]+1)-m)[[1]]
      return(meandf[everything.df$ICs == 0 & everything.df$Unstained == 0, x])}else{
        return(meandf[everything.df$ICs == 0 & everything.df$Unstained == 0, x])
      }}))
  
  if (refSample != 0){
    refSampleRow <- sapply(everything.df[names(everything.df) == strsplit(flname, " ")[[1]]],function(x){
      i <- grep(refSample, x)
      return(i)})
    refSampleRow <- refSampleRow[sapply(refSampleRow, function(x) length(x[1])>0 & x[1]>0)]
    x <- refSampleRow[!is.na(names(refSampleRow))][[1]]
    refSampleRow <- meandf[x,]
    refSampleRow <<- meandf[x,]}
  meandf[c.a & c.b,] <- as.data.frame(sapply(1:ncol(meandf.s), function(x) meandf[c.a & c.b, x]/refSampleRow[,x][[1]]*100))
  #Prepare data
  everything.df <- cbind(everything.df, meandf)
  everything.df$filename <- factor(everything.df$filename)
  rownames(everything.df) <- NULL
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
  if (read[lnNum] == "channels"){
    return(list("channels" = read[lnNum+1]))
  }
}
compute <- function(df, iv = 0, flowP = TRUE, verbose = FALSE, nomfi = FALSE){
  iv <- substitute(iv)
  if(flowP==T){
    df1 <- df[df$ICs == 0 & df$Unstained == 0,]
    df1 <- df1[names(df1) != "ICs" & names(df1) != "Unstained"]
    if (nomfi == T) df1 <- df1[grep("MFI", names(df1), invert = T)]
  }else {df1 <- df}
  numdf <- df1[,sapply(df1,is.double)]
  if (all(sapply(df1$filename, function(x) grepl("_",x)))) {
    l <- strsplit(levels(df1$filename),"_")[[1]][1:2]}else{l <- df1$filename[[1]]}
  fl <- paste(c(l, "Stats.csv"),collapse = "_")
  file.create(paste(c(l, "Stats.csv"),collapse = "_"))
  if (verbose == T) print(c("Output file name", fl))
  outFile <- file(fl, open = "w")
  if (is.numeric(iv) & iv == 0){
    warning("Automatic mode doesn't work yet.\nPlease specify an independant variables formula as \"iv\" variable in compute()")
    t <- grep("tissue|organ", casefold(names(df1)))
    c <- grep("cond|treat", casefold(names(df1)))
    g <- grep("geno", casefold(names(df1)))
    if(length(t)>0) df1[,t] <- factor(df1[,t])
    if(length(c)>0) df1[,c] <- factor(df1[,c])
    if(length(g)>0) df1[,g] <- factor(df1[,g])
    df1$filename <- factor(df1$filename)
    lst <- c(t,c,g)
    lst <- lst[sapply(seq_along(lst), function(x) length(levels(df1[,lst[x]])) > 1)]
    lst <- names(df1[,lst])
    stat <- sapply(1:ncol(numdf), function(x){
      p <- names(numdf[x])
      sapply(lst, function(i) {
        anova <- aov(df1[[x]] ~ df1[[i]], df1)
        if(verbose == T) print(c("Population ", p))
        if(verbose == T) print(c("Independant variable ", names(df1[i])))
        if(verbose == T) print(summary(a))
        #if(verbose == T)
        a <- summary(anova)
        d <- sapply(a, function(x){
          x <- as.data.frame(x)
          dv_col <- rep(dv,nrow(x))
          dv_col <- as.data.frame(dv_col)
          colnames(dv_col) <- "Population"
          x <- cbind(dv_col, x)
          if (verbose == T) print(x)
          if (verbose == T) print(typeof(x))
          write.csv2(x, outFile)
        })
        b <- TukeyHSD(anova)
        c <- sapply(b, function(x){
          if (verbose == T) print(x)
          x <- as.data.frame(x)
          dv_col <- rep(dv,nrow(x))
          dv_col <- as.data.frame(dv_col)
          colnames(dv_col) <- "Population"
          x <- cbind(dv_col, x)
          if (verbose == T) print(x)
          if (verbose == T) print(typeof(x))
          write.csv2(x, outFile)
        })
        x <- c(d,c)
        x <- as.data.frame(x)
        write.csv2(x, outFile)
      })
    })
    close(outFile)
  }else{
    f <- iv
    if (verbose == T) print(c("f ",f))
    if (verbose == T) print(typeof(f))
    stat <- sapply(1:ncol(numdf), function(x){
      dv <- names(numdf[x])
      if (verbose == T) print(c("dv ", dv))
      if (verbose == T) print(c("f ",f))
      f <- reformulate(as.character(f),dv)
      anova <- aov(f,df1)
      a <- summary(anova)
      d <- sapply(a, function(x){
        x <- as.data.frame(x)
        dv_col <- rep(dv,nrow(x))
        dv_col <- as.data.frame(dv_col)
        colnames(dv_col) <- "Population"
        x <- cbind(dv_col, x)
        if (verbose == T) print(x)
        if (verbose == T) print(typeof(x))
        write.csv2(x, outFile)
      })
      b <- TukeyHSD(anova)
      c <- sapply(b, function(x){
        if (verbose == T) print(x)
        x <- as.data.frame(x)
        dv_col <- rep(dv,nrow(x))
        dv_col <- as.data.frame(dv_col)
        colnames(dv_col) <- "Population"
        x <- cbind(dv_col, x)
        if (verbose == T) print(x)
        if (verbose == T) print(typeof(x))
        write.csv2(x, outFile)
      })
    })
    
    close(outFile)
  }
  return(paste(c("Output file ",fl, "is written to working directory", getwd()), collapse = " "))
}
plot_all <- function(df, ind_var1, facet_var = 0, flowP = TRUE, verbose = FALSE){
  require(ggplot2)
  olddir <- getwd()
  dir.create(paste(c(olddir,"plots"), collapse = "/"))
  setwd(paste(c(olddir,"plots"), collapse = "/"))
  if (flowP == TRUE){
    df <- df[df$ICs != 1 & df$Unstained != 1,]
    df <- df[names(df) != "ICs" & names(df) != "Unstained"]
  }
  numdf <- df[,sapply(df,is.double)]
  n <- colnames(numdf)
  print(c("Plots will be saved to ", getwd()))
  if (facet_var == 0){
    p <- sapply(1:ncol(numdf), function(x){
      a <- ggplot(df, aes_string(ind_var1, n[x])) +
        scale_color_manual(values = c("WT" = "darkred", "KO" = "deepskyblue4"))+
        geom_point(position = position_jitter(0.1))+
        labs(title = " ", y = gsub("_", " % of ",n[x]))+
        geom_boxplot(width=0.35,size=0.35, alpha=0.1)+
        #facet_wrap(~Organ, ncol=2)+
        #geom_path(data = data4, aes(x = x, y = y))+
        theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_blank(),plot.title = element_text(hjust = 0.5, face = "bold"),axis.ticks = element_line(linetype = "solid"), plot.background = element_rect(fill = NA),axis.text.y = element_text(face = "bold", size = 12),axis.title.x = element_blank(),axis.title.y = element_text(face = "bold", size = rel(1.3), angle = 90),axis.line = element_line(size = 1.5), panel.background = element_blank(),axis.text.x = element_text(size = 12, face = "bold", colour = "black"),legend.position = "none")+
        #annotate("text", x = 1.5, y = 8.5, label = "p=0.043", fontface = 2)+
        scale_x_discrete(limits = c("wt","ko"))
      ggsave(plot = a,filename=paste(n[x], ".svg",sep=""),width=10, height=8)
      if (verbose == T){
        print(a)
        return(a)
      }else{return(0)}
    })
  }else{
    p <- sapply(1:ncol(numdf), function(x){
      a <- ggplot(df, aes_string(ind_var1, n[x])) +
        geom_point(position = position_jitter(0.1))+
        labs(title = " ", y = gsub("_", " % of ",n[x]))+
        geom_boxplot(width=0.35,size=0.35, alpha=0.1)+
        facet_wrap(facet_var, ncol=2)+
        #geom_path(data = data4, aes(x = x, y = y))+
        theme(strip.text = element_text(face = "bold", size = 12), strip.background = element_blank(),plot.title = element_text(hjust = 0.5, face = "bold"),axis.ticks = element_line(linetype = "solid"), plot.background = element_rect(fill = NA),axis.text.y = element_text(face = "bold", size = 12),axis.title.x = element_blank(),axis.title.y = element_text(face = "bold", size = rel(1.3), angle = 90),axis.line = element_line(size = 1.5), panel.background = element_blank(),axis.text.x = element_text(size = 12, face = "bold", colour = "black"),legend.position = "none")+
        #annotate("text", x = 1.5, y = 8.5, label = "p=0.043", fontface = 2)+
        scale_x_discrete(limits = c("wt","ko"))
      ggsave(plot = a,filename=paste(n[x], ".svg",sep=""),width=10, height=8)
      if (verbose == T){
        print(a)
        return(a)
      }else{return(0)}
    })}
  setwd(olddir)
}
file <- "/Volumes/BigHDD/Inserm/Exps/RoraFoxp3_lung1/14-09-17/SPLEEN(Lung1)-14-09-17.xls"
df1 <- read_excel(file,
                       col_types = c("text", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric",
                                     "numeric", "numeric", "numeric"))
names(df1)[1] <- "filename"
names(df1) <- gsub("CD4 ", "CD4+", names(df1))
names(df1) <- gsub("CD8 ", "CD8+", names(df1))
names(df1) <- gsub("T-bet", "Tbet", names(df1))

df1.df <- popConstruct("pops.txt", df1, "m302")

df2 <- read_excel("/Volumes/BigHDD/Inserm/Exps/RoraFoxp3_lung1/14-09-17/LUNG(Lung1)_14-09-17.xls",
                   col_types = c("text", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric", "numeric",
                                 "numeric", "numeric"), na = "NA")
names(df2)[1] <- "filename"
names(df2) <- gsub("CD4 ", "CD4+", names(df2))
names(df2) <- gsub("CD8 ", "CD8+", names(df2))
names(df2) <- gsub("T-bet", "Tbet", names(df2))
df2.df <- popConstruct("pops.txt", df2, refSample = 0, refSampleRow)

df <- rbind(df1.df, df2.df)

compute(df, "Genotype*Treatment*Tissue")
plot_all(df, "Genotype", Tissue~Treatment)



#Digging syntax
#func <- function(a, var1 = 0) {
#    print(a)
#    x <- substitute(var1)
#    if(is.numeric(x) & x == 0){
#        print(is.numeric(x))
#        print(x == 0)
#    }
#    #a <- eval(substitute(alist(var1)))
#    #b <- eval(substitute(alist(var2)))
#    #print(a)
#    #print(b)
#    print(x)
#    #a <- eval(substitute(alist(...)))
#    }

#formula(1~2+3)
#a <- aov(formula(df1[[x]] ~ paste(names(df1[grep("Treatment|Tissue",names(df1))], collapse = ""))), df1)

###lmwrap <- function(d,y,x) {
#    ys <- deparse(substitute(y))
#    xs <- deparse(substitute(x))
#    f <- reformulate(xs,response=ys)
#    return(f)
#}
#mydata <- data.frame(X=1:10,Y=rnorm(10))
#a <- lmwrap(mydata,Y, X)
#reformulate("Genotype", names(df1[x]))
