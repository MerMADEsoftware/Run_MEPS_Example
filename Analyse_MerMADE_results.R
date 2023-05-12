###CREATING CONNECTIVITY MATRICES

#libraries needed
#library(igraph) #this package analyses network dynamics, which is exactly what we need!

#arguments used
#npatches : how many patches do you have in your seascape?
#out_conn : this is a TRUE/FALSE question, do you want to output the connectivity matrix?
#in_conn : a connectivity matrix to analyse, this is used if you want to read in a connectivity matrix and skip the creation steps
#out_outdeg : this is a TRUE/FALSE question, do you want to output the out degree information for each patch?
#out_indeg : this is a TRUE/FALSE question, do you want to output the in degree information for each patch?
#reps: how many replicates are in the file?
#NOTE: you can only choose one type of output, either connectivity matrix, out degree or in degree centrality.
#NOTE: if you want to speed up the process, you can run the function once to create the connectivity matrix, then
  #read it in as an argument and skip the creation section to then output the out or in degree centrality

conn_mat_reps<- function(npatches, out_conn, in_conn, out_outdeg, out_indeg, reps, simnum, peryear, inc_local){
  #this will hold the connectivity matrix
  connmat<- matrix(0, nrow=npatches, ncol=npatches)
  colnames(connmat)=seq(0, npatches-1,1)
  rownames(connmat)=seq(0, npatches-1,1)
  
  out_deg<-rep(0,npatches)
  in_deg<-rep(0,npatches)
  
  print("reading in file")
  indivs<- read.table(paste(simnum,"/Outputs/Individuals.txt", sep=""), header=T, as.is=T) #read in the file
  indivs_0<- indivs[indivs$Status!=0,] #status=0 is only original generations, remove
  indivs_sett<- indivs_0[!is.na(indivs_0$Ending_patch),] #we want only individuals that settled successfully
  if(inc_local==F){indivs_sett<- indivs_sett[indivs_sett$DistMoved>0,]}
  
  for(re in 0:(reps-1)){#for each replicate
    print(paste("doing replicate ", re, sep=""))
    which_rep<- indivs_sett[indivs_sett$Rep==re,]
    indivs_rep<- which_rep %>% distinct(Ind_ID, .keep_all = TRUE) #removing duplicates
    rep_con<- matrix(0, nrow=npatches, ncol=npatches)
    rep_out<- matrix(0, nrow=npatches, ncol=1)
    rep_in<- matrix(0, nrow=npatches, ncol=1)
    
    #i want the average PER YEAR connectivity matrix and centralities
    if(peryear==T){
      allyr_con=matrix(0, nrow=npatches, ncol=npatches)
      yr_out=matrix(0, nrow=npatches, ncol=1)
      yr_in=matrix(0, nrow=npatches, ncol=1)
      
      for(y in unique(indivs_rep$Year)){ #for every year
        yr_con<- matrix(0, nrow=npatches, ncol=npatches)
        yr<- indivs_rep[indivs_rep$Year==y,]
        for(r in 1:nrow(yr)){
          start<- yr$Starting_patch[r]
          end<- yr$Ending_patch[r]
          yr_con[start+1,end+1]=yr_con[start+1,end+1]+1 #because the patch names start at 0 and the row numbers start at 0
        } #fill the connectivity matrix
        #adjacency matrix for that year
        g1<- graph_from_adjacency_matrix(yr_con, mode=c("directed"), diag=F,
                                         weighted=T)
        #out centrality
        out_scores<- centr_degree(g1,mode="out")$res #get the out_degree centrality scores from the adjacency graph
        #in centrality
        in_scores<- centr_degree(g1,mode="in")$res #get the out_degree centrality scores from the adjacency graph
        
        #add them to the year totals
        if(y==0){
          yr_out=out_scores
          yr_in=in_scores
        }
        else{
          yr_out=cbind(yr_out, out_scores)
          yr_in=cbind(yr_in, in_scores)
        }
        allyr_con=allyr_con+yr_con
      }
      #take averages of the years for each rep
      rep_con=allyr_con/length(unique(indivs_sett$Year))
      rep_out=rowMeans(yr_out)
      rep_in=rowMeans(yr_in)
      
    }
    else{
      for(r in 1:nrow(indivs_rep)){ #for each individual in that replicate
        start<- indivs_rep$Starting_patch[r] #figure out where it started
        end<- indivs_rep$Ending_patch[r] #and where it ended
        #add 1 to the correct row and column in the connectivity matrix
        rep_con[start+1,end+1]=rep_con[start+1,end+1]+1 #because the patch names start at 0 and the row numbers start at 1, we need to do start+1 and end+1
      }
      g1<- graph_from_adjacency_matrix(rep_con, mode=c("directed"), diag=F,
                                       weighted=T)
      scores<- centr_degree(g1,mode="out")$res #get the out_degree centrality scores from the adjacency graph
      rep_out=scores #add it to the out_degree matrix
      scores<- centr_degree(g1,mode="in")$res #get the out_degree centrality scores from the adjacency graph
      rep_in=scores #add it to the in_degree matrix
      
    }
    connmat=connmat+rep_con
    if(re==0){
      out_deg=rep_out
      in_deg=rep_in
    }
    else{
      out_deg=cbind(out_deg, rep_out)
      in_deg=cbind(in_deg, rep_in)
    }
  }
  connmat=connmat/reps
  out_degree=rowMeans(out_deg)
  in_degree=rowMeans(in_deg)
  
  colnames(out_deg)=seq(0,reps-1,1)
  rownames(out_deg)=seq(0,npatches-1,1)
  colnames(in_deg)=seq(0,reps-1,1)
  rownames(in_deg)=seq(0,npatches-1,1)
  
  if(out_conn==T){return (connmat)}
  else if(out_outdeg==T){
    centralities<- cbind(seq(0,npatches-1,1),rowMeans(out_deg), rowMeans(in_deg))
    return(centralities)
  }
}

##################################

###CREATE A TABLE OF DISPERSAL METRICS
#this will tell you the fate of your dispersers

#arguments used
#rep : which replicate do you want to look at?
#year : which year do you want to look at?

disp_mets_reps<- function(nreps, year, simnum){
  metrics<- matrix(0, nrow=1, ncol=7)
  colnames(metrics)<- c("Mortality","Successful Dispersers", "Forced settlement: suitable", "Buffer capture", "Forced settlement: unsuitable","Mean distance (km)", "Max distance (km)")
  
  indivs<- read.table(paste(simnum, "/Outputs/Individuals.txt", sep=""), header=T, row.names=NULL, as.is=T)
  for(r in 0:(nreps-1)){
    
    
    indivs_rep<- indivs[indivs$Rep==rep,]
    yr<- indivs_rep[indivs_rep$Year<=year,]
    
    disp<- yr[yr$DistMoved!=0,] #use only dispersers
    
    #look at mortality
    i_mort<- disp[disp$Status==5,] #status=5 means they died due to dispersal related mortality
    perc_mort<- nrow(i_mort)/nrow(disp) #what percentage of dispersers met their end this way?
    
    #successful dispersers
    i_succ<- disp[disp$Status==3,] #chose to settle
    perc_succ<- nrow(i_succ)/nrow(disp)
    i_forced_succ<- disp[disp$Status==9,] #were forced to settle, got lucky with suitable habitat (passive disperser)
    perc_forced_suit<- nrow(i_forced_succ)/nrow(disp)
    #unsuccessful
    i_buff<- disp[disp$Status==12,] #settled because they were captured by buffer
    perc_buff<- nrow(i_buff)/nrow(disp)
    i_forced_die<- disp[disp$Status==4,] #were forced to settle, died on unsuitable habitat (passive disperser)
    perc_forced_unsuit<- nrow(i_forced_die)/nrow(disp)
    i_time<- disp[disp$Status==10,] #they ran out of time
    perc_time<- nrow(i_time)/nrow(disp)
    i_abs<- disp[disp$Status==8,] #they were absorbed by the seascape boundary and lost
    perc_abs<- nrow(i_abs)/nrow(disp)
    
    #to look at the distance that successful dispersers travelled (status=3)
    # mean_dist<- mean(i_succ$DistMoved, na.rm=T)
    mean_dist<- mean(disp$DistMoved, na.rm=T)
    # max_dist=max(i_succ$DistMoved, na.rm=T)
    max_dist=max(disp$DistMoved, na.rm=T)

    metrics[1,1]<- metrics[1,1] + perc_mort + perc_time + perc_abs #mortality is a sum of dispersal mortality, running out of time and being absorbed
    metrics[1,2]<- metrics[1,2]+ perc_succ
    metrics[1,3]<- metrics[1,3]+ perc_forced_suit
    metrics[1,4]<- metrics[1,4]+ perc_buff
    metrics[1,5]<- metrics[1,5]+ perc_forced_unsuit
    metrics[1,6]<- metrics[1,6]+ (mean_dist /1000) #turn into km, right now they are in metres
    metrics[1,7]<- metrics[1,7]+ (max_dist / 1000) #convert to km
  }
  
  metrics[1,1]<- metrics[1,1]/nreps
  metrics[1,2]<- metrics[1,2]/nreps
  metrics[1,3]<- metrics[1,3]/nreps
  metrics[1,4]<- metrics[1,4]/nreps
  metrics[1,5]<- metrics[1,5]/nreps
  metrics[1,6]<- metrics[1,6]/nreps
  metrics[1,7]<- metrics[1,7]/nreps
  
  colnames(metrics)<- c("Mortality","Successful Dispersers", "Forced settlement: suitable", "Buffer_capture", "Forced settlement: unsuitable","Mean distance (km)", "Max distance (km)")
  
  return(metrics)
}

#######################################################
#producing a heatmap of the connectivity matrix
#this function is edited from the gplots package function heatmap2()

new_heatmap2<- 
  function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
            distfun = dist, hclustfun = hclust, dendrogram = c("both", 
                                                               "row", "column", "none"), reorderfun = function(d, w) reorder(d, 
                                                                                                                             w), symm = FALSE, scale = c("none", "row", "column"), 
            na.rm = TRUE, revC = identical(Colv, "Rowv"), add.expr, breaks, 
            symbreaks = any(x < 0, na.rm = TRUE) || scale != "none", 
            col = "heat.colors", colsep, rowsep, sepcolor = "white", 
            sepwidth = c(0.05, 0.05), cellnote, notecex = 1, notecol = "cyan", 
            na.color = par("bg"), trace = c("column", "row", "both", 
                                            "none"), tracecol = "cyan", hline = median(breaks), vline = median(breaks), 
            linecol = tracecol, margins = c(5, 5), ColSideColors, RowSideColors, 
            cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
            labCol = NULL, srtRow = NULL, srtCol = NULL, adjRow = c(0, 
                                                                    NA), adjCol = c(NA, 0), offsetRow = 0.5, offsetCol = 0.5, 
            colRow = NULL, colCol = NULL, key = TRUE, keysize = 1.5, 
            density.info = c("histogram", "density", "none"), denscol = tracecol, 
            symkey = any(x < 0, na.rm = TRUE) || symbreaks, densadj = 0.25, 
            key.title = NULL, key.xlab = NULL, key.ylab = NULL, key.xtickfun = NULL, 
            key.ytickfun = NULL, key.par = list(), main = NULL, xlab = NULL, 
            ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL, extrafun = NULL, 
            ...) {
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    retval <- list()
    scale <- if (symm && missing(scale)) 
      "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
      col <- get(col, mode = "function")
    if (!missing(breaks) && any(duplicated(breaks))) 
      stop("breaks may not contain duplicate values")
    if (!missing(breaks) && (scale != "none")) 
      warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
              "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || any(is.na(Rowv))) 
      Rowv <- FALSE
    if (is.null(Colv) || any(is.na(Colv))) 
      Colv <- FALSE
    else if (all(Colv == "Rowv")) 
      Colv <- Rowv
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
      stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
      stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
      stop("`margins' must be a numeric vector of length 2")
    if (missing(cellnote)) 
      cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
      if (((is.logical(Rowv) && !isTRUE(Rowv)) || (is.null(Rowv))) && 
          (dendrogram %in% c("both", "row"))) {
        warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        if (dendrogram == "both") 
          dendrogram <- "column"
        else dendrogram <- "none"
      }
    }
    if (!inherits(Colv, "dendrogram")) {
      if (((is.logical(Colv) && !isTRUE(Colv)) || (is.null(Colv))) && 
          (dendrogram %in% c("both", "column"))) {
        warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        if (dendrogram == "both") 
          dendrogram <- "row"
        else dendrogram <- "none"
      }
    }
    if (inherits(Rowv, "dendrogram")) {
      ddr <- Rowv
      rowInd <- order.dendrogram(ddr)
      if (length(rowInd) > nr || any(rowInd < 1 | rowInd > 
                                     nr)) 
        stop("Rowv dendrogram doesn't match size of x")
      if (length(rowInd) < nr) 
        nr <- length(rowInd)
    }
    else if (is.integer(Rowv)) {
      distr <- distfun(x)
      hcr <- hclustfun(distr)
      ddr <- as.dendrogram(hcr)
      ddr <- reorderfun(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd)) 
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Rowv)) {
      Rowv <- rowMeans(x, na.rm = na.rm)
      distr <- distfun(x)
      hcr <- hclustfun(distr)
      ddr <- as.dendrogram(hcr)
      ddr <- reorderfun(ddr, Rowv)
      rowInd <- order.dendrogram(ddr)
      if (nr != length(rowInd)) 
        stop("row dendrogram ordering gave index of wrong length")
    }
    else if (!isTRUE(Rowv)) {
      rowInd <- nr:1
      ddr <- as.dendrogram(hclust(dist(diag(nr))))
    }
    else {
      rowInd <- nr:1
      ddr <- as.dendrogram(Rowv)
    }
    if (inherits(Colv, "dendrogram")) {
      ddc <- Colv
      colInd <- order.dendrogram(ddc)
      if (length(colInd) > nc || any(colInd < 1 | colInd > 
                                     nc)) 
        stop("Colv dendrogram doesn't match size of x")
      if (length(colInd) < nc) 
        nc <- length(colInd)
    }
    else if (identical(Colv, "Rowv")) {
      if (nr != nc) 
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      if (exists("ddr")) {
        ddc <- ddr
        colInd <- order.dendrogram(ddc)
      }
      else colInd <- rowInd
    }
    else if (is.integer(Colv)) {
      distc <- distfun(if (symm) 
        x
        else t(x))
      hcc <- hclustfun(distc)
      ddc <- as.dendrogram(hcc)
      ddc <- reorderfun(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd)) 
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (isTRUE(Colv)) {
      Colv <- colMeans(x, na.rm = na.rm)
      distc <- distfun(if (symm) 
        x
        else t(x))
      hcc <- hclustfun(distc)
      ddc <- as.dendrogram(hcc)
      ddc <- reorderfun(ddc, Colv)
      colInd <- order.dendrogram(ddc)
      if (nc != length(colInd)) 
        stop("column dendrogram ordering gave index of wrong length")
    }
    else if (!isTRUE(Colv)) {
      colInd <- 1:nc
      ddc <- as.dendrogram(hclust(dist(diag(nc))))
    }
    else {
      colInd <- 1:nc
      ddc <- as.dendrogram(Colv)
    }
    retval$rowInd <- rowInd
    retval$colInd <- colInd
    retval$call <- match.call()
    x <- x[rowInd, colInd]
    x.unscaled <- x
    cellnote <- cellnote[rowInd, colInd]
    if (is.null(labRow)) 
      labRow <- if (is.null(rownames(x))) 
        (1:nr)[rowInd]
    else rownames(x)
    else labRow <- labRow[rowInd]
    if (is.null(labCol)) 
      labCol <- if (is.null(colnames(x))) 
        (1:nc)[colInd]
    else colnames(x)
    else labCol <- labCol[colInd]
    if (!is.null(colRow)) 
      colRow <- colRow[rowInd]
    if (!is.null(colCol)) 
      colCol <- colCol[colInd]
    if (scale == "row") {
      retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
      x <- sweep(x, 1, rm)
      retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
      x <- sweep(x, 2, rm)
      retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
      if (missing(col) || is.function(col)) 
        breaks <- 16
      else breaks <- length(col) + 1
    }
    if (length(breaks) == 1) {
      if (!symbreaks) 
        breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm), 
                      length = breaks)
      else {
        extreme <- max(abs(x), na.rm = TRUE)
        breaks <- seq(-extreme, extreme, length = breaks)
      }
    }
    nbr <- length(breaks)
    ncol <- length(breaks) - 1
    if (is(col, "function")) 
      col <- col(ncol)
    min.breaks <- min(breaks)
    max.breaks <- max(breaks)
    x[x < min.breaks] <- min.breaks
    x[x > max.breaks] <- max.breaks
    if (missing(lhei) || is.null(lhei)) 
      lhei <- c(keysize, 4)
    if (missing(lwid) || is.null(lwid)) 
      lwid <- c(keysize, 4)
    if (missing(lmat) || is.null(lmat)) {
      lmat <- rbind(4:3, 2:1)
      if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors) != 
            nc) 
          stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
                        1)
        lhei <- c(lhei[1], 0.2, lhei[2])
      }
      if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors) != 
            nr) 
          stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                                             1), 1), lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
      }
      lmat[is.na(lmat)] <- 0
    }
    if (length(lhei) != nrow(lmat)) 
      stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
    if (length(lwid) != ncol(lmat)) 
      stop("lwid must have length = ncol(lmat) =", ncol(lmat))
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    plot.index <- 1
    if (!missing(RowSideColors)) {
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
      plot.index <- plot.index + 1
    }
    if (!missing(ColSideColors)) {
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
      plot.index <- plot.index + 1
    }
    par(mar = c(margins[1], 5, 0, margins[2])) #changed from par(mar = c(margins[1], 0, 0, margins[2]))
    x <- t(x)
    cellnote <- t(cellnote)
    if (revC) {
      iy <- nr:1
      if (exists("ddr")) 
        ddr <- rev(ddr)
      x <- x[, iy]
      cellnote <- cellnote[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
            c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
          breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
      retval$rowDendrogram <- ddr
    if (exists("ddc")) 
      retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!invalid(na.color) & any(is.na(x))) {
      mmat <- ifelse(is.na(x), 1, NA)
      image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
            col = na.color, add = TRUE)
    }
    if (is.null(srtCol) && is.null(colCol)) 
      axis(1, 1:nc, labels = labCol, las = 2, line = -0.5 + 
             offsetCol, tick = 0, cex.axis = cexCol, hadj = adjCol[1], 
           padj = adjCol[2])
    else {
      if (is.null(srtCol) || is.numeric(srtCol)) {
        if (missing(adjCol) || is.null(adjCol)) 
          adjCol = c(1, NA)
        if (is.null(srtCol)) 
          srtCol <- 90
        xpd.orig <- par("xpd")
        par(xpd = NA)
        xpos <- axis(1, 1:nc, labels = rep("", nc), las = 2, 
                     tick = 0)
        text(x = xpos, y = par("usr")[3] - (1 + offsetCol) * 
               strheight("M"), labels = labCol, adj = adjCol, 
             cex = cexCol, srt = srtCol, col = colCol)
        par(xpd = xpd.orig)
      }
      else warning("Invalid value for srtCol ignored.")
    }
    if (is.null(srtRow) && is.null(colRow)) {
      axis(2, iy, labels = labRow, las = 2, line = -0.5 + offsetRow, 
           tick = 0, cex.axis = cexRow, hadj = adjRow[1], padj = adjRow[2]) #changed from axis(4,...)
    }
    else {
      if (is.null(srtRow) || is.numeric(srtRow)) {
        xpd.orig <- par("xpd")
        par(xpd = NA)
        ypos <- axis(4, iy, labels = rep("", nr), las = 2, 
                     line = -0.5, tick = 0)
        text(x = par("usr")[2] + (1 + offsetRow) * strwidth("M"), 
             y = ypos, labels = labRow, adj = adjRow, cex = cexRow, 
             srt = srtRow, col = colRow)
        par(xpd = xpd.orig)
      }
      else warning("Invalid value for srtRow ignored.")
    }
    if (!is.null(xlab)) 
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    if (!is.null(ylab)) 
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
      eval(substitute(add.expr))
    if (!missing(colsep)) 
      for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, 
                                xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 
                                  1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
      for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
                                                        1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
                                                                                                         1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
                                col = sepcolor, border = sepcolor)
    min.scale <- min(breaks)
    max.scale <- max(breaks)
    x.scaled <- scale01(t(x), min.scale, max.scale)
    if (trace %in% c("both", "column")) {
      retval$vline <- vline
      vline.vals <- scale01(vline, min.scale, max.scale)
      for (i in 1:length(colInd)) {
        if (!is.null(vline)) {
          abline(v = i - 0.5 + vline.vals, col = linecol, 
                 lty = 2)
        }
        xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
        xv <- c(xv[1], xv)
        yv <- 1:length(xv) - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (trace %in% c("both", "row")) {
      retval$hline <- hline
      hline.vals <- scale01(hline, min.scale, max.scale)
      for (i in 1:length(rowInd)) {
        if (!is.null(hline)) {
          abline(h = i - 0.5 + hline.vals, col = linecol, 
                 lty = 2)
        }
        yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
        yv <- rev(c(yv[1], yv))
        xv <- length(yv):1 - 0.5
        lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
      }
    }
    if (!missing(cellnote)) 
      text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
           col = notecol, cex = notecex)
    plot.index <- plot.index + 1
    par(mar = c(margins[1], 0, 0, 0))
    if (dendrogram %in% c("both", "row")) {
      flag <- try(plot.dendrogram(ddr, horiz = TRUE, axes = FALSE, 
                                  yaxs = "i", leaflab = "none"))
      if ("try-error" %in% class(flag)) {
        cond <- attr(flag, "condition")
        if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
          stop("Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
      }
    }
    else plot.new()
    par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
    if (dendrogram %in% c("both", "column")) {
      flag <- try(plot.dendrogram(ddc, axes = FALSE, xaxs = "i", 
                                  leaflab = "none"))
      if ("try-error" %in% class(flag)) {
        cond <- attr(flag, "condition")
        if (!is.null(cond) && conditionMessage(cond) == "evaluation nested too deeply: infinite recursion / options(expressions=)?") 
          stop("Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option(\"expressions\"=...).")
      }
    }
    else plot.new()
    if (!is.null(main)) 
      title(main, cex.main = 1.5 * op[["cex.main"]])
    if (key) {
      mar <- c(5, 4, 2, 1)
      if (!is.null(key.xlab) && is.na(key.xlab)) 
        mar[1] <- 2
      if (!is.null(key.ylab) && is.na(key.ylab)) 
        mar[2] <- 2
      if (!is.null(key.title) && is.na(key.title)) 
        mar[3] <- 1
      par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
      if (length(key.par) > 0) 
        do.call(par, key.par)
      tmpbreaks <- breaks
      if (symkey) {
        max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
        min.raw <- -max.raw
        tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
        tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
      }
      else {
        min.raw <- min.breaks
        max.raw <- max.breaks
      }
      z <- seq(min.raw, max.raw, by = min(diff(breaks)/100))
      image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks, 
            xaxt = "n", yaxt = "n")
      par(usr = c(0, 1, 0, 1))
      if (is.null(key.xtickfun)) {
        lv <- pretty(breaks)
        xv <- scale01(as.numeric(lv), min.raw, max.raw)
        xargs <- list(at = xv, labels = lv)
      }
      else {
        xargs <- key.xtickfun()
      }
      xargs$side <- 1
      do.call(axis, xargs)
      if (is.null(key.xlab)) {
        if (scale == "row") 
          key.xlab <- "Row Z-Score"
        else if (scale == "column") 
          key.xlab <- "Column Z-Score"
        else key.xlab <- "Value"
      }
      if (!is.na(key.xlab)) {
        mtext(side = 1, key.xlab, line = par("mgp")[1], padj = 0.5, 
              cex = par("cex") * par("cex.lab"))
      }
      if (density.info == "density") {
        dens <- density(x, adjust = densadj, na.rm = TRUE, 
                        from = min.scale, to = max.scale)
        omit <- dens$x < min(breaks) | dens$x > max(breaks)
        dens$x <- dens$x[!omit]
        dens$y <- dens$y[!omit]
        dens$x <- scale01(dens$x, min.raw, max.raw)
        lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
              lwd = 1)
        if (is.null(key.ytickfun)) {
          yargs <- list(at = pretty(dens$y)/max(dens$y) * 
                          0.95, labels = pretty(dens$y))
        }
        else {
          yargs <- key.ytickfun()
        }
        yargs$side <- 2
        do.call(axis, yargs)
        if (is.null(key.title)) 
          key.title <- "Color Key\nand Density Plot"
        if (!is.na(key.title)) 
          title(key.title)
        par(cex = 0.5)
        if (is.null(key.ylab)) 
          key.ylab <- "Density"
        if (!is.na(key.ylab)) 
          mtext(side = 2, key.ylab, line = par("mgp")[1], 
                padj = 0.5, cex = par("cex") * par("cex.lab"))
      }
      else if (density.info == "histogram") {
        h <- hist(x, plot = FALSE, breaks = breaks)
        hx <- scale01(breaks, min.raw, max.raw)
        hy <- c(h$counts, h$counts[length(h$counts)])
        lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
              col = denscol)
        if (is.null(key.ytickfun)) {
          yargs <- list(at = pretty(hy)/max(hy) * 0.95, 
                        labels = pretty(hy))
        }
        else {
          yargs <- key.ytickfun()
        }
        yargs$side <- 2
        do.call(axis, yargs)
        if (is.null(key.title)) 
          key.title <- "Color Key\nand Histogram"
        if (!is.na(key.title)) 
          title(key.title)
        par(cex = 0.5)
        if (is.null(key.ylab)) 
          key.ylab <- "Count"
        if (!is.na(key.ylab)) 
          mtext(side = 2, key.ylab, line = par("mgp")[1], 
                padj = 0.5, cex = par("cex") * par("cex.lab"))
      }
      else {
        if (is.null(key.title)) 
          key.title <- "Color Key"
        if (!is.na(key.title)) 
          title(key.title)
      }
      if (trace %in% c("both", "column")) {
        vline.vals <- scale01(vline, min.raw, max.raw)
        if (!is.null(vline)) {
          abline(v = vline.vals, col = linecol, lty = 2)
        }
      }
      if (trace %in% c("both", "row")) {
        hline.vals <- scale01(hline, min.raw, max.raw)
        if (!is.null(hline)) {
          abline(v = hline.vals, col = linecol, lty = 2)
        }
      }
    }
    else {
      par(mar = c(0, 0, 0, 0))
      plot.new()
    }
    retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)], 
                                    high = retval$breaks[-1], color = retval$col)
    retval$layout <- list(lmat = lmat, lhei = lhei, lwid = lwid)
    if (!is.null(extrafun)) 
      extrafun()
    invisible(retval)
  }

