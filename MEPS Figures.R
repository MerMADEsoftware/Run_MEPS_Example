#R Script for recreating the figures shown in MEPS paper: 
#"Spatial Spatial patterns of within-stock connectivity provide novel insights for fisheries management

#libraries to be installed and read in
library(igraph)
library(rgl)
library(raster)
library(dplyr)
library(sp)
library(viridis)
library(gplots)
library(gtools)

#set your working directory
#setwd()

#read in the functions we need (included in the github repository)
source("Analyse_MerMADE_results.R")

#####################################################################
#Recreating Figure 2A: seascape with patches and movement tracks
res=1500 #res : resolution of the seascape, this is the horizontal resolution (in metres) in the x-y direction
min_x=428430 #min_x : minimum x coordinate (in metres)
min_y=6123500 #min_y : minimum y coordinate (in metres)
nr=245 #nr : number of rows in your file
nc=198 #nc : number of columns in your file
max_x=min_x+nc*res
max_y=min_y+nr*res
patches<- as.matrix(read.table("area4_patches.txt", header=F)) #this file is included on the github repository, make sure you've downloaded it to your working directory
patch_raster<- raster(patches, xmn=min_x, xmx=max_x, ymn=min_y, ymx=max_y )
landscape<- as.matrix(read.table("area4_habs.txt", header=F))

#standardise the minimum values. -9999 is a placeholder value in MerMADE but interferes with plotting
for(r in 1:nrow(landscape)){
  for(c in 1:ncol(landscape)){
    if(is.na(landscape[r,c]) || landscape[r,c]==-9999){
      landscape[r,c]=0
    }
  }
}
land_raster<- raster(landscape, xmn=min_x, xmx=max_x, ymn=min_y, ymx=max_y )
cols<- c("white", "darkgoldenrod2", "gray89")
plot(land_raster, axes=F, box=F, asp=1, yaxs="i", xaxs="i", col=cols, legend=F)


##add tracks
track_file<- "Test1/Outputs/Indiv_Mov_rep0.txt"
no_col<- max(count.fields(track_file, sep="\t"))
indiv_mov<- read.table(track_file, skip=1, sep="\t", fill=T, col.names=1:no_col) #this step can take a minute because it's a large file

rep_data<- indiv_mov[indiv_mov[,1]==0,]  
year_data<- rep_data[rep_data[,2]==0,]
nindivs=(unique(year_data[,4]))
colours<- rainbow(length(nindivs))

for(i in 1:length(nindivs)){
  #to plot 10% of tracks
  to_plot=runif(1)
  if(to_plot<0.1){
    indiv_data<- year_data[year_data[,4]==nindivs[i],]
    xy<- cbind(indiv_data[1,-(1:6)][!is.na(indiv_data[1,-(1:6)])],indiv_data[2,-(1:6)][!is.na(indiv_data[2,-(1:6)])])
    pts <- sp::SpatialPoints(xy)
    sp_line <- as(pts,"SpatialLines")
    plot(sp_line, add=T)
  }
}

#######################################################
#Recreating Figure 3:graphic network representation

npatches=43
out_conn=T
in_conn=NULL #because we want it to read in a connectivity matrix we already have
out_outdeg=F #true because we want it to output out-degree
out_indeg=F
nreps=20 #reps: how many replicates are in the file?
simnum="Test1"
peryear=T
inc_local=F #this will disregard individuals that did not disperse due to emigration probability

connmat<- conn_mat_reps(npatches, out_conn, in_conn, out_outdeg, out_indeg, nreps, simnum, peryear, inc_local)
#write.table(format(connmat, digits=3,scientific=F), file="connmat.txt",col.names = T, row.names = T, sep="\t", quote=F)

#read in the vertex coordinates so they are placed as we did in the paper
vert_coords<- as.matrix(read.table("network_vertex_coords.txt"))

###trying to reduce the number of vertices to get a clearer picture
freq_conn<- connmat

for(r in 1:nrow(freq_conn)){
  for(c in 1:ncol(freq_conn)){
    if(connmat[r,c]<0.125){
      freq_conn[r,c]=0
    }
    if(r==c){
      freq_conn[r,c]=0
    }
  }
}
g1<- graph_from_adjacency_matrix(freq_conn, mode=c("directed"), diag=F, weighted=T)
#jpeg("Figure 3.jpeg", width = 8, height = 9, units = 'in', res = 300)
plot(g1, layout=vert_coords, edge.width=(E(g1)$weight/20), 
     vertex.color="white", edge.arrow.size=0.25, vertex.size=10, ylim=c(-.6,.65))
#dev.off()

###################################################
#Figures 4 and 5 were produce in GIS but here is code to analyse the metrics and an approximation of the figures

#we want it to output the out-degree centrality (how many other patches is each patch supplying)
npatches=43
out_conn=F
in_conn=NULL
out_outdeg=T #true because we want it to output out-degree
out_indeg=F
peryear=T
nreps=20
inc_local=F

out_degree<- as.matrix( conn_mat_reps(npatches, out_conn, in_conn, out_outdeg, out_indeg, nreps, simnum, peryear))

#add patch sizes for reference
p_sizes<- read.table("area4_sizes.txt", header=T)
out_degree<- cbind(out_degree, p_sizes[,2])
colnames(out_degree)<- c("Patch", "Origin", "Destination", "Patch Size")
#write.table(out_degree, "centralities.txt", col.names=T, sep="\t", quote=F)

#look at top out patches
top_out<- out_degree[order(out_degree[,2], decreasing=T),c(1,2,4)]
#look at top in patches
top_in<- out_degree[order(out_degree[,3], decreasing=T),c(1,3,4)]


############################
#what about self-recruitment?
patch_self<- matrix(0,nrow=npatches, ncol=4)
colnames(patch_self)<- c("Patch", "Total_recruit", "Total_self", "Prop_self")
patch_self[,1]<- seq(0,(npatches-1),1)
#connmat<- read.table("connectivity_matrix.txt", header=T, row.names=T, sep="\t")
for(p in 1:npatches){
  patch_self[p,2]=sum(connmat[,p]) #take the total of each column to calculate all incoming recruits to that patch
  patch_self[p,3]=connmat[p,p] #take only the diagonal to find the self-recruits
}
patch_self[,4]=patch_self[,3]/patch_self[,2] #take the proportion of self-recruits vs total recruits

################################
#plot your connectivity measures on a 2D seascape
res=1500 #res : resolution of the seascape, this is the horizontal resolution (in metres) in the x-y direction
min_x=428430 #min_x : minimum x coordinate (in metres)
min_y=6123500 #min_y : minimum y coordinate (in metres)
nr=245 #nr : number of rows in your file
nc=198 #nc : number of columns in your file
max_x=min_x+nc*res
max_y=min_y+nr*res
patches<- as.matrix(read.table("area4_patches.txt", header=F))
patch_raster<- raster(patches, xmn=min_x, xmx=max_x, ymn=min_y, ymx=max_y )
outs<- matrix(NA, nrow=nrow(patches), ncol=ncol(patches))
ins<- matrix(NA, nrow=nrow(patches), ncol=ncol(patches))
selfs<- matrix(NA, nrow=nrow(patches), ncol=ncol(patches))

for(r in 1:nrow(patches)){
  for(c in 1:ncol(patches)){
    if(!is.na(patches[r,c]) && patches[r,c]!=-9999){
      pa<- patches[r,c]
      # ins[r,c]=in_degree[pa+1,2]
      # outs[r,c]=out_degree[pa+1,2]      
      ins[r,c]=out_degree[pa+1,3]
      outs[r,c]=out_degree[pa+1,2]
      selfs[r,c]=patch_self[pa+1,4]
    }
  }
}
ins_raster<- raster(ins, xmn=min_x, xmx=max_x, ymn=min_y, ymx=max_y )
outs_raster<- raster(outs, xmn=min_x, xmx=max_x, ymn=min_y, ymx=max_y )
selfs_raster<- raster(selfs, xmn=min_x, xmx=max_x, ymn=min_y, ymx=max_y )

landscape<- as.matrix(read.table("area4_habs.txt", header=F))

#standardise the minimum values. -9999 is a placeholder value in MerMADE but interferes with plotting
for(r in 1:nrow(landscape)){
  for(c in 1:ncol(landscape)){
    if(is.na(landscape[r,c]) || landscape[r,c]==-9999){
      landscape[r,c]=0
    }
    if(landscape[r,c]!=2){landscape[r,c]=0}
  }
}
land_raster<- raster(landscape, xmn=min_x, xmx=max_x, ymn=min_y, ymx=max_y )

#plot your destination centrality
plot(land_raster, axes=F, box=F, asp=1, yaxs="i", xaxs="i", col=c("white", "gray89"), main="Destination centrality", legend=F)
plot(ins_raster, col=rev(heat.colors(10)), asp=1, add=T)

#plot your origin centrality
plot(land_raster, axes=F, box=F, asp=1, yaxs="i", xaxs="i", col=c("white", "gray89"), main="Origin centrality", legend=F)
plot(outs_raster, col=rev(heat.colors(10)), asp=1, add=T)

#plot your self-recruitment 
plot(land_raster, axes=F, box=F, asp=1, yaxs="i", xaxs="i", col=c("white", "gray89"), main="Self-recruitment", legend=F)
plot(selfs_raster, col=rev(heat.colors(10)), asp=1, add=T)

#############################
#Recreating Figure 6: overall population sizes with single vs repeated depletion
#jpeg("Figure 6.jpeg", width = 8, height = 4, units = 'in', res = 300) #if you wanted to save a jpeg
tests<- c(1,2,3,4,5,6)
time=50
reps=20
dir="test"
layout(matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE))
par( oma= c(1,3,1,0), mai= c(.2,.2,.2,.1))
plot(0, type="l", xlim=c(0,time), ylim=c(1500, 4000), xlab="Time", ylab="Population size", xaxt="n")
mtext("A", side=3, line=1,adj=0, cex=1.2)
colours<- c("black", "blue", "red", "green", "orange", "darkorchid3", "yellow")
for(i in 1:length(tests)){
  p_ext_pop<- read.table(paste(dir,tests[i], "/Outputs/Population.txt", sep=""), header=T, as.is=T)
  tot_series<-c()
  reps=max(p_ext_pop$Rep)+1
  for(t in 0:(time-1)){
    p_t<- p_ext_pop[p_ext_pop$Year==t,]
    p_tot=0;
    for(r in 0:(reps-1)){
      p_r<- p_t[p_t$Rep==r,]
      p_tot<- p_tot+ sum(p_r$NInd)
    }
    p_mean=p_tot/reps
    tot_series<- c(tot_series, p_mean)
  }
  lines(x=seq(0,time-1, 1), y=tot_series, col=colours[i], lwd=1.5)
}
abline(v=20, lty="dashed")
par(mai=c(.2,.1,.2,.1))
plot(0, type="l", xlim=c(18,time-1), ylim=c(1500, 3000), xlab=" ", ylab=" ", yaxt="n", xaxt="n")
mtext("B", side=3, line=1,adj=0, cex=1.2)

for(i in 1:length(tests)){
  p_ext_pop<- read.table(paste(dir,tests[i], "/Outputs/Population.txt", sep=""), header=T, as.is=T)
  tot_series<-c()
  reps=max(p_ext_pop$Rep)+1
  for(t in 0:(time-1)){
    p_t<- p_ext_pop[p_ext_pop$Year==t,]
    p_tot=0;
    for(r in 0:(reps-1)){
      p_r<- p_t[p_t$Rep==r,]
      p_tot<- p_tot+ sum(p_r$NInd)
    }
    p_mean=p_tot/reps
    tot_series<- c(tot_series, p_mean)
  }
  lines(x=seq(18,time-1, 1), y=tot_series[19:length(tot_series)], col=colours[i], lwd=1.5)
}
abline(v=20, lty="dashed")


tests<- c(1,7,8,9,10,11)
par(mai=c(.2,.2,.2,.1))
plot(0, type="l", xlim=c(0,time), ylim=c(1500, 4000), xlab="Time", ylab="Population size")
mtext("                                       Population Size", side=2, line=3, cex=1.1)
mtext("                                                                    Time (years)", side=1, line=3, cex=1.1)
mtext("C", side=3, line=1,adj=0, cex=1.2)
colours<- c("black", "blue", "red", "green", "orange", "darkorchid3", "yellow")
for(i in 1:length(tests)){
  p_ext_pop<- read.table(paste(dir,tests[i], "/Outputs/Population.txt", sep=""), header=T, as.is=T)
  tot_series<-c()
  reps=max(p_ext_pop$Rep)+1
  for(t in 0:(time-1)){
    p_t<- p_ext_pop[p_ext_pop$Year==t,]
    p_tot=0;
    for(r in 0:(reps-1)){
      p_r<- p_t[p_t$Rep==r,]
      p_tot<- p_tot+ sum(p_r$NInd)
    }
    p_mean=p_tot/reps
    tot_series<- c(tot_series, p_mean)
  }
  lines(x=seq(0,time-1, 1), y=tot_series, col=colours[i], lwd=1.5)
}
#abline(v=20, lty="dashed")
ext_years<- seq(20, 49,2)
for(ab in 1:length(ext_years)){
  abline(v=ext_years[ab], lty="dashed")
}
par(mai=c(.2, .1, .2, .1))
plot(0, type="l", xlim=c(18,time-1), ylim=c(1500, 3000), xlab=" ", ylab=" ", yaxt="n")
mtext("D", side=3, line=1,adj=0, cex=1.2)

for(i in 1:length(tests)){
  p_ext_pop<- read.table(paste(dir,tests[i], "/Outputs/Population.txt", sep=""), header=T, as.is=T)
  tot_series<-c()
  reps=max(p_ext_pop$Rep)+1
  
  for(t in 0:(time-1)){
    p_t<- p_ext_pop[p_ext_pop$Year==t,]
    p_tot=0;
    for(r in 0:(reps-1)){
      p_r<- p_t[p_t$Rep==r,]
      p_tot<- p_tot+ sum(p_r$NInd)
    }
    p_mean=p_tot/reps
    tot_series<- c(tot_series, p_mean)
  }
  lines(x=seq(18,time-1, 1), y=tot_series[19:length(tot_series)], col=colours[i], lwd=1.5)
}
#abline(v=20, lty="dashed")
ext_years<- seq(20, 49,2)
for(ab in 1:length(ext_years)){
  abline(v=ext_years[ab], lty="dashed")
}

par( mai= c(0,0,.5,0))
plot(0, type="l", xlim=c(0,100), ylim=c(0,10), axes=F, xlab=" ", ylab=" ")
legend(x=0, y=10, legend=c("Control", "Patch 26", "Patch 39", "Patch 3", "Patch 14", "Patch 32"), col=colours, 
       title= "Depletion at:", lty="solid", lwd=2, ncol=6, cex=1.1)

#dev.off()
##################################################
#Figure 7: patch level effects of single vs repeated depletion
patches<-c(26,39,3,14,32) #patches targeted
tests_r<- c(7,8,9,10,11) #repeated depletion tests
tests_s<- c(2,3,4,5,6) #single depletion tests
tests_c<-c(1) #control
reps=20
time=50
#jpeg("Figure 7.jpeg", width = 8, height = 4, units = 'in', res = 300) #if you want to save a jpeg
types<- c("solid", "dashed", "dotted")
layout(matrix(c(1,2,3,4,5,6), nrow = 2, ncol = 3, byrow = TRUE))
par( oma= c(1,2,0,0), mai= c(.5,.4,.3,.3))
#par(mfrow=c(3,3))
for(p in 1:length(patches)){
  if(p==1){
    par(mai=c(0.2,0.4,0.2,0))
    plot(0, type="l", xlim=c(10,time), ylim=c(0, 600), xlab=" ", ylab=" ", xaxt="n",
         main=paste("Patch", patches[p], sep=""))
  }
  else if(p==2){
    par(mai=c(0.2,0.2,0.2,0))
    plot(0, type="l", xlim=c(10,time), ylim=c(0, 400), xlab=" ", ylab=" ", xaxt="n", yaxt="n", 
         main=paste("Patch", patches[p], sep=""))
  }
  else if(p==3){
    par(mai=c(0.2,0.2,0.2,0.2))
    plot(0, type="l", xlim=c(10,time), ylim=c(0, 600), xlab=" ", ylab=" ", yaxt="n", 
         main=paste("Patch", patches[p], sep=""))
  }
  else if(p==4){
    par(mai=c(0.4,0.4,0.2,0))
    plot(0, type="l", xlim=c(10,time), ylim=c(0, 150), xlab=" ", ylab=" ", 
         main=paste("Patch", patches[p], sep=""))
    mtext("                                            Number of individuals in patch", side=2, line=3, las=0, padj=0, cex=1.2)
    
  }
  else if(p==5){
    par(mai=c(0.4,0.2,0.2,0))
    plot(0, type="l", xlim=c(10,time), ylim=c(0, 80), xlab=" ", ylab=" ", yaxt="n", 
         main=paste("Patch", patches[p], sep=""))
    mtext("Time (years)", side=1, line=3, cex=1.5)
  }
  else{
    plot(0, type="l", xlim=c(10,time), ylim=c(0, 20), xlab="Time (years)", ylab=" ", 
         main=paste("Patch", patches[p], sep=""))
  }
  
  to_do<- c(tests_c, tests_s[p], tests_r[p])
  
  for(i in 1:length(to_do)){
    p_ext_pop<- read.table(paste(dir,to_do[i], "/Outputs/Population.txt", sep=""), header=T, as.is=T)
    tot_series<-c()
    
    for(t in 10:(time-1)){
      p_t<- p_ext_pop[p_ext_pop$Year==t,]
      p_tot<- p_t[p_t$Patch_ID==patches[p],]
      if(nrow(p_tot)==0){
        #print(paste("in year ", t , " there were no indivs in patch " , patches[p] , sep=""))
        p_tot=0
      }
      else{
        p_tot=mean(p_tot$NInd)
      }
      tot_series<- c(tot_series, p_tot)
    }
    lines(x=seq(10,time-1, 1), y=tot_series,lty=types[i])
  }
  
}
par( mai= c(0,0,0,0))

plot(0, type="l", xlim=c(0,50), ylim=c(0,10), axes=F, xlab=" ", ylab=" ")
legend(x=3, y=8, legend=c("Control", "Single Depletion", "Repeated Depletion"), 
       lty=types, ncol=1, cex=1.2)

#############################################################
###Recreating Appendix figure A1

#jpeg("Figure A1a.jpeg", width = 8, height = 8, units = 'in', res = 300)
lmat = rbind(c(0,4,0),c(2,1,3)) #4=key, 2=row dendrogram, 1= heatmap, 3=col dendrogram
lwid=c(.5,4,.5)
lhei=c(1,4)
breaks1<- c(0, 1, 5, 10, 25, 50, 75, 100, 500, 1000)

#remove the self-recruitment diagonal from connectivity matrix for clarity
connmat_nodiag<- connmat
for(r in 1:nrow(connmat)){
  for(c in 1:ncol(connmat)){
    if(r==c){connmat_nodiag[r,c]=0}
  }
}

#creating Figure A1A: per year connectivity
new_heatmap2(connmat_nodiag,scale="none", dendrogram="none", lmat=lmat, lwid=lwid,
             lhei=lhei, col=c("white", rev(gray.colors(8))), breaks=breaks1, key.ylab=NA, key.title="", tracecol=NA,
             Colv=NA, Rowv=NA , key=T)
#dev.off()
#creating FigureA1B: 50-year cumulative connectivity
#jpeg("Figure A1b.jpeg", width = 8, height = 8, units = 'in', res = 300)

#calculate the connectivity matrix of total movement across the 50-year sim, not a per year average
npatches=43
out_conn=T
in_conn=NULL 
out_outdeg=F 
out_indeg=F
nreps=20 
simnum="Test1"
peryear=F #because we want to produce the heatmap over the whole simulation time period
inc_local=F

connmat_total<- conn_mat_reps(npatches, out_conn, in_conn, out_outdeg, out_indeg, nreps, simnum, peryear, inc_local)
connmat_nodiag<- connmat_total
for(r in 1:nrow(connmat_total)){
  for(c in 1:ncol(connmat_total)){
    if(r==c){connmat_nodiag[r,c]=0}
  }
}

new_heatmap2(connmat_nodiag,scale="none", dendrogram="none", lmat=lmat, lwid=lwid,
             lhei=lhei, col=c("white", rev(gray.colors(8))), breaks=breaks1, key.ylab=NA, key.title="", tracecol=NA,
             Colv=NA, Rowv=NA , key=T)


