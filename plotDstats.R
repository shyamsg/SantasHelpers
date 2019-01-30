###################
# Function to plot D stats
###################

plotDstats_H1Var <- function(allstats, h2, h3, h2name, h3name, outgroupname, layoutMatrix, xlimits) {
    ################ Get subset of stats 
  stats = subset(allstats, (H2==h2 | H1==h2) & H3==h3)
  flipped = which(stats$H1 == h2)
  stats$Dstat[flipped] = -stats$Dstat[flipped]
  stats$jackEst[flipped] = -stats$jackEst[flipped]
  stats$Z[flipped] = -stats$Z[flipped]
  temp = stats$H2[flipped]
  stats$H2[flipped] = h2
  stats$H1[flipped] = temp
  
  # # start plot
  # lo<-matrix(0,40,80)
  # lo[1:2,5:78]<-1
  # lo[2:8,27:44]<-2
  # lo[2:8,45:62]<-3
  # lo[2:8,63:80]<-4
  layout(layoutMatrix)
  
  par(mar=c(0,0,0,0))
  plot(1, col="white", xaxt="n", yaxt="n", bty="n")
  par(mar=c(3.6, 2, 2, 1))
  m<-matrix(c(c(1,1,1,1), c(1,1,1,2), c(1,1,3,3), c(1,4,4,4)), 4,4,byrow=TRUE)
  labels = c(h2name, "H1", h3name, outgroupname)
  hc <- hclust(dist(m))
  hc$labels<-c(labels)
  hc$order<-rev(hc$order)
  plot(hc, hang=-1, yaxt="n", main="D < 0", ylab="", xlab="", sub="")
  segments(1, .5, 3, .5, col="orchid3")
  
  plot(hc, hang=-1, yaxt="n", main="D = 0", ylab="", xlab="", sub="")
  
  plot(hc, hang=-1, yaxt="n", main="D > 0", ylab="", xlab="", sub="")
  segments(2, .5, 3, .5, col="orchid3")
  
  l2Int = qnorm(0.0005, lower.tail=F)
  limits1 = aes(xmax = Dstat + SE, xmin= Dstat - SE) 
  limits2 = aes(xmax = Dstat + l2Int*SE, xmin= Dstat - l2Int*SE)
  p<-ggplot(stats, aes(y=reorder(H1, Dstat), x=Dstat), order=stats$Dstat)+
    geom_point() +
    geom_errorbarh(limits1, height=0.5, lwd=0.5)+
    geom_errorbarh(limits2, height=0.3, lwd=0.5)+
    geom_vline(aes(xintercept=0), alpha=0.5, colour="orange", linetype=1)+
    ylab("H2")+
    xlab("D")+
    xlim(xlimits) +
    theme(panel.background=element_rect(fill='white', colour="lightsteelblue"), panel.grid.major.y=element_line(colour="lightsteelblue", size=0.2, linetype="dotted"))
  
  our.first.vp <- viewport(x = 0.5, y = 0, height = .84, width = 1, just = c("centre", "bottom"))
  
  pushViewport(our.first.vp)
  
  print(p, newpage=FALSE)
  #message(length(unique(Dstat$H2)))
  upViewport(1)
}

plotDstats_H3Var <- function(allstats, h2, h1, h2name, h1name, outgroupname, layoutMatrix, xlimits) {
  ################ Get subset of stats 
  stats = subset(allstats, (H2==h2 & H1==h1) | (H2==h1 & H1==h2) )
  flipped = which(stats$H1 == h2)
  stats$Dstat[flipped] = -stats$Dstat[flipped]
  stats$jackEst[flipped] = -stats$jackEst[flipped]
  stats$Z[flipped] = -stats$Z[flipped]
  temp = stats$H2[flipped]
  stats$H2[flipped] = h2
  stats$H1[flipped] = temp
  
  # # start plot
  # lo<-matrix(0,40,80)
  # lo[1:2,5:78]<-1
  # lo[2:8,27:44]<-2
  # lo[2:8,45:62]<-3
  # lo[2:8,63:80]<-4
  layout(layoutMatrix)
  
  par(mar=c(0,0,0,0))
  plot(1, col="white", xaxt="n", yaxt="n", bty="n")
  par(mar=c(3.6, 2, 2, 1))
  m<-matrix(c(c(1,1,1,1), c(1,1,1,2), c(1,1,3,3), c(1,4,4,4)), 4,4,byrow=TRUE)
  labels = c(h2name, h1name, "H3", outgroupname)
  hc <- hclust(dist(m))
  hc$labels<-c(labels)
  hc$order<-rev(hc$order)
  plot(hc, hang=-1, yaxt="n", main="D < 0", ylab="", xlab="", sub="")
  segments(1, .5, 3, .5, col="orchid3")
  
  plot(hc, hang=-1, yaxt="n", main="D = 0", ylab="", xlab="", sub="")
  
  plot(hc, hang=-1, yaxt="n", main="D > 0", ylab="", xlab="", sub="")
  segments(2, .5, 3, .5, col="orchid3")
  
  l2Int = qnorm(0.0005, lower.tail=F)
  limits1 = aes(xmax = Dstat + SE, xmin= Dstat - SE) 
  limits2 = aes(xmax = Dstat + l2Int*SE, xmin= Dstat - l2Int*SE)
  p<-ggplot(stats, aes(y=reorder(H3, Dstat), x=Dstat), order=stats$Dstat)+
    geom_point() +
    geom_errorbarh(limits1, height=0.5, lwd=0.5)+
    geom_errorbarh(limits2, height=0.3, lwd=0.5)+
    geom_vline(aes(xintercept=0), alpha=0.5, colour="orange", linetype=1)+
    ylab("H2")+
    xlab("D")+
    xlim(xlimits) +
    theme(panel.background=element_rect(fill='white', colour="lightsteelblue"), panel.grid.major.y=element_line(colour="lightsteelblue", size=0.2, linetype="dotted"))
  
  our.first.vp <- viewport(x = 0.5, y = 0, height = .84, width = 1, just = c("centre", "bottom"))
  
  pushViewport(our.first.vp)
  
  print(p, newpage=FALSE)
  #message(length(unique(Dstat$H2)))
  upViewport(1)
}