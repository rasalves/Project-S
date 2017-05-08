library(spdep)
require(maptools)

Corner_text <- function(text, location="topright"){
legend(location,legend=text, bty ="n", pch=NA) 
}
	  
bh <- readShapePoly(system.file(paste("etc/shapes/","boston_tracts.shp",sep=""),
      package="spdep")[1])	 
stepat = (max(bh$LON)-min(bh$LON))/3

#CONSTRUCTING THE REGION CLUSTERS

var1 = rep(0,length(bh$LON))

var1[bh$LON<(min(bh$LON)+stepat)]=rnorm(sum(bh$LON<(min(bh$LON)+stepat)),0,1)
var1[bh$LON>(min(bh$LON)+stepat+stepat)]=rnorm(sum(bh$LON>(min(bh$LON)+stepat+stepat)),10,1)
var1[var1==0]=rnorm(sum(var1==0),20,1)

varnorm = ((var1-min(var1))/(max(var1)-min(var1)))
par(mar=c(0,0,0,0))
plot(bh, border="white", col=rgb(1-varnorm,1-varnorm,1-varnorm))
Corner_text(text="GENERATED REGION")

nPrune = 2


rndt.bh <- rndtree(nb.w,varnorm)
res_ms_2=mySkater(rndt.bh ,varnorm,nPrune)
par(mar=c(0,0,0,0))
pl=plot(bh, border="white", col=res_ms_2$regions+1,xlab = "RDN TREE - SKT")
Corner_text(text="RDN TREE - SKT")

mst.bh <- kruskal(nb.w,varnorm)
res_ms_1=mySkater(mst.bh,varnorm,nPrune)
par(mar=c(0,0,0,0))
plot(bh, border="white", col=(res_ms_1$regions+1),xlab = "MS TREE - SKT")
Corner_text(text="MS TREE - SKT")

rndt.bh <- rndtree(nb.w,varnorm)
res_ms_2=mySkater_rnd(rndt.bh ,varnorm,nPrune)
par(mar=c(0,0,0,0))
plot(bh, border="white", col=res_ms_2$regions+1,xlab = "RDN TREE - RND SKT")
Corner_text(text="RDN TREE - RND SKT")

mst.bh <- kruskal(nb.w,varnorm)
res_ms_1=mySkater_rnd(mst.bh,varnorm,nPrune)
par(mar=c(0,0,0,0))
plot(bh, border="white", col=res_ms_1$regions+1,xlab = "MS TREE - RDN SKT")
Corner_text(text="MS TREE - RDN SKT")



