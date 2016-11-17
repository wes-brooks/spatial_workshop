#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MN Lung Cancer Example in Section 9.7
#
#
library('maps', 'maptools', 'RColorBrewer')
lung <- read.table("data/MN/MNlung.txt", header=TRUE, sep="\t")
radon <- read.table("data/MN/MNradon.txt", header=TRUE)
Obs <- apply(cbind(lung[,3], lung[,5]), 1, sum)
Exp <- apply(cbind(lung[,4], lung[,6]), 1, sum)
rad.avg <- rep(0, length(lung$X))
for(i in 1:length(lung$X)) {
	rad.avg[i]<-mean(radon[radon$county==i,2])
}
x <- rad.avg
#
rad.avg[26] <- 0
rad.avg[63] <- 0
x[26] <- NA
x[63] <- NA
#
poismod <- glm(Obs~offset(log(Exp))+x,family="poisson")
quasimod <- glm(Obs~offset(log(Exp))+x,family=quasipoisson(link="log"))
#
# Use inla to fit non-spatial model and ICAR+non-spatial model, both including radon 
# as a loglinear covariate. 
#
source("http://www.math.ntnu.no/inla/givemeINLA.R")
inla.upgrade()
library(INLA)
#
# We give two different ways of obtaining shapefiles
#
# with projection
library(rgdal)
county.shp <- readOGR(dsn=getwd(),"tl_2009_27_county00")
proj4string(county.shp)
par(mfrow=c(1,1))
plot(county.shp) # Just plots boundaries
#
#
county.shp <-readShapePoly("tl_2009_27_county00.shp", 
proj4string=CRS("+proj=longlat +ellps=GRS80 +datum=NAD27"))
plot(county.shp) # Just plots boundaries
#
# get data ready
#
dat.inla <- data.frame(Obs=Obs, Exp=Exp, x=x)
dat.inla$region.struct <- 1:nrow(dat.inla)
dat.inla$region.unstruct <- 1:nrow(dat.inla)
# Poisson model with non-spatial only
formula <- Obs ~ x + f(region.unstruct,model="iid",param=c(1,0.140))
#
poismod.inla0 <- inla(formula, offset=I(log(Exp)), family="poisson", data=dat.inla,
                 control.predictor=list(compute=TRUE))
summary(poismod.inla0)
#
# Prior on the association parameter: 95% interval of the relative risk associated with
# 1-unit change in radon is [0.1,10].
#
upper <- 10
W <- (log(upper)/1.96)^2
#
# Poisson model with ICAR spatial random effects and non-spatial random effects, and
# including radon as a loglinear covariate
#
formula <- Obs ~ x + f(region.struct,model="besag",graph.file="MN_county.graph", 
           param=c(0.21,0.140)) + f(region.unstruct,model="iid",param=c(1,0.140))
# 
# Fit Poisson model:
#
poismod.inla <- inla(formula, offset=I(log(Exp)), family="poisson", data=dat.inla,
               control.predictor=list(compute=TRUE),control.fixed=list(mean=c(0),prec=c(1/W)))
summary(poismod.inla)
#
# spatial random effects
#
u <- poismod.inla$summary.random$region.struct$mean
sd.u <- sd(u)
sd.u
#
# non-spatial random effects
#
v <- poismod.inla$summary.random$region.unstruct$mean
hyper.param <- inla.hyperpar(poismod.inla)
summary(hyper.param)
sigma.region.unstruct <- sqrt(1/48.29) # This is the median of sigma_epi
#
# ratio of spatial to total variation
#
sd.u^2/(sd.u^2 + sigma.region.unstruct^2) 
#
plotvar <- u
summary(plotvar)
# symbol plot -- equal-interval class intervals
nclr <- 8
plotclr <- brewer.pal(nclr,"Greys")
brks <- quantile(plotvar,probs=seq(0,1,1/(nclr)))
colornum <- findInterval(plotvar, brks, all.inside=T)
colcode <- plotclr[colornum]
#
# Fig 9.5(b)
#
pdf("MN-Spatial.pdf")
plot(county.shp, col=colcode, lwd=0.5, axes=F)
legend(-92, 46.5, legend=leglabs(round(brks,digits=2)),fill=plotclr,cex=0.8,bty="n")
dev.off()
#
# Now non-spatial
#
plotvar <- v
summary(plotvar)
# symbol plot -- equal-interval class intervals
nclr <- 8
plotclr <- brewer.pal(nclr,"Greys")
brks <- quantile(plotvar,probs=seq(0,1,1/(nclr)))
colornum <- findInterval(plotvar, brks, all.inside=T)
colcode <- plotclr[colornum]
#
# Fig 9.5(a)
#
pdf("MN-NonSpatial.pdf")
plot(county.shp, col=colcode, lwd=0.5, axes=F)
legend(-92, 46.5, legend=leglabs(round(brks,digits=2)),fill=plotclr,cex=0.8,bty="n")
dev.off()



