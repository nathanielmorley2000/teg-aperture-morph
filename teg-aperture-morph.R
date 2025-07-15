# load libraries
library("Momocs")
library(dplyr)
library(ggplot2)
library(borealis)
library(dispRity)

# load characteristics list
characteristics <- read.csv("Data/TegulaCharacteristicsList.csv")

# import jpgs 
silhouettes.jpgs<-list.files("Data/Silhouettes/", full.names=T)
returns <- import_jpg(silhouettes.jpgs, auto.notcentered = T)
silhouettes<-Out(returns)
pile(silhouettes)

# minimize size differences
dorsal.smooth <- coo_smooth(silhouettes, 100)
pile(dorsal.smooth)

dorsal.scale<-coo_scale(dorsal.smooth)
pile(dorsal.scale)

dorsal.center<-coo_center(dorsal.scale)
pile(dorsal.center)

dorsal.align <- coo_align(dorsal.center)
pile(dorsal.align)

dorsal.slide<-(coo_slidedirection(dorsal.align,direction="up"))

pile(dorsal.slide)

calibrate_harmonicpower_efourier(dorsal.slide)
efou_dorsal<-efourier(dorsal.slide,norm=F,nb.h=17)

pca_dorsal <- PCA(efou_dorsal,  center=T)
summary(pca_dorsal)
PC_contrib_dorsal <- PCcontrib(pca_dorsal, nax=c(1:4), sd.r = c(-2, -1, -0.5, 0, 0.5, 1, 2), gap = 1)
PC_contrib_dorsal$gg + geom_polygon(fill="slategrey", col="black") + ggtitle("PC LOADINGS") + coord_flip()

plot.new()
plot_PCA(pca_dorsal, labelpoints = T, zoom = 1)  




# begin eliptical fourier analysis
PrasTeg50 <- import_jpg1("Data/PrasTeg50.jpg")
PrasTeg50og <- import_jpg("Data/PrasTeg50.jpg", auto.notcentered=TRUE, fun.notcentered=NULL,
                          threshold=0.5)
coo<-Out(PrasTeg50)
coo$fac<-c("PrasTeg50")
FF <-efourier(coo_scale(coo_center(coo)))
