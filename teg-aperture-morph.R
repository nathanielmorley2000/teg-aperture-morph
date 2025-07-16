# load libraries
library("Momocs")
library(dplyr)
library(ggplot2)
#library(borealis)
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

i = 1
while (i < length(dorsal.align) + 1) {
  dorsal.align$ldk[i] <- coo_ldk(dorsal.align[i], 1)
  i = i + 1
}

dorsal.slide<-coo_slide(dorsal.align, ldk = 1)
pile(dorsal.slide)

align.out <- Out(dorsal.slide)
coo.until <- coo_untiltx(dorsal.slide)
pile(coo.until)

coo_until <- coo_untiltx(dorsal.slide, ldk=1)
pile(coo_until)


