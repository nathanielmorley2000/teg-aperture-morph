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

# Create a vector to store indices of the selected points
selected_points <- data.frame(c(i, x, y))

for (i in 1:length(silhouettes)) {
  plot(silhouettes[i])
  cat("Click on the homologous point for shape", i, "\n")
  pt <- locator(1)
  coords <- silhouettes[i]$coo
  distances <- sqrt((coords[,1] - pt$x)^2 + (coords[,2] - pt$y)^2)
  selected_points[x, i] <- pt$x #which.min(distances)
}


# Loop through each shape in the Momocs Out object
for (i in 1:length(silhouettes)) {
  plot(silhouettes[i], main = paste("Click homologous point for shape", i))
  
  # Manually click the homologous point
  pt <- locator(1)
  
  # Extract the (x, y) coordinates of the current silhouette
  coords <- silhouettes[i]$coo
  
  # Calculate distances from the clicked point to all outline points
  distances <- sqrt((coords[, 1] - pt$x)^2 + (coords[, 2] - pt$y)^2)
  
  # Store the index of the closest point
  selected_points[i] <- which.min(distances)
}


# minimize size differences
dorsal.smooth <- coo_smooth(silhouettes, 100)
pile(dorsal.smooth)

dorsal.scale<-coo_scale(dorsal.smooth)
pile(dorsal.scale)

dorsal.center<-coo_center(dorsal.scale)
pile(dorsal.center)

dorsal.align <- coo_align(dorsal.center)
pile(dorsal.align)

ldk <- coo_ldk(dorsal.align, 1)

i = 1
while (i < length(dorsal.align) + 1) {
  ldk[i] <- coo_ldk(dorsal.align[i], 1)
  i = i + 1
}

i = 1
while (i < length(dorsal.align) + 1) {
  dorsal.slide[i] <- coo_slide(dorsal.align[i], id = 1, ldk = ldk[i])
  i = i + 1
}



dorsal.slide <- coo_slide(dorsal.align[1], id = 1, ldk = ldk)
pile(dorsal.slide)



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
