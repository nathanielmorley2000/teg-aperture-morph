# load libraries
library("StereoMorph")

digitizeImage(image.file='Data/PrasiolaHighEnergy/', shapes.file='Data/Shapes/')














library("Momocs")
#library(dplyr)
#library(ggplot2)
#library(borealis)
#library(dispRity)

# load characteristics list
characteristics <- read.csv("Data/TegulaCharacteristicsList.csv")

# import jpgs 
silhouettes.jpgs<-list.files("Data/PrasHighEnergy_Open/", full.names=T)
returns <- import_Conte(silhouettes.jpgs)

returns <- import_jpg("Data/PrasHighEnergy_Aligned/PrasTeg50.jpg")
silhouettes<-Coo(returns)
pile(silhouettes)

# smooth, scale, and center silhouettes
silhouettes.smooth <- coo_smooth(silhouettes, 1000)
silhouettes.scale<-coo_scale(silhouettes.smooth)
silhouettes.center<-coo_center(silhouettes.scale)

# define landmark
silhouettes.ldk <- def_ldk(silhouettes.center, 1)
stack(silhouettes.ldk)

# set landmark as starting point
slide.silhouettes <- coo_slide(silhouettes.ldk,ldk = 1)
stack(slide.silhouettes)











dorsal.align <- coo_align(dorsal.center)
pile(dorsal.align)

ldk.version <- def_ldk(dorsal.align, 1)
stack(ldk.version)



i = 1
while (i < length(dorsal.align) + 1) {
  dorsal.align$ldk[i] <- coo_ldk(dorsal.align[i], 1)
  i = i + 1
}

stack(dorsal.align)

prebot <- silhouettes %>% coo_center %>% coo_scale %>%
  coo_align %>% coo_slidedirection("right")
prebot %>% stack # some dephasing remains
prebot %>% coo_slidedirection("right") %>% coo_untiltx() %>% stack # much better



landmarked <- silhouettes %>% coo_center %>% coo_untiltx(ldk=1)







dorsal.smooth <- coo_smooth(silhouettes, 1000)
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

prebot <- dorsal.align %>% coo_slidedirection("right")
prebot %>% stack # some dephasing remains
prebot %>% coo_slidedirection("right") %>% coo_untiltx() %>% stack # much better


dorsal.align %>% coo_untiltx(ldk=1) %>% stack



align.out <- Out(dorsal.align)

dorsal.slide<-coo_slide(align.out, ldk = 1)
pile(dorsal.slide)

align.out <- Out(dorsal.slide)
coo.until <- coo_untiltx(dorsal.slide)
pile(coo.until)

coo_until <- coo_untiltx(dorsal.slide, ldk=1)
pile(coo_until)


