# load package
library("StereoMorph")

# launch StereoMorph application
digitizeImages(image.file='Data/PrasiolaHighEnergy/', shapes.file='Data/Shapes/', curves.ref = "curves.txt")
