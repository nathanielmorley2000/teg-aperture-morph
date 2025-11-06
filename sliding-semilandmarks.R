# Call packages
library("geomorph")
library("stringr")
library("ggplot2")

############################## Generalized Procrustes Analysis ##############################

# Upload TPS file with semi-landmark data
landmarks <- readland.tps(file = "Data/TegulaOutlines.TPS", specID = "imageID", negNA = FALSE,
                          readcurves = TRUE, warnmsg = TRUE)

# Change specimen names to exclude directory information
two.d.landmarks <- two.d.array(landmarks)
rownames(two.d.landmarks) <- stringr::str_sub(rownames(two.d.landmarks), start = 78L)
named.landmarks <- arrayspecs(two.d.landmarks, 90, 2)

# Upload specimen characteristics data
characteristics <- read.csv("Data/TegulaCharacteristicsList.csv", row.names = 1)

# Subset characteristics to only include specimens that were used in analysis
specific.characteristics <- characteristics[row.names(characteristics) %in% rownames(two.d.landmarks),]

# Upload sliders file
curves <- as.matrix(read.csv("Data/Sliders.csv", header=T))

# Run Generalized Procrustes Analysis by using bending energy for sliding
Procrustes <- gpagen(named.landmarks,
                     curves = curves, 
                     ProcD = FALSE, print.progress = FALSE)
plot(Procrustes)
summary(Procrustes)

# Save Procrustes coordinates results as CSV file
two.d.Procrustes <- two.d.array(Procrustes$coords)
write.csv(two.d.Procrustes, "Results/ProcrustesCoords.csv")

# Save GPA results to dataframe
GPA.Results <- geomorph.data.frame(Procrustes, 
                                   location = as.factor(specific.characteristics$Location),
                                   energy = as.factor(specific.characteristics$Energy.Setting),
                                   repairs = as.factor(specific.characteristics$Repair.Scars..0.1.),
                                   height = as.numeric(specific.characteristics$Height..mm.))

#############################################################################################

############################## Ordinations and Statistical Analysis ##############################

# Produce simple PCA (no size correction)
PCA <- gm.prcomp(GPA.Results$coords)
plot(PCA, col = GPA.Results$location)
legend("bottomright", legend = levels(GPA.Results$location), 
       col = 1:3, pch = 16)

# Produce simple allometry PCA
fit <- procD.lm(coords ~ log(height), data = GPA.Results,
                print.progress = FALSE)
plotAllometry(fit, size = GPA.Results$height, logsz = TRUE,
              method = "RegScore", pch = 19)


# Produce group allometries
fit2 <- procD.lm(coords ~ log(height) + location, data = GPA.Results,
                print.progress = FALSE)
fit3 <- procD.lm(coords ~ log(height) * location, data = GPA.Results,
                 print.progress = FALSE)

plot2 <- plotAllometry(fit2, size = GPA.Results$height, logsz = TRUE, method = "RegScore",
              pch = 19, col = as.numeric(GPA.Results$location))
legend("bottomright", legend = levels(GPA.Results$location), 
       col = 1:3, pch = 16)

morph <- morphol.disparity(fit2)
summary(morph)


# find how shape changes along PC axes
ref <- mshape(Procrustes$coords)
plotRefToTarget(PCA$shapes$shapes.comp1$min, ref, method="TPS") #PC1 minimum value
plotRefToTarget(PCA$shapes$shapes.comp1$max, ref, method="TPS") #PC1 maximum value

# Mathers = 1 = black
# Prasiola = 2 = red
# Strawberry = 3 = green

picknplot.shape(plot2)

##################################################################################################