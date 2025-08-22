# Call packages
library("geomorph")
library("stringr")
#library("ggplot2")

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

# Save GPA results to dataframe
GPA.Results <- geomorph.data.frame(Procrustes, 
                                   location = as.factor(specific.characteristics$Location),
                                   energy = as.factor(specific.characteristics$Energy.Setting),
                                   repairs = as.factor(specific.characteristics$Repair.Scars..0.1.))


# Create PCA plot
PCA <- gm.prcomp(GPA.Results$coords)
PCA.Energy <- plot(PCA, main = "PCA", col = GPA.Results$energy, theme = classic)
shapeHulls(PCA.Energy, groups = GPA.Results$energy,
           group.cols = 1:3, 
           group.lwd = c(3,3,3))
legend("bottomright", c("High", "Low", "Moderate"), 
       col = 1:3, lwd = 2)


pc.plot <- plot(pleth.phylo, phylo = TRUE)
gp <- factor(c(rep(1, 5), rep(2, 4)))
shapeHulls(pc.plot, groups = gp, group.cols = 1:2, 
           group.lwd = rep(2, 2), group.lty = c(2, 1))
legend("topright", c("P. cinereus clade", "P. hubrichti clade"), 
       col = 1:2, lwd = 2, lty = c(2, 1))











plot(PCA, main = "PCA", col = GPA.Results$repairs)


PCA.Plot <- plot(PCA, main = "PCA")
picknplot.shape(PCA.Plot)


summary(PCA)
plot(PCA, main = "PCA")
plot(PCA, main = "PCA", flip = 1) # flip the first axis
plot(PCA, main = "PCA", axis1 = 3, axis2 = 4) # change PCs viewed

