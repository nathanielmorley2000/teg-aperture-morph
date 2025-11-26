# Call packages
library("geomorph")
library("stringr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("gridExtra")

#############################################################################################
############################## Generalized Procrustes Analysis ##############################
#############################################################################################

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

##################################################################################################
############################## Ordinations and Statistical Analysis ##############################
##################################################################################################

# Produce simple PCA (no size correction)
PCA <- gm.prcomp(GPA.Results$coords)
summary(PCA)

# First three PCs explain the most variation. Create dataframe with first three PCs
PCA_Results <- data.frame(Location = specific.characteristics$Location,
                          RepairScars = specific.characteristics$Repair.Scars..0.1.,
                          height = specific.characteristics$Height..mm.,
                          PC1 = PCA$x[,1],
                          PC2 = PCA$x[,2],
                          PC3 = PCA$x[,3])
# Make Scree Plot
## Calculate loadings
prop <- data.frame(as.list(PCA$d/sum(PCA$d))) %>%
  pivot_longer(cols = 1:89)

## Produce and preview plot
scree_plot <- ggplot(data = prop[1:9,], aes(x = name, y = value)) +
  geom_col(fill = "blue3") +
  xlab("Principal Component") +
  ylab("Proportion of Variance Explained") +
  ggtitle("A") +
  theme_test()
scree_plot


# PC1 and PC2
## Create function to calculate the convex hull for each group
find_hull1 <- function(df) df[chull(df$PC1, df$PC2), ]

## Apply the function to your PCA scores grouped by your factor
PCA_Hulls1 <- PCA_Results %>%
  group_by(Location) %>%
  do(find_hull1(.))

## Plot PC Axes 1 and 2
PCA_Plot1 <- ggplot() +
  geom_polygon(data = PCA_Hulls1, aes(x = PC1,
                                      y = PC2,
                                      fill = Location), alpha = 0.2) + # Convex hulls 
  geom_point(data = PCA_Results, aes(x = PC1,
                                     y = PC2,
                                     color = Location,
                                     shape = Location), size = 3) + # Points
  xlab("PC1: 32.9%") +
  ylab("PC2: 15.2%") +
  ggtitle("B") +
  theme_test()
PCA_Plot1


# PC1 and PC3
## Create function to calculate the convex hull for each group
find_hull2 <- function(df) df[chull(df$PC1, df$PC3), ]

## Apply the function to your PCA scores grouped by your factor
PCA_Hulls2 <- PCA_Results %>%
  group_by(Location) %>%
  do(find_hull2(.))

## Plot PC Axes 1 and 3
PCA_Plot2 <- ggplot() +
  geom_polygon(data = PCA_Hulls2, aes(x = PC1,
                                      y = PC3, 
                                      fill = Location), alpha = 0.2) + # Convex hulls 
  geom_point(data = PCA_Results, aes(x = PC1,
                                     y = PC3,
                                     color = Location,
                                     shape = Location), size = 3) + # Points
  xlab("PC1: 32.9%") +
  ylab("PC3: 11.7%") +
  ggtitle("C") +
  theme_test()


# Compile master plot with PCs and scree plot
master_plot <- grid.arrange(scree_plot, PCA_Plot1, PCA_Plot2, layout_matrix = matrix(c(1, 2, 1, 3), nrow = 2))
ggsave("Results/PCA.png", master_plot, height = 7, width = 9, units = "in")


# Create tps deformation plot to see how shape changes along axes
## Set reference shape
ref <- mshape(Procrustes$coords)

## Save as svg
png(filename = "Results/tpsDeformations.png", height = 700, width = 500, units = "px", res = 150)

## Set up plot layout
par(mfrow=c(3,2), family="serif", mai=c(0.1,0.1,0.1,0.1))

## Produce tps deformations
plotRefToTarget(PCA$shapes$shapes.comp1$min, ref, method="TPS") #PC1 minimum value
plotRefToTarget(PCA$shapes$shapes.comp1$max, ref, method="TPS") #PC1 maximum value
plotRefToTarget(PCA$shapes$shapes.comp2$min, ref, method="TPS") #PC2 minimum value
plotRefToTarget(PCA$shapes$shapes.comp2$max, ref, method="TPS") #PC2 maximum value
plotRefToTarget(PCA$shapes$shapes.comp3$min, ref, method="TPS") #PC3 minimum value
plotRefToTarget(PCA$shapes$shapes.comp3$max, ref, method="TPS") #PC3 maximum value

## Turn graphics device off
dev.off()





# Check assumptions of PC1 for one-way ANOVA
ggplot() +
  geom_boxplot(data = PCA_Results, aes (x = Location, y = PC1, fill = Location))

# Too many outliers. Perform Kruskal-Wallis test
kruskal.test(PC1 ~ Location, data = PCA_Results)
pairwise.wilcox.test(PCA_Results$PC1, PCA_Results$Location,
                     p.adjust.method = "holm")




# Compile master plot with PCs and scree plot
master_tps <- grid.arrange(PC1Min, PC1Max, PC2Min, PC2Max, PC3Min, PC3Max,
                           layout_matrix = matrix(c(1, 2, 3, 4, 5, 6), nrow = 2))



PC1Min













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