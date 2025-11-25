# Call packages
library("geomorph")
library("stringr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("gridExtra")

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
summary(PCA)

# First three PCs explain the most variation. Create dataframe with first three PCs
PCA_Results <- data.frame(Location = specific.characteristics$Location,
                          RepairScars = specific.characteristics$Repair.Scars..0.1.,
                          height = specific.characteristics$Height..mm.,
                          PC1 = PCA$x[,1],
                          PC2 = PCA$x[,2],
                          PC3 = PCA$x[,3])

# Make scree plot
prop <- data.frame(as.list(PCA$d/sum(PCA$d))) %>%
  pivot_longer(cols = 1:89)

scree_plot <- ggplot(data = prop[1:9,], aes(x = name, y = value)) +
  geom_col(fill = "blue3") +
  xlab("Principal Component") +
  ylab("Proportion of Variance Explained") +
  ggtitle("Scree Plot") +
  theme_test()
scree_plot

# create function to calculate the convex hull for each group
find_hull <- function(df) df[chull(df$PC1, df$PC2), ]

# apply the function to your PCA scores grouped by your factor
PCA_Hulls2 <- PCA_Results %>%
  group_by(Location) %>%
  do(find_hull(.))

# Plot PC Axes 1 and 2
PCA_Plot2 <- ggplot() +
  geom_polygon(data = PCA_Hulls2, aes(x = PC1,
                                     y = PC2, 
                                     fill = Location), alpha = 0.2) + # Convex hulls 
  geom_point(data = PCA_Results, aes(x = PC1,
                                   y = PC2,
                                   color = Location,
                                   shape = Location), size = 3) + # Points
  xlab("PC1: 32.9%") +
  ylab("PC2: 15.2%") +
  theme_test()
PCA_Plot2


# create function to calculate the convex hull for each group
find_hull <- function(df) df[chull(df$PC1, df$PC3), ]

# apply the function to your PCA scores grouped by your factor
PCA_Hulls3 <- PCA_Results %>%
  group_by(Location) %>%
  do(find_hull(.))

# Plot PC Axes 1 and 3
PCA_Plot3 <- ggplot() +
  geom_polygon(data = PCA_Hulls3, aes(x = PC1,
                                      y = PC3, 
                                      fill = Location), alpha = 0.2) + # Convex hulls 
  geom_point(data = PCA_Results, aes(x = PC1,
                                     y = PC3,
                                     color = Location,
                                     shape = Location), size = 3) + # Points
  xlab("PC1: 32.9%") +
  ylab("PC3: 11.7%") +
  theme_test()

# Compile master plot with PCs and scree plot
master_plot <- grid.arrange(scree_plot, PCA_Plot2, PCA_Plot3, layout_matrix = matrix(c(1, 2, 1, 3), nrow = 2))
ggsave("PCA.png", master_plot, height = 4.5, width = 9, units = "in")



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