# Call packages
library("geomorph")
library("RRPP")
library("stringr")
library("dplyr")
library("tidyr")
library("tibble")
library("ggplot2")
library("gridExtra")
library("plotly")
library("htmlwidgets")



# MORPHOMETRICS ================================================================

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



# MOPRHOLOGICAL VARIATION ======================================================

## Make PCA ====================================================================

#Produce simple PCA (no size correction)
PCA <- gm.prcomp(GPA.Results$coords)
summary(PCA)

# First three PCs explain the most variation. Create dataframe with first three PCs
PCA_Results <- data.frame(Location = specific.characteristics$Location,
                          RepairScars = specific.characteristics$Repair.Scars..0.1.,
                          height = specific.characteristics$Height..mm.,
                          PC1 = PCA$x[,1],
                          PC2 = PCA$x[,2],
                          PC3 = PCA$x[,3])



## Plot ordination =============================================================

# Create three-dimensional PCA plot
interactive_PCA <- plot_ly(data = PCA_Results,
                           x = ~PC1,
                           y = ~PC2,
                           z = ~PC3,
                           color = ~Location,
                           type = "scatter3d",
                           mode = "markers")

# Save 3D plot as interactive HTML file
saveWidget(as_widget(interactive_PCA), "Interactive_PCA.html")

# Save 3D plot as static PDF
save_image(interactive_PCA, "Static_PCA.pdf")



## Perform statistics ==========================================================

# Perform a MANOVA on the three-dimensional PCA
res.main <- manova(cbind(PC1, PC2, PC3) ~ Location, data = PCA_Results)
MANOVA_Table <- summary(res.main) # Significant omnibus results (p << 0.001)
MANOVA_Table

# Perform ANOVA to identify which axes are significant
ANOVA_Table <- summary.aov(res.main) # Significant results on PC1 (p << 0.001) and PC3 (p << 0.001)


TukeyHSD(summary.aov(res.main))



as.matrix(PCA_Results[1:30,4:6])









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


## TPS deformations ============================================================
# Set reference shape
ref <- mshape(Procrustes$coords)

# Save as svg
png(filename = "Results/tpsDeformations.png", height = 700, width = 500, units = "px", res = 150)

# Set up plot layout
par(mfrow=c(3,2), family="serif", mai=c(0.1,0.1,0.1,0.1))

# Produce tps deformations
plotRefToTarget(ref, PCA$shapes$shapes.comp1$min, method="TPS") #PC1 minimum value
plotRefToTarget(ref, PCA$shapes$shapes.comp1$max, method="TPS") #PC1 maximum value
plotRefToTarget(ref, PCA$shapes$shapes.comp2$min, method="TPS") #PC2 minimum value
plotRefToTarget(ref, PCA$shapes$shapes.comp2$max, method="TPS") #PC2 maximum value
plotRefToTarget(ref, PCA$shapes$shapes.comp3$min, method="TPS") #PC3 minimum value
plotRefToTarget(ref, PCA$shapes$shapes.comp3$max, method="TPS") #PC3 maximum value

# Turn graphics device off
dev.off()

# Perform Kruskal-Wallis test to test significance along PC1
kruskal.test(PC1 ~ Location, data = PCA_Results)
pairwise.wilcox.test(PCA_Results$PC1, PCA_Results$Location,
                     p.adjust.method = "holm")

# Perform Kruskal-Wallis test to test significance along PC2
kruskal.test(PC2 ~ Location, data = PCA_Results)

# Perform Kruskal-Wallis test to test significance along PC3
kruskal.test(PC3 ~ Location, data = PCA_Results)
pairwise.wilcox.test(PCA_Results$PC3, PCA_Results$Location,
                     p.adjust.method = "holm")




# ALLOMETRIC VARIATION =========================================================

# Produce common allometry model
commonAllometry <- procD.lm(coords ~ log(height) + location, data = GPA.Results,
                print.progress = FALSE)

# Produce unique allometries model
uniqueAllometry <- procD.lm(coords ~ log(height) * location, data = GPA.Results,
                 print.progress = FALSE)

# Compare models to find which is most likely
modComp <- model.comparison(commonAllometry, uniqueAllometry, 
                            type = "logLik")
modCompData <- summary(modComp)
write.csv(modCompData, "Results/modCompData.csv")

# Model output for common allometry model (most likely)
commonAllometryData <- summary(commonAllometry)$table
write.csv(commonAllometryData, "Results/commonAllometryData.csv", row.names = TRUE)

# Create allometry plot using default graphics
allometryPlot <- plotAllometry(commonAllometry, size = GPA.Results$height, logsz = TRUE, method = "RegScore",
                               pch = 19, col = as.numeric(GPA.Results$location))

# Create data frame with allometry plot for ggplot
allometryPlotData <- data.frame(Predictor = allometryPlot[["plot_args"]][["x"]],
                                Regression.Scores = allometryPlot[["plot_args"]][["y"]]) %>%
  rownames_to_column(var = "Specimen") %>%
  add_column(Location = specific.characteristics$Location, .after = 1)
  
# Create nice allometry plot
allometryLinePlot <- ggplot(data = allometryPlotData, aes(x = Predictor,
                                                          y = Regression.Scores,
                                                          colour = Location)) +
  geom_smooth(method = "lm", se = FALSE, aes(linetype = Location)) +
  xlim(2.50, 3.50) +
  ylim(-0.15, 0.15) +
  xlab("log(Spire Height [mm])") +
  ylab("Standardized Shape Scores") +
  ggtitle("A") +
  theme_test() +
  theme(legend.position = "bottom",
        legend.title=element_blank())

allometryScatterPlot <- ggplot(data = allometryPlotData, aes(x = Predictor,
                                                          y = Regression.Scores,
                                                          colour = Location,
                                                          shape = Location)) +
  geom_point(size = 3) +
  xlim(2.50, 3.50) +
  ylim(-0.15, 0.15) +
  xlab("log(Spire Height [mm])") +
  ylab("Standardized Shape Scores") +
  ggtitle("B") +
  theme_test() +
  theme(legend.position = "bottom",
        legend.title=element_blank())

niceAllometryPlot <- grid.arrange(allometryLinePlot, allometryScatterPlot, layout_matrix = matrix(c(1, 2), nrow = 1))
niceAllometryPlot

ggsave("Results/AllometryResiduals.png", niceAllometryPlot,
       width = 10, height = 5, units = "in")
