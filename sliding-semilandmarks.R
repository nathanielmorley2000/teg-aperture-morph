# Call packages
library("geomorph")
library("RRPP")
library("stringr")
library("tibble")
library("ggplot2")
library("gridExtra")
library("plotly")
library("htmlwidgets")
library("openxlsx")
library("broom")
library("dplyr")



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

# Produce simple PCA (no size correction)
PCA <- gm.prcomp(GPA.Results$coords)

# Create eigenvalues table for principal components and save to CSV file
PCA_summary <- summary(PCA)
PCA_summary <- as.data.frame(t(PCA_summary$PC.summary)) %>%
  rownames_to_column(var = "Principal.Component") %>%
  mutate(Principal.Component = str_replace_all(Principal.Component, "Comp", "PC")) %>%
  mutate(across(-Principal.Component, ~ signif(., 3)))
write.csv(PCA_summary, "Results/PrincipalComponents.csv", row.names = FALSE)

# First three PCs explain the most variation. Create dataframe with first three PCs
PCA_Results <- data.frame(Location = specific.characteristics$Location,
                          RepairScars = specific.characteristics$Repair.Scars..0.1.,
                          height = specific.characteristics$Height..mm.,
                          PC1 = PCA$x[,1],
                          PC2 = PCA$x[,2],
                          PC3 = PCA$x[,3])

# Reorder plot to align with gradients
PCA_Results <- PCA_Results %>%
  mutate(Location = factor(Location, levels = c("Strawberry", "Mathers", "Prasiola"))) %>%
  mutate(RepairScars = ifelse(RepairScars == 1, "Present", "Absent"))



## Plot ordination =============================================================

### Study locality =============================================================

# Create three-dimensional PCA plot
interactive_PCA <- plot_ly(data = PCA_Results,
                           x = ~PC1,
                           y = ~PC2,
                           z = ~PC3,
                           color = ~Location,
                           colors = c("#f9766e", "#000000", "#00bdbf"),
                           symbol = ~Location,
                           symbols = c(18, 19, 15),
                           marker = list(size = 5),
                           type = "scatter3d",
                           mode = "markers") %>%
  layout(legend = list(orientation = "h",
                       xanchor = "center",
                       x = 0.5,
                       itemsizing = "constant"), 
         scene = list(xaxis = list(title = "PC1 (32.9%)",
                                   dtick = 0.05,
                                   tickformat = ".2f",
                                   autorange = "reversed"),
                      yaxis = list(title = "PC2 (15.2%)",
                                   nticks = 5,
                                   tickformat = ".2f",
                                   autorange = "reversed"),
                      zaxis = list(title = "PC3 (11.7%)",
                                   nticks = 5,
                                   tickformat = ".2f"),
                      camera = list(eye = list(x = -1.25, y = 2.00, z = 0.75))))
interactive_PCA

# Save 3D plot as interactive HTML file
saveWidget(as_widget(interactive_PCA), "Results/MorphologyResults/StudyLocality/Interactive_PCA.html")

# Add title to PCA for compound figure
static_PCA <- interactive_PCA %>% layout(title = list(text = "A",
                                                      x = 0.2,
                                                      y = 0.8))

# Save 3D plot as static SVG file
orca(static_PCA, "Results/MorphologyResults/StudyLocality/Static_PCA.svg", scale = 1.5)

# Plot PC1 individually
PC1_Plot <- ggplot(PCA_Results, aes(x = PC1, y = "", color = Location, shape = Location)) +
  geom_point(position = position_jitter(width = 0, height = 0.15), size = 4) +
  scale_x_continuous(breaks = seq(-0.10, 0.15, by = 0.05)) +
  scale_color_manual(values = c("Strawberry" = "#f9766e", "Mathers" = "#000000", "Prasiola" = "#00bdbf")) +
  scale_shape_manual(values = c("Strawberry" = 18, "Mathers" = 19, "Prasiola" = 15)) +
  labs(title = "B", x = "PC1 (32.9%)", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.y.left = element_blank(),
        axis.ticks.y = element_blank())
ggsave("Results/MorphologyResults/StudyLocality/PC1.svg", PC1_Plot, width = 10, height = 2) # Save 1D plot

# Plot PC2 individually
PC2_Plot <- ggplot(PCA_Results, aes(x = PC2, y = "", color = Location, shape = Location)) +
  geom_point(position = position_jitter(width = 0, height = 0.15), size = 4) +
  scale_x_continuous(breaks = seq(-0.10, 0.10, by = 0.05)) +
  scale_color_manual(values = c("Strawberry" = "#f9766e", "Mathers" = "#000000", "Prasiola" = "#00bdbf")) +
  scale_shape_manual(values = c("Strawberry" = 18, "Mathers" = 19, "Prasiola" = 15)) +
  labs(title = "C", x = "PC2 (15.2%)", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.y.left = element_blank(),
        axis.ticks.y = element_blank())
ggsave("Results/MorphologyResults/StudyLocality/PC2.svg", PC2_Plot, width = 10, height = 2) # Save 1D plot

# Plot PC3 individually
PC3_Plot <- ggplot(PCA_Results, aes(x = PC3, y = "", color = Location,  shape = Location)) +
  geom_point(position = position_jitter(width = 0, height = 0.15), size = 4) +
  scale_x_continuous(breaks = seq(-0.10, 0.10, by = 0.05)) +
  scale_color_manual(values = c("Strawberry" = "#f9766e", "Mathers" = "#000000", "Prasiola" = "#00bdbf")) +
  scale_shape_manual(values = c("Strawberry" = 18, "Mathers" = 19, "Prasiola" = 15)) +
  labs(title = "D", x = "PC3 (11.7%)", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.y.left = element_blank(),
        axis.ticks.y = element_blank())
ggsave("Results/MorphologyResults/StudyLocality/PC3.svg", PC3_Plot, width = 10, height = 2) # Save 1D plot



### Repair scars ===============================================================

# Create three-dimensional PCA plot
interactive_PCA <- plot_ly(data = PCA_Results,
                           x = ~PC1,
                           y = ~PC2,
                           z = ~PC3,
                           color = ~RepairScars,
                           colors = c("#f9766e", "#00bdbf"),
                           symbol = ~RepairScars,
                           symbols = c(19, 15),
                           marker = list(size = 5),
                           type = "scatter3d",
                           mode = "markers") %>%
  layout(legend = list(orientation = "h",
                       xanchor = "center",
                       x = 0.5,
                       itemsizing = "constant"), 
         scene = list(xaxis = list(title = "PC1 (32.9%)",
                                   dtick = 0.05,
                                   tickformat = ".2f",
                                   autorange = "reversed"),
                      yaxis = list(title = "PC2 (15.2%)",
                                   nticks = 5,
                                   tickformat = ".2f",
                                   autorange = "reversed"),
                      zaxis = list(title = "PC3 (11.7%)",
                                   nticks = 5,
                                   tickformat = ".2f"),
                      camera = list(eye = list(x = -1.25, y = 2.00, z = 0.75))))
interactive_PCA

# Save 3D plot as interactive HTML file
saveWidget(as_widget(interactive_PCA), "Results/MorphologyResults/RepairScars/Interactive_PCA.html")

# Add title to PCA for compound figure
static_PCA <- interactive_PCA %>% layout(title = list(text = "A",
                                                      x = 0.2,
                                                      y = 0.8))

# Save 3D plot as static SVG file
orca(static_PCA, "Results/MorphologyResults/RepairScars/Static_PCA.svg", scale = 1.5)

# Plot PC1 individually
PC1_Plot <- ggplot(PCA_Results, aes(x = PC1, y = "", color = RepairScars, shape =RepairScars)) +
  geom_point(position = position_jitter(width = 0, height = 0.15), size = 4) +
  scale_x_continuous(breaks = seq(-0.10, 0.15, by = 0.05)) +
  scale_color_manual(values = c("Absent" = "#f9766e", "Present" = "#00bdbf")) +
  scale_shape_manual(values = c("Absent" = 19, "Present" = 15)) +
  labs(title = "B", x = "PC1 (32.9%)", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.y.left = element_blank(),
        axis.ticks.y = element_blank())
ggsave("Results/MorphologyResults/RepairScars/PC1.svg", PC1_Plot, width = 10, height = 2) # Save 1D plot

# Plot PC2 individually
PC2_Plot <- ggplot(PCA_Results, aes(x = PC2, y = "", color = RepairScars, shape = RepairScars)) +
  geom_point(position = position_jitter(width = 0, height = 0.15), size = 4) +
  scale_x_continuous(breaks = seq(-0.10, 0.10, by = 0.05)) +
  scale_color_manual(values = c("Absent" = "#f9766e",  "Present" = "#00bdbf")) +
  scale_shape_manual(values = c("Absent" = 19, "Present" = 15)) +
  labs(title = "C", x = "PC2 (15.2%)", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.y.left = element_blank(),
        axis.ticks.y = element_blank())
ggsave("Results/MorphologyResults/RepairScars/PC2.svg", PC2_Plot, width = 10, height = 2) # Save 1D plot

# Plot PC3 individually
PC3_Plot <- ggplot(PCA_Results, aes(x = PC3, y = "", color = RepairScars,  shape = RepairScars)) +
  geom_point(position = position_jitter(width = 0, height = 0.15), size = 4) +
  scale_x_continuous(breaks = seq(-0.10, 0.10, by = 0.05)) +
  scale_color_manual(values = c("Absent" = "#f9766e", "Present" = "#00bdbf")) +
  scale_shape_manual(values = c("Absent" = 19, "Present" = 15)) +
  labs(title = "D", x = "PC3 (11.7%)", y = "") +
  theme_classic() +
  theme(legend.position = "none",
        axis.line.y.left = element_blank(),
        axis.ticks.y = element_blank())
ggsave("Results/MorphologyResults/RepairScars/PC3.svg", PC3_Plot, width = 10, height = 2) # Save 1D plot



## TPS deformations ============================================================

# Set reference shape
ref <- mshape(Procrustes$coords)

# Produce tps deformations
## PC1 Min
svg(filename = "Results/MorphologyResults/TPS_Deformations/PC1min.svg")
plotRefToTarget(ref, PCA$shapes$shapes.comp1$min, method="TPS") 
dev.off()

## PC1 Max
svg(filename = "Results/MorphologyResults/TPS_Deformations/PC1max.svg")
plotRefToTarget(ref, PCA$shapes$shapes.comp1$max, method="TPS")
dev.off()

## PC2 Min
svg(filename = "Results/MorphologyResults/TPS_Deformations/PC2min.svg")
plotRefToTarget(ref, PCA$shapes$shapes.comp2$min, method="TPS")
dev.off()

## PC2 Max
svg(filename = "Results/MorphologyResults/TPS_Deformations/PC2max.svg")
plotRefToTarget(ref, PCA$shapes$shapes.comp2$max, method="TPS")
dev.off()

## PC3 Min
svg(filename = "Results/MorphologyResults/TPS_Deformations/PC3min.svg")
plotRefToTarget(ref, PCA$shapes$shapes.comp3$min, method="TPS")
dev.off()

## PC3 Max
svg(filename = "Results/MorphologyResults/TPS_Deformations/PC3max.svg")
plotRefToTarget(ref, PCA$shapes$shapes.comp3$max, method="TPS")
dev.off()



## Procrustes ANOVA ============================================================

# Create two-factor Procrustes ANOVA with interaction term
interactive <- procD.lm(coords ~ location * repairs, data = GPA.Results, RRPP = TRUE, print.progress = FALSE)
summary_interactive <- summary(interactive, formula = FALSE)

# Compile and save Results into one file
statistical_results <- list("Two-factor ANOVA" = tidy(summary_interactive$table))
write.xlsx(statistical_results, file = "Results/MorphologyResults/StatisticalResults.xlsx")



# ALLOMETRIC VARIATION =========================================================

# Produce simple allomtery model
simpleAllom <- procD.lm(coords ~ log(Csize), data = GPA.Results, print.progress = FALSE)

# Produce common allometry models
commonAllom1 <- procD.lm(coords ~ log(Csize) + location, data = GPA.Results,
                         print.progress = FALSE)
commonAllom2 <- procD.lm(coords ~ log(Csize) + repairs, data = GPA.Results,
                         print.progress = FALSE)
commonAllom3 <- procD.lm(coords ~ log(Csize) + (location + repairs), data = GPA.Results,
                         print.progress = FALSE)
commonAllom4 <- procD.lm(coords ~ log(Csize) + (location * repairs), data = GPA.Results,
                         print.progress = FALSE)

# Produce unique allometries model
uniqueAllom1 <- procD.lm(coords ~ log(Csize) * location, data = GPA.Results,
                 print.progress = FALSE)
uniqueAllom2 <- procD.lm(coords ~ log(Csize) * repairs, data = GPA.Results,
                         print.progress = FALSE)
uniqueAllom3 <- procD.lm(coords ~ log(Csize) * (location + repairs), data = GPA.Results,
                         print.progress = FALSE)
uniqueAllom4 <- procD.lm(coords ~ log(Csize) * (location * repairs), data = GPA.Results,
                         print.progress = FALSE)

# Compare models to find which is most likely
modComp <- model.comparison(simpleAllom, commonAllom1, commonAllom2, commonAllom3, commonAllom4,
                            uniqueAllom1, uniqueAllom2, uniqueAllom3, uniqueAllom4,
                            type = "logLik")
modCompData <- summary(modComp)
write.csv(modCompData, "Results/AllometryResults/modCompData.csv")

# Model output for common allometry model against predation history (most likely)
simpleAllometryData <- summary(simpleAllom)$table
write.csv(simpleAllometryData, "Results/AllometryResults/simpleAllometryData.csv", row.names = TRUE)

# Create allometry plot using default graphics
allometryPlot <- plotAllometry(simpleAllom, size = GPA.Results$Csize, logsz = TRUE, method = "RegScore",
                               pch = 19)

# Create data frame with allometry plot for ggplot
allometryPlotData <- data.frame(Predictor = allometryPlot[["plot_args"]][["x"]],
                                Regression.Scores = allometryPlot[["plot_args"]][["y"]]) %>%
  rownames_to_column(var = "Specimen")
  
# Create nice allometry plot
allometryLinePlot <- ggplot(data = allometryPlotData, aes(x = Predictor,
                                                          y = Regression.Scores)) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_point(size = 3) +
  xlab("log(Centroid Size)") +
  ylab("Standardized Shape Scores") +
  theme_test() +
  theme(legend.position = "right")

ggsave("Results/AllometryResults/allometryLinePlot.svg", allometryLinePlot)

