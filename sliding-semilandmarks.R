# Call packages
library("geomorph")
library("stringr")

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






# Using Procrustes Distance for sliding
A <- gpagen(landmarks,
            curves = curves,
            ProcD = TRUE, print.progress = TRUE)

# Using bending energy for sliding
B <- gpagen(landmarks,
            curves = curves,
            ProcD = FALSE, print.progress = FALSE)

summary(A)
summary(B)


plotTangentSpace(A, axis1 = 1, axis2 = 2, warpgrids = TRUE, mesh = NULL, label = FALSE,
                 groups = NULL, legend = FALSE)

# Convert semi-landmark data into a 2-d array
twod_landmarks <- two.d.array(landmarks)


