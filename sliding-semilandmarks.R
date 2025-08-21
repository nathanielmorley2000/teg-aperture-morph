# Call packages
library("geomorph")

# Upload TPS file with semi-landmark data
landmarks <- readland.tps(file = "Data/TegulaOutlines - Copy.TPS", specID = "ID", negNA = FALSE,
                          readcurves = TRUE, warnmsg = TRUE)

# Upload specimen characteristics data
characteristics <- read.csv("Data/TegulaCharacteristicsList.csv")

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


