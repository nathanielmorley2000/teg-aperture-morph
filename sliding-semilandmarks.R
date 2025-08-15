# Call packages
library("geomorph")

# Upload TPS file with semi-landmark data
landmarks <- readland.tps(file = "Data/TegulaOutlines.TPS", specID = "ID", negNA = FALSE,
                          readcurves = FALSE, warnmsg = TRUE)

# Upload specimen characteristics data
characteristics <- read.csv("Data/TegulaCharacteristicsList.csv")

# Upload sliders file
sliders <- read.csv("Data/Sliders.csv", header = FALSE)
