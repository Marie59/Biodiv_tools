#Rscript

##########################
##    Leaf area index   ##
##########################

#####Packages : leafR
#               mapview

#####Load arguments

library(leafR)
library(mapview)

#####Import data
#setwd("/home/pndb/Scriptr_MJ/data/L1A/L1A/localisation/shape/")

#table1 <- read.table("data_FrenchBBS.tabular", sep = "\t", dec = ".", header = TRUE)
normlas.file = system.file("extdata", "lidar_example.laz", package="leafR")
#colnames(table1) <- c("site", "annee", "espece", "abond")

#####Your analysis

####Leaf area index####

# Calculate LAD from voxelization
VOXELS_LAD = lad.voxels(normlas.file,
                        grain.size = 2)

# Calculate the LAD profile
lad_profile = lad.profile(VOXELS_LAD)

lidar.lai = lai(lad_profile); lidar.lai
understory.lai = lai(lad_profile, min = 1, max = 5); understory.lai

# relative LAD PROFILE
relative.lad_profile = lad.profile(VOXELS_LAD, relative = TRUE)

#understory relative LAI (% of total LAI)
relative.understory.lai = lai(relative.lad_profile, min = 1, max = 5); relative.understory.lai

#mapping fr visualization
#Calculate LAI derived from LAD profile
lidar.lai = lai(lad_profile)
VOXELS_LAD.5 = lad.voxels(normlas.file,
                        grain.size = 5)

#Map using relative values (%)
relative.lai_raster = lai.raster(VOXELS_LAD.5, relative.value = lidar.lai)
mapview(relative.lai_raster)
