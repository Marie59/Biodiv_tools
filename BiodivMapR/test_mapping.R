library(mapview)
library(leafpop)
library(RColorBrewer)
#####ALPHA
alpha_path <- "/home/pndb/RESULTS/S2A_T33NUD_20180104_Subset/SPCA/ALPHA/Shannon_2_Fullres"
alpha_file <- raster(alpha_path)

alphamean_path <- "/home/pndb/RESULTS/S2A_T33NUD_20180104_Subset/SPCA/ALPHA/Shannon_2_MeanFilter_Fullres"
alphamean_file <- raster(alphamean_path)


mapview(alpha_file, layer.name = "Estimated Shannon index", col.regions = brewer.pal(7, "Dark2")) +   
  mapview(alphamean_file, legend = FALSE, col.regions = brewer.pal(7, "Dark2")) 



#####BETA
beta_path <- "/home/pndb/RESULTS/S2A_T33NUD_20180104_Subset/SPCA/BETA/BetaDiversity_BCdiss_PCO_2"
beta_file <- raster(beta_path)

betafull_path <- "/home/pndb/RESULTS/S2A_T33NUD_20180104_Subset/SPCA/BETA/BetaDiversity_BCdiss_PCO_2_Fullres"
betafull_file <- raster(betafull_path)

mapview(beta_file, layer.name = "PCoA", col.regions = c("red", "blue", "green")) + 
  mapview(betafull_file, legend = FALSE, col.regions = c("red", "blue", "green"))





alpha_path2 <-  "/home/pndb/Scriptr_MJ/data/S2A_T33NUD_Plots/Forest_HighDiversity.shp"
forest_highdivt <- read_sf(alpha_path2)

alpha_path3 <-  "/home/pndb/Scriptr_MJ/data/S2A_T33NUD_Plots/Forest_LowDiversity.shp"
forest_lowdivt <- read_sf(alpha_path3)

alpha_path4 <-  "/home/pndb/Scriptr_MJ/data/S2A_T33NUD_Plots/Forest_MediumDiversity.shp"
forest_mediumdivt <- read_sf(alpha_path4)

alpha_path5 <-  "/home/pndb/Scriptr_MJ/data/S2A_T33NUD_Plots/LowVegetation.shp"
lowveget <- read_sf(alpha_path5)



diversity_level <- left_join(forest_highdivt, forest_lowdivt, by = "id")

mapping <- mapview(forest_highdivt, col.regions = "purple") + 
  mapview(forest_lowdivt, col.regions = "blue") + 
  mapview(forest_mediumdivt, col.regions = "green") + 
  mapview(lowveget, col.regions = "yellow")
mapping
