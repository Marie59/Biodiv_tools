library(obisindicators)
library(dplyr)
library(dggridR)
library(sf)
library(ggplot2)

#Get biological occurrences
#Use the 1 million records subsampled from the full OBIS dataset
occ <- occ_1M # occ_1M OR occ_SAtlantic

#Create a discrete global grid
#Create an ISEA discrete global grid of resolution 9 using the dggridR package:

dggs <- dgconstruct(projection = "ISEA", topology = "HEXAGON", res = 9)

#Then assign cell numbers to the occurrence data
occ$cell <- dgGEO_to_SEQNUM(dggs, occ$decimalLongitude, occ$decimalLatitude)[["seqnum"]]

#Calculate indicators
#The following function calculates the number of records, species richness, Simpson index, Shannon index, Hurlbert index (n = 50), and Hill numbers for each cell.

#Perform the calculation on species level data
idx <- calc_indicators(occ)

#dd cell geometries to the indicators table (idx)
grid <- dgcellstogrid(dggs, idx$cell) %>% 
  st_wrap_dateline() %>% 
  rename(cell = seqnum) %>% 
  left_join(
    idx,
    by = "cell")

#Plot maps of indicators
#Let’s look at the resulting indicators in map form.
# ES(50)
es_50 <- gmap_indicator(grid, "es", label = "ES(50)")
es_50
# Shannon index
shannon_map <- gmap_indicator(grid, "shannon", label = "Shannon index")
shannon_map

# Number of records, log10 scale, Robinson projection (default)
records_map_robin <- gmap_indicator(grid, "n", label = "# of records", trans = "log10")
records_map_robin

# Number of records, log10 scale, Geographic projection
records_map_geo <- gmap_indicator(grid, "n", label = "# of records", trans = "log10", crs=4326)
records_map_geo

# Simpson index
simpson_map <- gmap_indicator(grid, "simpson", label = "Simpson index")
simpson_map

# maxp
maxp_map <- gmap_indicator(grid, "maxp", label = "maxp index")
maxp_map
