#Rscript

###########################################
##    Mapping alpha and beta diversity   ##
###########################################

#####Packages : obisindicators
#               dplyr
#               sf
#               ggplot2
#               rnaturalearth
#               rnaturalearthdata
#               viridis
#               dggridr

remotes::install_github("r-barnes/dggridR")
library(magrittr)
#####Load arguments

args <- commandArgs(trailingOnly = TRUE)

# url for the S2 subset

if (length(args) < 1) {
    stop("This tool needs at least 1 argument")
}else{
    raster <- args[1]
    hr <- args[2]
    crs <- as.numeric(args[3])
    source(args[4])
    source(args[5])
    source(args[6])
}   

if (hr == "false") {
  hr <- FALSE
}else{
  hr <- TRUE
}

if (crs == "0") {
   crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
}
#####Import data

#Get biological occurrences
#Use the 1 million records subsampled from the full OBIS dataset
occ <- read.table(raster, sep = ",", dec = ".", header = hr, fill = TRUE, encoding = "UTF-8") # occ_1M OR occ_SAtlantic
#colnames(occ) <- c("X", "decimalLongitude", "decimalLatitude", "species", "date_year", "records", "cell")
#Create a discrete global grid
#Create an ISEA discrete global grid of resolution 9 using the dggridR package:

dggs <- dggridR::dgconstruct(projection = "ISEA", topology = "HEXAGON", res = 9)

#Then assign cell numbers to the occurrence data
occ$cell <- dggridR::dgGEO_to_SEQNUM(dggs, occ$decimalLongitude, occ$decimalLatitude)[["seqnum"]]

#Calculate indicators
#The following function calculates the number of records, species richness, Simpson index, Shannon index, Hurlbert index (n = 50), and Hill numbers for each cell.

#Perform the calculation on species level data
idx <- calc_indicators(occ)
write.table(idx, file = "Index.csv", sep = ",", dec = ".", na = " ", col.names = T, row.names = F, quote = FALSE)


#dd cell geometries to the indicators table (idx)
grid <- dggridR::dgcellstogrid(dggs, idx$cell) %>% 
  sf::st_wrap_dateline() %>% 
  dplyr::rename(cell = seqnum) %>% 
  dplyr::left_join(
    idx,
    by = "cell")

#Plot maps of indicators
#Let’s look at the resulting indicators in map form.
# ES(50)
es_50_map <- gmap_indicator(grid, "es", label = "ES(50)")
es_50 <- ggplot2::ggsave("ES_50.png", es_50_map, scale = 0.38, width = 12, height = 7, units = "in", dpi = 300, limitsize = TRUE)

# Shannon index
shannon_map <- gmap_indicator(grid, "shannon", label = "Shannon index")
shannon <- ggplot2::ggsave("Shannon_index.png", shannon_map, scale = 0.38, width = 12, height = 7, units = "in", dpi = 300, limitsize = TRUE)


# Number of records, log10 scale, Geographic projection
records_map <- gmap_indicator(grid, "n", label = "# of records", trans = "log10", crs = crs)
records <- ggplot2::ggsave("Records.png", records_map, scale = 0.38, width = 12, height = 7, units = "in", dpi = 300, limitsize = TRUE)

# Simpson index
simpson_map <- gmap_indicator(grid, "simpson", label = "Simpson index")
simpson <- ggplot2::ggsave("Simpson_index.png", simpson_map, scale = 0.38, width = 12, height = 7, units = "in", dpi = 300, limitsize = TRUE)

# maxp
maxp_map <- gmap_indicator(grid, "maxp", label = "maxp index")
maxp <- ggplot2::ggsave("Maxp.png", maxp_map, scale = 0.38, width = 12, height = 7, units = "in", dpi = 300, limitsize = TRUE)

#Mapping
es_50
shannon
simpson
maxp
records
