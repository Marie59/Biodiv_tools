#Rscript

############################################
##  Validate ISO 19139 metadata documen   ##
############################################

#####Packages : ncdf4,
#               geometa,
#               httr
library(ncdf4)
library(XML)
#####Load arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("This tool needs at least 1 argument")
}else{
    input_data <- args[1]
}   

##------------------------------------------##
##      Read ISO 19139 from a file or url   ##
##------------------------------------------##

# Test depuis catalogue Indores http://indores-tmp.in2p3.fr/geonetwork/srv/fre/catalog.search#/metadata/112ebeea-e79c-422c-8a43-a5a8323b446b
# <!--ISO 19139 XML compliance: NO-->
input_data <- xml2::read_xml(input_data)

tmp <- tempfile(fileext = ".xml")

xml2::write_xml(input_data, tmp)

md <- geometa::readISO19139(tmp)


# validate iso
cat("\nValidation of metadata according to ISO 19139\n\n ", md$validate(), file = "Metadata_validation.txt", fill = 1, append = TRUE)
#validation <- capture.output(md[["validate"]])

#write.csv(validation, "Metadata_validation.csv", na = " ", row.names = F, col.names = F, quote = FALSE)
