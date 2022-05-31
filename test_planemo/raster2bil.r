#Rscript

###################################
##    Converting raster to BIL   ##
###################################
#####Packages : stringr
#               raster
#               xfun


#####Load arguments

args <- commandArgs(trailingOnly = TRUE)


# url for the S2 subset

if (length(args) < 1) {
    stop("This tool needs at least 1 argument")
}else{
    data_raster <- args[1]
    sensor <- as.character(args[2])
}   



#' Reads ENVI hdr file
#'
#' @param HDRpath Path of the hdr file
#'
#' @return list of the content of the hdr file
#' @export

read_ENVI_header <- function(HDRpath) {
  # header <- paste(header, collapse = "\n")
  if (!grepl(".hdr$", HDRpath)) {
    stop("File extension should be .hdr")
  }
  HDR <- readLines(HDRpath)
  ## check ENVI at beginning of file
  if (!grepl("ENVI", HDR[1])) {
    stop("Not an ENVI header (ENVI keyword missing)")
  } else {
    HDR <- HDR [-1]
  }
  ## remove curly braces and put multi-line key-value-pairs into one line
  HDR <- gsub("\\{([^}]*)\\}", "\\1", HDR)
  l <- grep("\\{", HDR)
  r <- grep("\\}", HDR)

  if (length(l) != length(r)) {
    stop("Error matching curly braces in header (differing numbers).")
  }

  if (any(r <= l)) {
    stop("Mismatch of curly braces in header.")
  }

  HDR[l] <- sub("\\{", "", HDR[l])
  HDR[r] <- sub("\\}", "", HDR[r])

  for (i in rev(seq_along(l))) {
    HDR <- c(
      HDR [seq_len(l [i] - 1)],
      paste(HDR [l [i]:r [i]], collapse = "\n"),
      HDR [-seq_len(r [i])]
    )
  }

  ## split key = value constructs into list with keys as names
  HDR <- sapply(HDR, split_line, "=", USE.NAMES = FALSE)
  names(HDR) <- tolower(names(HDR))

  ## process numeric values
  tmp <- names(HDR) %in% c(
    "samples", "lines", "bands", "header offset", "data type",
    "byte order", "default bands", "data ignore value",
    "wavelength", "fwhm", "data gain values"
  )
  HDR [tmp] <- lapply(HDR [tmp], function(x) {
    as.numeric(unlist(strsplit(x, ",")))
  })

  return(HDR)
}


#' writes ENVI hdr file
#'
#' @param HDR list. content to be written
#' @param HDRpath character. Path of the hdr file
#'
#' @return None
#' @importFrom stringr str_count
#' @export

write_ENVI_header <- function(HDR, HDRpath) {
  h <- lapply(HDR, function(x) {
    if (length(x) > 1 || (is.character(x) && stringr::str_count(x, '\\w+') > 1)) {
      x <- paste0('{', paste(x, collapse = ','), '}')
    }
    # convert last numerics
    x <- as.character(x)
  })
  writeLines(c('ENVI', paste(names(HDR), h, sep = ' = ')), con = HDRpath)
  return(invisible())
}


#' get hdr name from image file name, assuming it is BIL format
#'
#' @param ImPath character. ath of the image
#' @param showWarnings boolean. set TRUE if warning because HDR does not exist
#'
#' @return corresponding hdr
#' @import tools
#' @export

get_HDR_name <- function(ImPath,showWarnings=TRUE) {
  if (xfun::file_ext(ImPath) == "") {
    ImPathHDR <- paste(ImPath, ".hdr", sep = "")
  } else if (xfun::file_ext(ImPath) == "bil") {
    ImPathHDR <- gsub(".bil", ".hdr", ImPath)
  } else if (xfun::file_ext(ImPath) == "zip") {
    ImPathHDR <- gsub(".zip", ".hdr", ImPath)
  } else {
    ImPathHDR <- paste(tools::file_path_sans_ext(ImPath), ".hdr", sep = "")
  }
  if (showWarnings==TRUE){
    if (!file.exists(ImPathHDR)) {
      message("WARNING : COULD NOT FIND HDR FILE")
      print(ImPathHDR)
      message("Process may stop")
    }
  }
  return(ImPathHDR)
}

#' converts a raster into BIL format as expected by biodivMapR codes
#'
#' @param Raster_Path character. Full path for the raster to be converted
#' @param Sensor character. Name of the sensor. a .hdr template for the sensor should be provided in extdata/HDR
#' @param Convert_Integer boolean. Should data be converted into integer ?
#' @param Multiplying_Factor numeric. Multiplying factor (eg convert real reflectance values between 0 and 1 into integer between 0 and 10000).
#' @param Output_Dir character. Path to output directory.
#' @param Multiplying_Factor_Last numeric. Multiplying factor for last band.
#' @param Mask boolean is the file a mask file with 0s and 1s only?
#'
#' @return Output_Path path for the image converted into ENVI BIL format
#' @importFrom raster brick writeRaster
#' @importFrom tools tools::file_path_sans_ext
#' @export
raster2BIL <- function(Raster_Path, Sensor = "unknown", Output_Dir = FALSE, Convert_Integer = TRUE,
                       Multiplying_Factor = 1, Multiplying_Factor_Last = 1, Mask = FALSE) {

  # get directory and file name of original image
  Input_File <- basename(Raster_Path)
  Input_Dir <- dirname(Raster_Path)
  # define path where data will be stored
  if (Output_Dir == FALSE) {
    Output_Path <- file.path(Input_Dir, "biodivMapR_Convert_BIL", tools::file_path_sans_ext(Input_File))
  } else {
    dir.create(Output_Dir, showWarnings = FALSE, recursive = TRUE)
    Output_Path <- file.path(Output_Dir, tools::file_path_sans_ext(Input_File))
  }
  message("The converted file will be written in the following location:")
  print(Output_Path)
  Output_Dir <- dirname(Output_Path)
  dir.create(Output_Dir, showWarnings = FALSE, recursive = TRUE)

  # apply multiplying factors
  message("reading initial file")
  Output_Img <- Multiplying_Factor * raster::brick(Raster_Path)
  Last_Band_Name <- Output_Img@data@names[length(Output_Img@data@names)]
  Output_Img[[Last_Band_Name]] <- Multiplying_Factor_Last * Output_Img[[Last_Band_Name]]

  # convert into integer if requested or if is Mask
  if (Convert_Integer == TRUE | Mask == TRUE) {
    Output_Img <- round(Output_Img)
  }

  # write raster
  message("writing converted file")
  if (Convert_Integer == TRUE & Mask == FALSE) {
    r <- raster::writeRaster(Output_Img, filename = Output_Path, format = "EHdr", overwrite = TRUE, datatype = "INT2S")
  } else if (Convert_Integer == TRUE & Mask == TRUE) {
    r <- raster::writeRaster(Output_Img, filename = Output_Path, format = "EHdr", overwrite = TRUE, datatype = "INT1U")
  } else {
    r <- raster::writeRaster(Output_Img, filename = Output_Path, format = "EHdr", overwrite = TRUE)
  }
  hdr(r, format = "ENVI")

  # remove unnecessary files
  File2Remove <- paste(Output_Path, ".aux.xml", sep = "")
  File2Remove2 <- paste(tools::file_path_sans_ext(Output_Path), ".aux.xml", sep = "")
  file.remove(File2Remove)

  File2Remove <- paste(Output_Path, ".prj", sep = "")
  File2Remove2 <- paste(tools::file_path_sans_ext(Output_Path), ".prj", sep = "")
  if (file.exists(File2Remove)) {
    file.remove(File2Remove)
  } else if (file.exists(File2Remove2)) {
    file.remove(File2Remove2)
  }

  File2Remove <- paste(Output_Path, ".sta", sep = "")
  File2Remove <- paste(tools::file_path_sans_ext(Output_Path), ".sta", sep = "")
  if (file.exists(File2Remove)) {
    file.remove(File2Remove)
  } else if (file.exists(File2Remove2)) {
    file.remove(File2Remove2)
  }

  File2Remove <- paste(Output_Path, ".stx", sep = "")
  File2Remove2 <- paste(tools::file_path_sans_ext(Output_Path), ".stx", sep = "")
  if (file.exists(File2Remove)) {
    file.remove(File2Remove)
  } else if (file.exists(File2Remove2)) {
    file.remove(File2Remove2)
  }

  File2Rename <- paste(tools::file_path_sans_ext(Output_Path), ".hdr", sep = "")
  File2Rename2 <- paste(Output_Path, ".hdr", sep = "")
  if (file.exists(File2Rename)) {
    file.rename(from = File2Rename, to = File2Rename2)
  }

  # change dot into underscore
  Output_Path_US <- file.path(
    dirname(Output_Path),
    gsub(basename(Output_Path), pattern = "[.]", replacement = "_")
  )
  if (!Output_Path_US == Output_Path) {
    file.rename(from = Output_Path, to = Output_Path_US)
  }

  Output_Path_US_HDR <- paste0(Output_Path_US, ".hdr")
  if (!Output_Path_US_HDR == paste0(Output_Path, ".hdr")) {
    file.rename(from = paste0(Output_Path, ".hdr"), to = Output_Path_US_HDR)
    ### UTILITY?? ###
    file.rename(from = Output_Path_US_HDR, to = Output_Path_US_HDR)
  }

  if (!Sensor == "unknown") {
    HDR_Temp_Path <- system.file("extdata", "HDR", paste0(Sensor, ".hdr"), package = "biodivMapR")
    if (file.exists(HDR_Temp_Path)) {
      message("reading header template corresponding to the sensor located here:")
      print(HDR_Temp_Path)
      # get raster template corresponding to the sensor
      HDR_Template <- read_ENVI_header(HDR_Temp_Path)
      # get info to update hdr file
      # read hdr
      HDR_input <- read_ENVI_header(get_HDR_name(Output_Path))
      if (!is.null(HDR_Template$wavelength)) {
        HDR_input$wavelength <- HDR_Template$wavelength
      }
      if (!is.null(HDR_Template$`sensor type`)) {
        HDR_input$`sensor type` <- HDR_Template$`sensor type`
      }
      if (!is.null(HDR_Template$`band names`)) {
        HDR_input$`band names` <- HDR_Template$`band names`
      }
      if (!is.null(HDR_Template$`wavelength units`)) {
        HDR_input$`wavelength units` <- HDR_Template$`wavelength units`
      }
      # define visual stretch in the VIS domain
      HDR_input$`default stretch` <- "0 1000 linear"
      # write corresponding hdr file
      write_ENVI_header(HDR_input, get_HDR_name(Output_Path))
    } else if (!file.exists(HDR_Temp_Path)) {
      message("Header template corresponding to the sensor expected to be found here")
      print(HDR_Temp_Path)
      message("please provide this header template in order to write info in HDR file")
      print(get_HDR_name(Output_Path))
      message("or manually add wavelength location in HDR file, if relevant")
    }
  } else if (Sensor == "unknown") {
    message("please make sure that the following header file contains information required")
    print(get_HDR_name(Output_Path))
    message("or manually add wavelength location in HDR file, if relevant")
  }
  return(Output_Path)
}


raster2bil <- raster2BIL(Raster_Path = data_raster, Sensor = sensor)

cat("\nwrite files. \n--> \"", paste(raster2bil, "\"\n", sep = ""), file = "Metadata.hdr", sep = "", append = TRUE)



