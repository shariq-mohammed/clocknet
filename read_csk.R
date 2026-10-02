# =====================================================================
# Script Name:  read_csk.R
# Author:       Adlin Pinheiro
# Date:         09/28/2026
# Version:      1.0
#
# Description:
# Parses a dCDT drawing file (.csk XML format) and extracts pen stroke 
# data into a tabular format. Hierarchically structured XML is flattened 
# into a point-level data frame. Coordinates are normalized to [0,1] 
# based on the bounding box of drawn symbols (excluding NCD markers).
# Separate scaling is applied to COMMAND and COPY drawings.
#
# INPUTS:
#   xml_file (character): Path to a .csk file (XML format) containing 
#       dCDT pen stroke data
#
# OUTPUTS:
#   A data frame with one row per pen point (see header for full spec)
#
# ARGUMENTS:
#   xml_file (character): Path to .csk XML file
#   keep_drawings (character vector): Drawing types to retain; 
#       default c("COMMAND", "COPY"); set to NULL to keep all
#
# =====================================================================

library(stringr)
library(xml2)

read_csk <- function(xml_file, keep_drawings = c("COMMAND", "COPY")) {
  
  # Validate input file
  if (!file.exists(xml_file)) {
    stop(paste0("File not found: ", xml_file))
  }
  
  xml1 <- read_xml(xml_file)
  point_nodes <- xml_find_all(xml1, ".//point")
  
  # Validate XML structure
  if (length(point_nodes) == 0) {
    stop(paste0("No point nodes found in XML file: ", xml_file))
  }
  
  # Helper function to safely extract attributes
  safe_xml_attr <- function(nodes, xpath, attr_name) {
    result <- xml_attr(xml_find_first(nodes, xpath), attr_name)
    # Ensure NA values are explicit, not factor levels
    result
  }
  
  # Extract attributes with defensive checks
  drawing_vals <- safe_xml_attr(point_nodes, "ancestor::drawing", "type")
  symbollabel_vals <- safe_xml_attr(point_nodes, "ancestor::symbol", "label")
  symboltype_vals <- safe_xml_attr(point_nodes, "ancestor::symbol", "type")
  strokelabel_vals <- safe_xml_attr(point_nodes, "ancestor::stroke", "label")
  hooklettype_vals <- safe_xml_attr(point_nodes, "ancestor::hooklet", "hooklettype")
  
  x_vals <- as.numeric(xml_attr(point_nodes, "x"))
  y_vals <- as.numeric(xml_attr(point_nodes, "y"))
  pressure_vals <- as.numeric(xml_attr(point_nodes, "pressure"))
  timestamp_vals <- as.numeric(xml_attr(point_nodes, "timestamp"))
  
  # Create data frame with explicit NA handling for factors
  df <- data.frame(
    drawing = factor(drawing_vals, exclude = NA),
    symbollabel = factor(symbollabel_vals, exclude = NA),
    symboltype = factor(symboltype_vals, exclude = NA),
    strokelabel = factor(strokelabel_vals, exclude = NA),
    hooklettype = factor(hooklettype_vals, exclude = NA),
    x = x_vals,
    y = y_vals,
    pressure = pressure_vals,
    timestamp = timestamp_vals,
    stringsAsFactors = FALSE
  )
  
  # Remove rows where critical attributes are NA
  critical_na <- is.na(df$drawing) | is.na(df$x) | is.na(df$y)
  if (any(critical_na)) {
    warning(paste0("Removing ", sum(critical_na), " points with missing drawing/x/y attributes"))
    df <- df[!critical_na, ]
  }
  
  # Warnings for missing drawing types
  if (!is.null(keep_drawings)) {
    for (d in keep_drawings) {
      if (!(d %in% as.character(df$drawing))) {
        message(paste0("Warning: no ", tolower(d), " clock for ", xml_file))
      }
    }
  }
  
  if (nrow(df) == 0) {
    warning(paste0("No data remaining after filtering for ", xml_file))
    df$scaled_x <- numeric()
    df$scaled_y <- numeric()
    df$image_scale <- numeric()
    return(df)
  }
  
  # Scale coordinates by drawing level
  scale_one_drawing <- function(dat) {
    
    dat <- droplevels(dat)
    drawing_name <- unique(as.character(dat$drawing))
    
    # Initialize scaled coordinates and image scale column
    dat$scaled_x <- NA_real_
    dat$scaled_y <- NA_real_
    dat$image_scale <- NA_real_
    
    # Use non-NCD points for scaling
    index <- as.character(dat$symbollabel) != "NCD" &
      !is.na(dat$x) &
      !is.na(dat$y)
    
    if (any(index)) {
      
      x_range <- range(dat$x[index], na.rm = TRUE)
      y_range <- range(dat$y[index], na.rm = TRUE)
      
      x_diff <- diff(x_range)
      y_diff <- diff(y_range)
      image_scale <- max(x_diff, y_diff)
      
      # Store image scale as a column for all rows in this drawing
      dat$image_scale <- image_scale
      
      # Check for valid image scale
      if (!is.finite(image_scale) || image_scale <= 0) {
        message(paste0(
          "Warning: could not scale ", tolower(drawing_name),
          " clock (image_scale = ", image_scale, ") for ", xml_file
        ))
        return(dat)
      }
      
      x_mid <- mean(x_range)
      y_mid <- mean(y_range)
      
      dat$scaled_x[index] <- ((dat$x[index] - x_mid) / image_scale) + 0.5
      dat$scaled_y[index] <- ((dat$y[index] - y_mid) / image_scale) + 0.5
      
    } else {
      message(paste0(
        "Warning: no valid non-NCD points to scale for ", tolower(drawing_name),
        " clock in ", xml_file
      ))
    }
    
    return(dat)
  }
  
  # Scale separately by drawing level (COMMAND and COPY scaled independently)
  out_list <- lapply(split(df, df$drawing, drop = TRUE), scale_one_drawing)
  
  # Combine into one data frame
  out_df <- do.call(rbind, out_list)
  rownames(out_df) <- NULL
  
  return(out_df)
}

# Example usage:
# file <- "/path/to/dcdt_example.csk"
# out <- read_csk(file)