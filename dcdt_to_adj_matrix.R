# =====================================================================
# Script Name:  dcdt_to_adj_matrix.R
# Author:       Adlin Pinheiro
# Date:         09/28/2026
# Version:      1.0
#
# Description:
# Computes a kernel-weighted adjacency matrix for a dCDT network by:
# (1) extracting centroids for 19 symbol nodes (14 digits/hands + 
#     4 clock-face quadrants + 1 clock-face center)
# (2) computing pairwise Euclidean distances between centroids
# (3) applying a user-specified kernel transformation to distances
#
# Missing symbols (not drawn on the test) appear as rows/columns of NAs.
# The resulting matrix is suitable for graph construction and network 
# analysis of dCDT drawings.
#
# INPUTS:
#   dcdt_df (data.frame): Parsed dCDT pen stroke data with columns for 
#       (x, y) coordinates, symbollabel, symboltype, and drawing type
#
# OUTPUTS:
#   A 19×19 matrix of kernel-transformed edge weights (or raw distances 
#   if return.dist.matrix=TRUE). Row and column names are symboltype labels.
#   Missing nodes are represented by NA rows/columns in the distance matrix
#   and 0 in the weighted adjacency matrix.
#
# ARGUMENTS:
#   condition (character): "COMMAND" or "COPY"; specifies which test to process
#   scale (logical): if TRUE, uses scaled_x/scaled_y; if FALSE, uses raw x/y 
#        coordinates           
#   return.dist.matrix (logical): if TRUE, returns pairwise Euclidean distances; 
#       if FALSE, returns kernel-transformed weights
#   kernel (character): kernel type - one of "exponential", "gaussian", 
#       "inverse", or "linear"
#   lambda (numeric): decay rate for exponential kernel (default: 1.0);
#       higher values = faster decay with distance
#   sigma (numeric): bandwidth for Gaussian kernel (default: 0.5);
#       lower values = sharper decay; typically in (0.2, 1.0)
#
# DEPENDENCIES:
#   - spatstat (convexhull, ppp, owin, centroid.owin)
#   - dplyr (group_by, summarize)
#   - geometry (convex hull computations)
#   - source files: cf_source.R, symbol_centroid.R
#
# =====================================================================

library(spatstat) 
library(parallel)
library(dplyr)
library(geometry)

source("/restricted/projectnb/dnpstats/adlinp/Project 2/Scripts/Paper GitHub/cf_source.R")
source("/restricted/projectnb/dnpstats/adlinp/Project 2/Scripts/Paper GitHub/symbol_centroid.R")
source("/restricted/projectnb/dnpstats/adlinp/Project 2/Scripts/Paper GitHub/read_csk.R")

dcdt_to_adj_matrix <- function(dcdt_df, 
                                condition  = c("COMMAND", "COPY"),
                                scale = TRUE, return.dist.matrix = FALSE,
                                kernel = c("exponential", "gaussian", "inverse", "linear"),
                                lambda = 1,
                                sigma = 0.5) {
  
  if (is.null(dcdt_df) || nrow(dcdt_df) == 0) return(NULL)
  
  condition <- match.arg(condition)
  
  kernel <- match.arg(kernel)
  
  dcdt_df <- dcdt_df[dcdt_df$drawing == condition, ]  
  
  if (is.null(dcdt_df) || nrow(dcdt_df) == 0) return(NULL)
  
  source <- cf_source(dcdt_df)
  
  # Get centroids for all relevant symbol labels except CF
  symbols <- c("ONE", "TWO", "THREE", "FOUR", "FIVE", "SIX", "SEVEN", "EIGHT", 
               "NINE", "TEN", "ELEVEN", "TWELVE", "HOUR_HAND", "MINUTE_HAND")
  
  if (!any(levels(dcdt_df$symboltype) %in% symbols) & 
      !any(levels(dcdt_df$symboltype) == "CLOCKFACE")) {
    return(NULL)  
  }
  
  # Subset to just digits and hands
  digits_dcdt_df <- dcdt_df[dcdt_df$symboltype %in% symbols, ]
  
  # Initialize empty data frame
  digit_centroids <- data.frame(symbollabel = character(),
                                symboltype = character(),
                                centroid_x = numeric(),
                                centroid_y = numeric(),
                                stringsAsFactors = FALSE)
  
  # Get data frame of centroid coordinates for each symbol label
  if (nrow(digits_dcdt_df) > 0) { 
    
    unique_labels <- as.character(unique(droplevels(digits_dcdt_df$symbollabel)))
    
    for (i in seq_along(unique_labels)) {
      
      result <- symbol_centroid(digits_dcdt_df,
                                symbollabel = unique_labels[i],
                                scale = scale)
      
      digit_centroids <- rbind(digit_centroids, result)
      
    } # Ends unique labels loop
  } # Ends if nrow > 0
  
  # Get centroid coordinates all 5 CF elements (center, Q1-Q4)
  # Append to the results of the other clock symbols
  if (source %in% c("CF", "CF'")) {
    
    cf <- dcdt_df[dcdt_df$symbollabel == source, ]
    digit_cf_centroids <- rbind(digit_centroids, symbol_centroid(cf,
                                                              symbollabel = source,
                                                              scale = scale))
  } else { # If no "CF", return the center coordinates of the entire drawing
    
    if (scale) {
      ch <- convexhull(ppp(dcdt_df$scaled_x, 
                           dcdt_df$scaled_y, 
                           window=owin(c(0,1), c(0,1)),
                           checkdup=FALSE))
      
    } else {
      ch <- convexhull(ppp(dcdt_df$x, 
                           dcdt_df$y, 
                           window=owin(c(min(dcdt_df$x), max(dcdt_df$x)),
                                       c(min(dcdt_df$y), max(dcdt_df$y))),
                           checkdup=FALSE))
    }
    
    centroid <- centroid.owin(ch)
    
    digit_cf_centroids <- rbind(digit_centroids, data.frame(symbollabel = "CF_CENTER",
                                                         symboltype = "CLOCKFACE_CENTER",
                                                         centroid_x = centroid$x,
                                                         centroid_y = centroid$y))
  }
  
  # Make sure the symboltype is a factor with 19 levels
  digit_cf_centroids$symboltype <- factor(digit_cf_centroids$symboltype,
                                       levels = c(symbols, "CLOCKFACE_Q1", "CLOCKFACE_Q2",
                                                  "CLOCKFACE_Q3", "CLOCKFACE_Q4", 
                                                  "CLOCKFACE_CENTER"))
  
  digit_cf_centroids$centroid_x <- as.numeric(digit_cf_centroids$centroid_x)
  digit_cf_centroids$centroid_y <- as.numeric(digit_cf_centroids$centroid_y)
  
  # If any symbol label is repeated, get the average x and y coordinate
  # Removes column symbollabel
  digit_cf_centroids <- digit_cf_centroids %>%
    group_by(symboltype) %>%
    summarize(centroid_x = mean(centroid_x, na.rm = TRUE),
              centroid_y = mean(centroid_y, na.rm = TRUE))
  
  # Check for missing symbol factor levels
  missing_levels <- setdiff(levels(digit_cf_centroids$symboltype), unique(digit_cf_centroids$symboltype))
  
  # If there are missing symbols, add them to the data frame of centroid 
  # coordinates with NA values
  if (length(missing_levels) > 0) {
    missing_dcdt_df <- data.frame(symboltype = missing_levels, 
                                  centroid_x = NA, 
                                  centroid_y = NA)
    
    # Add the missing rows to the original data frame
    digit_cf_centroids <- rbind(digit_cf_centroids, missing_dcdt_df)
  }
  
  # Sort rows
  digit_cf_centroids <- digit_cf_centroids[order(digit_cf_centroids$symboltype), ]
  
  # Obtain distance matrix (matrix of pairwise distances between each symbol label)
  # Distance will be NA if a node isn't drawn on the test
  dist_mat <- as.matrix(dist(as.matrix(digit_cf_centroids[, c("centroid_x", "centroid_y")])))
  colnames(dist_mat) <- digit_cf_centroids$symboltype
  rownames(dist_mat) <- digit_cf_centroids$symboltype
  
  if (return.dist.matrix) {
    return(dist_mat)
  }
  
  # Convert distances to edge weights and output network adjacency matrix
  
  if (kernel == "exponential") {
    
    out_mat <- exp(-lambda * dist_mat) 
    
  } else if (kernel == "linear") {
    
    
    if (scale) {
      # Max possible distance is sqrt(2) in a 1x1 square
      out_mat <- sqrt(2) - dist_mat
    } else {

      x_range <- max(dcdt_df$x, na.rm = TRUE) - min(dcdt_df$x, na.rm = TRUE)
      y_range <- max(dcdt_df$y, na.rm = TRUE) - min(dcdt_df$y, na.rm = TRUE)
      max_dist <- sqrt(x_range^2 + y_range^2)
      out_mat <- max_dist - dist_mat
      
    }
    
  } else if (kernel == "gaussian") {
    
    out_mat <- exp(- (dist_mat^2) / (2 * sigma^2))
    
  } else if (kernel == "inverse") {
    out_mat <-  1/(1 + dist_mat)  
  }
  
  out_mat[is.na(out_mat)] <- 0        # Missing distances are set to 0 (no edge)
  diag(out_mat) <- 0                  # Self-loops are 0
  
  return(out_mat)
}

# Example usage:
# file <- "/path/to/dcdt_example.csk"
# dcdt_network <- dcdt_to_adj_matrix(read_csk(file))