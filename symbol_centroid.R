# =====================================================================
# Script Name:  symbol_centroid.R
# Author:       Adlin Pinheiro
# Date:         09/28/2026
# Version:      1.0
#
# Description:
# Computes the spatial centroid(s) of a dCDT symbol using the convex hull 
# of drawn points. For non-clock-face symbols (digits, hands), returns a 
# single centroid. For the clock face (CF or CF'), subdivides the drawing
# into 4 angular quadrants (NE, NW, SW, SE based on centroid), computes 
# a centroid for each quadrant, and also returns the overall clock-face 
# center. Used within dcdt_to_adj_matrix() to extract 19-node network 
# node positions for pairwise distance calculations.
#
# INPUTS:
#   df (data.frame): Parsed dCDT pen stroke data with columns:
#       symbollabel, symboltype, x, y, scaled_x, scaled_y (i.e.
#       output from read_csk.R)
#
# OUTPUTS:
#   A data frame with one row per centroid:
#       - symbollabel (character): Symbol identifier (e.g., "1", "CF_Q1", 
#         "CF_CENTER")
#       - symboltype (character): Symbol type (e.g., "ONE", "CLOCKFACE_Q1",
#         "CLOCKFACE_CENTER")
#       - centroid_x, centroid_y (numeric): Centroid coordinates; 
#         NA_real_ if insufficient points
#
# ARGUMENTS:
#   df (data.frame): dCDT stroke data
#   symbollabel (character): Symbol to process (e.g., "1", "2", "CF", "CF'")
#   scale (logical): if TRUE, uses scaled_x/scaled_y; if FALSE, uses raw x/y
#
# CENTROID CALCULATION:
#   (1) Subset df to symbollabel; extract x, y coordinates
#   (2) If <3 points or insufficient variation: return (NA, NA)
#   (3) Create owin window (either [0,1]×[0,1] for scaled or data bounds)
#   (4) Create ppp (point pattern) and compute convex hull
#   (5) Extract area-weighted centroid via centroid.owin()
#
# CLOCK-FACE SPECIAL HANDLING (symbollabel in "CF", "CF'"):
#   (1) Compute overall centroid (cen_x, cen_y)
#   (2) For each point, compute angle θ from centroid:
#       θ = 180 - (180/π) * atan2(y - cen_y, x - cen_x)
#       (0° = E, 90° = N, 180° = W, 270° = S)
#   (3) Assign quadrants: Q1 [0°-90°], Q2 [90°-180°], Q3 [180°-270°], 
#       Q4 [270°-360°]
#   (4) Compute centroid for each present quadrant
#   (5) Return quadrant centroids + overall clock-face center (5 rows total)
#
# DEPENDENCIES:
#   - spatstat (owin, ppp, convexhull, centroid.owin)
#   - Used within: dcdt_to_adj_matrix.R
#
# RETURNS:
#   - 1 row for non-clock-face symbols (digit or hand)
#   - Up to 5 rows for clock face (4 quadrants + 1 center); fewer if 
#     quadrants are empty
#   - All centroid coordinates are NA_real_ if insufficient data
# =====================================================================

library(spatstat)

symbol_centroid <- function(df, 
                            symbollabel,
                            scale = TRUE) {
  
  # Subset to the particular clock symbol
  symbol_dcdt_df <- df[df$symbollabel == symbollabel, ]
  
  if (scale) {
    x <- symbol_dcdt_df$scaled_x
    y <- symbol_dcdt_df$scaled_y
  } else {
    x <- symbol_dcdt_df$x
    y <- symbol_dcdt_df$y
  }
  
  # Helper function to calculate centroid
  get_centroid <- function(x, y, scale = TRUE) {
    
    # If too few points for convex hull, use mean coordinate
    if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
      return(c(NA_real_, NA_real_))
    }
    
    if (scale) {
      win <- owin(c(0, 1), c(0, 1))
    } else {
      win <- owin(
        c(min(x), max(x)),
        c(min(y), max(y))
      )
    }
    
    cen <- centroid.owin(
      convexhull(
        ppp(
          x,
          y,
          window = win,
          checkdup = FALSE
        )
      )
    )
    
    c(as.numeric(cen[1]), as.numeric(cen[2]))
  }
  
  # Get centroid of the symbol
  centroid <- get_centroid(x, y, scale = scale)
  cen_x <- centroid[1]
  cen_y <- centroid[2]
  
  # If symbol label is the clock face, separate the clock face into 4 quadrants
  if (symbollabel %in% c("CF", "CF'")) {
    
    symbol_dcdt_df$angle <- 180 - (180 / pi) * atan2(
      y - cen_y,
      x - cen_x
    )
    
    symbol_dcdt_df$quadrant <- NA_integer_
    
    symbol_dcdt_df$quadrant[
      symbol_dcdt_df$angle >= 0 & symbol_dcdt_df$angle < 90
    ] <- 1
    
    symbol_dcdt_df$quadrant[
      symbol_dcdt_df$angle >= 90 & symbol_dcdt_df$angle < 180
    ] <- 2
    
    symbol_dcdt_df$quadrant[
      symbol_dcdt_df$angle >= 180 & symbol_dcdt_df$angle < 270
    ] <- 3
    
    symbol_dcdt_df$quadrant[
      symbol_dcdt_df$angle >= 270 & symbol_dcdt_df$angle <= 360
    ] <- 4 
    
    quads <- sort(unique(na.omit(symbol_dcdt_df$quadrant)))
    
    # Initialize empty data frame
    cf_quads <- data.frame(
      symbollabel = character(),
      symboltype = character(),
      centroid_x = numeric(),
      centroid_y = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(quads)) {
      
      q <- quads[i]
      quad_dcdt_df <- symbol_dcdt_df[symbol_dcdt_df$quadrant == q, ]
      
      if (scale) {
        xq <- quad_dcdt_df$scaled_x
        yq <- quad_dcdt_df$scaled_y
      } else {
        xq <- quad_dcdt_df$x
        yq <- quad_dcdt_df$y
      }
      
      cen <- get_centroid(xq, yq, scale = scale)
      
      cf_quads <- rbind(
        cf_quads,
        data.frame(
          symbollabel = paste0("CF_Q", q),
          symboltype = paste0("CLOCKFACE_Q", q),
          centroid_x = cen[1],
          centroid_y = cen[2],
          stringsAsFactors = FALSE
        )
      )
    }
    
    # Append clock-face center
    cf_center_row <- data.frame(
      symbollabel = "CF_CENTER",
      symboltype = "CLOCKFACE_CENTER",
      centroid_x = cen_x,
      centroid_y = cen_y,
      stringsAsFactors = FALSE
    )
    
    out <- rbind(cf_quads, cf_center_row)
    
  } else {
    
    out <- data.frame(
      symbollabel = symbollabel,
      symboltype = as.character(symbol_dcdt_df$symboltype)[1],
      centroid_x = cen_x,
      centroid_y = cen_y,
      stringsAsFactors = FALSE
    )
  }
  
  return(out)
}

# Example usage:
# file <- "/path/to/dcdt_example.csk"
# symbol_centroid(read_csk(file), symbollabel = "2", scale = TRUE)
# symbol_centroid(read_csk(file), symbollabel = "CF", scale = TRUE)
