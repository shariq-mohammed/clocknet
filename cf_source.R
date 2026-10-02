# =====================================================================
# Script Name:  cf_source.R
# Author:       Adlin Pinheiro
# Date:         09/28/2026
# Version:      1.0
#
# Description:
# Selects which clock face to use for feature extraction when multiple 
# clock faces (CF and/or CF') are drawn on a dCDT. Uses spatial reasoning 
# (convex hull containment) to determine whether digits/hands are drawn 
# within CF, CF', both, or neither. Returns the most appropriate clock 
# face based on spatial overlap. For use in dcdt_to_adj_matrix.R to 
# determine which CF to use for network creation.
#
# INPUTS:
#   dcdt_df (data.frame): Parsed dCDT pen stroke data from read_csk()
#       with columns: drawing, symbollabel, symboltype, scaled_x, scaled_y
#
# OUTPUTS:
#   A character/logical value indicating which clock face to use:
#       - NULL: No useable data (all points are noise/unclassified)
#       - NA: Cannot determine; either no CF drawn, or ambiguous overlap
#       - "CF": First drawn clock face selected
#       - "CF'": Second drawn clock face selected
#
# ARGUMENTS:
#   dcdt_df (data.frame): dCDT stroke data containing all drawings
#   condition (character): "COMMAND" or "COPY"; specifies which test to process
#
# LOGIC:
#   (1) Filter to specified condition (COMMAND or COPY)
#   (2) If no CF present: return NA
#   (3) If only CF or only CF' present: return that one
#   (4) If both CF and CF' present:
#       a. Check for spatial overlap between CF and CF' point sets
#       b. If overlapping: return "CF'" (most recent)
#       c. If non-overlapping:
#          - Extract digits/hands (exclude CF, CF', NOISE, NCD)
#          - Check which CF contains the digits/hands via convex hull
#          - If only CF has internal points: return "CF"
#          - If only CF' has internal points: return "CF'"
#          - If both have internal points: compare count of distinct elements
#            and return the CF with more elements (ties → "CF'")
#          - If neither has internal points: return NA with warning
#   (5) Return NA with warning if ambiguous or no valid data
#
# DEPENDENCIES:
#   - geometry (convhulln, inhulln)
#
# NOTES:
#   - Does not handle scenario with 3+ clock faces
#   - Returns NA (not an error) for ambiguous cases with warnings
#   - Uses convex hull containment; sensitive to extreme outlier points
#   - Scaled coordinates [0,1]×[0,1] assumed for spatial comparisons
# =====================================================================

library(geometry) 

cf_source <- function(dcdt_df, condition = c("COMMAND", "COPY")) {
  
  condition <- match.arg(condition)
  
  dcdt_df <- dcdt_df[dcdt_df$drawing == condition, ]  
  
  if (nrow(dcdt_df) == 0) return(NULL)
  
  symbollabels <- levels(droplevels(dcdt_df$symbollabel))
  
  which_cf <- NA  # Default initialization
  
  # Check if all points are noise/unclassified
  if (all(dcdt_df$symboltype %in% c("NOT_CLOCK_DATA", "UNCLASSIFIED", "NOISE"))) {
    which_cf <- NULL  # No useable data
    
    # No CF drawn
  } else if (!"CF" %in% symbollabels) {
    which_cf <- NA
    
    # CF only (no CF')
  } else if ("CF" %in% symbollabels & !"CF'" %in% symbollabels) {
    which_cf <- "CF"
    
    # CF' only (no CF) - may not ever exist
  } else if (!"CF" %in% symbollabels & "CF'" %in% symbollabels) {
    which_cf <- "CF'"
    
    # Both CF and CF' exist - decide which to use
  } else if ("CF" %in% symbollabels & "CF'" %in% symbollabels) {
    
    tryCatch({
      
      # Extract CF and CF' point coordinates
      cf_dcdt_df <- subset(dcdt_df, dcdt_df$symbollabel == "CF")
      cfp_dcdt_df <- subset(dcdt_df, dcdt_df$symbollabel == "CF'")
      
      cf <- as.matrix(cf_dcdt_df[, c("scaled_x", "scaled_y")])
      cfp <- as.matrix(cfp_dcdt_df[, c("scaled_x", "scaled_y")])
      
      # Compute convex hulls
      cf.ch <- convhulln(cf)
      cfp.ch <- convhulln(cfp)
      
      # Check for overlapping points between CF and CF'
      cfp_in_cf <- inhulln(cf.ch, cfp)
      cf_in_cfp <- inhulln(cfp.ch, cf)
      
      # If CF and CF' do not overlap
      if (sum(cfp_in_cf) == 0 & sum(cf_in_cfp) == 0) {
        
        # Restrict to digits and hands (exclude clock faces and noise)
        nocf <- dcdt_df[!(dcdt_df$symbollabel %in% c("CF", "CF'", "NOISE", "NCD")), ]
        
        # No digits/hands drawn
        if (nrow(nocf) == 0) {
          warning("Both CF and CF' drawn with no internal symbols; cannot determine which clock to use")
          which_cf <- NA
          
          # Check if coordinates are all NA
        } else if (all(is.na(nocf$scaled_x)) | all(is.na(nocf$scaled_y))) {
          warning("No valid coordinates in internal symbols for CF/CF' comparison")
          which_cf <- NA
          
        } else {
          
          # Check which CF contains the digits/hands
          in_cfpoints <- inhulln(cf.ch, as.matrix(nocf[, c("scaled_x", "scaled_y")]))
          in_cfppoints <- inhulln(cfp.ch, as.matrix(nocf[, c("scaled_x", "scaled_y")]))
          
          # CF has internal points, CF' does not
          if (sum(in_cfpoints) > 0 & sum(in_cfppoints) == 0) {
            which_cf <- "CF"
            
            # CF' has internal points, CF does not
          } else if (sum(in_cfpoints) == 0 & sum(in_cfppoints) > 0) {
            which_cf <- "CF'"
            
            # Both CF and CF' have internal points - compare number of distinct elements
          } else if (sum(in_cfpoints) > 0 & sum(in_cfppoints) > 0) {
            
            points_in_cf <- nocf[in_cfpoints, ]
            points_in_cfp <- nocf[in_cfppoints, ]
            
            n_in_cf_elements <- nlevels(droplevels(points_in_cf$symbollabel))
            n_in_cfp_elements <- nlevels(droplevels(points_in_cfp$symbollabel))
            
            if (n_in_cf_elements > n_in_cfp_elements) {
              which_cf <- "CF"
            } else if (n_in_cf_elements < n_in_cfp_elements) {
              which_cf <- "CF'"
            } else {
              # Equal number of elements - use second drawn clock
              which_cf <- "CF'"
            }
            
            # Neither CF nor CF' contains internal points
          } else if (sum(in_cfpoints) == 0 & sum(in_cfppoints) == 0) {
            warning("Digits/hands exist but are outside both CF and CF'")
            which_cf <- NA
          }
        }
        
        # CF and CF' overlap - use most recent (CF')
      } else if (sum(cfp_in_cf) > 0 | sum(cf_in_cfp) > 0) {
        which_cf <- "CF'"
      }
      
    }, error = function(e) {
      warning("Error during CF/CF' selection: ", conditionMessage(e))
      which_cf <<- NULL
    })
  }
  
  return(which_cf)
}

# Example usage:
# file <- "/path/to/dcdt_example.csk"
# cf_source(read_csk(file), condition = "COMMAND")
# cf_source(read_csk(file), condition = "COPY")
