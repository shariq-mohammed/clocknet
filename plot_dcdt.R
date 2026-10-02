# =====================================================================
# Script Name:  plot_dcdt.R
# Author:       Adlin Pinheiro
# Date:         09/28/2026
# Version:      1.0
#
# Description:
# Visualizes a single dCDT drawing (COMMAND or COPY) as a point pattern.
# Extracts pen stroke coordinates, creates a spatial point pattern (ppp),
# and rotates by -90 degrees to correct drawing orientation (clock face
# with 12 o'clock at top). Optionally displays raw or scaled coordinates.
#
# INPUTS:
#   dcdt_df (data.frame): Parsed dCDT pen stroke data from read_csk()
#       with columns: drawing, x, y, scaled_x, scaled_y, etc.
#
# OUTPUTS:
#   None (plots to current graphics device)
#
# ARGUMENTS:
#   dcdt_df (data.frame): dCDT data frame
#   condition (character): "COMMAND" or "COPY"; specifies which test to plot
#   scale (logical): if TRUE, plots scaled coordinates [0,1]×[0,1];
#       if FALSE, plots raw pixel coordinates with automatic window bounds
#   rotation_angle (numeric): Rotation angle in radians; default -pi/2 
#       (90° clockwise) to align clock face with 12 o'clock at top.
#       Use 0 for no rotation, pi/2 for 90° counterclockwise, etc.
#   main (character): Plot title; passed to plot(); default NULL
#   ... : Additional arguments passed to plot() (e.g., pch, cex, col)
#
# DEPENDENCIES:
#   - spatstat.geom (ppp, owin, rotate.ppp)
#   - source file: read_csk.R
#
# NOTES:
#   - If condition does not exist in dcdt_df, empty plot is returned
#   - Window bounds are [0,1]×[0,1] for scaled; data-dependent for raw
# =====================================================================

plot_dcdt <- function(dcdt_df, condition = c("COMMAND", "COPY"), 
                             scale = TRUE, main = NULL,
                             rotation_angle = -pi/2,
                             ...) {
  
  dcdt_df <- dcdt_df[dcdt_df$drawing == condition, ]
  
  if (scale == FALSE){
    x <- dcdt_df$x
    y <- dcdt_df$y
    pp <- spatstat.geom::rotate.ppp(spatstat.geom::ppp(x, y, 
                                                       window=owin(c(min(x), max(x)),
                                                                   c(min(y), max(y))),
                                                       checkdup=FALSE), 
                                    angle=rotation_angle,
                                    center="centroid")
  } else {
    x <- dcdt_df$scaled_x
    y <- dcdt_df$scaled_y
    pp <- spatstat.geom::rotate.ppp(spatstat.geom::ppp(x, y, 
                                                       window=owin(c(0, 1),
                                                                   c(0, 1)),
                                                       checkdup=FALSE), 
                                    angle=rotation_angle, 
                                    center="centroid")
  }
  
  plot(pp, main = NULL, ...)
  
}

# Example usage:
# file <- "/path/to/dcdt_example.csk"
# plot_dcdt(read_csk(file))


