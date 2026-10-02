# =====================================================================
# Script Name:  plot_clocknet.R
# Author:       Adlin Pinheiro
# Date:         09/28/2026
# Version:      1.0
#
# Description:
# qgraph-based network visualization for the 19-node dCDT network. Creates 
# a fixed clock-based spatial layout and plots a 19 x 19 adjacency matrix
# from dcdt_to_adj_matrix().
#
# INPUTS:
#   adj_mat:
#       A 19 x 19 adjacency matrix with rows/columns ordered as:
#       digits 1-12, HH, MH, CF Q1-CF Q4, and CF CENTER.
#
#   ...:
#       Additional arguments passed directly to qgraph().
#
# OUTPUTS:
#   A qgraph network plot using the fixed dCDT spatial layout.
#
# LAYOUT:
#   Clock digits (1-12): Arranged on circle of radius 0.75, with 12 at top
#       - Positions: angles [0°, 30°, ..., 330°] in standard orientation
#       - Then reordered so 12 is at angle 90° (12 o'clock position)
#   Hour hand (HH): (-0.3, 0.2)
#   Minute hand (MH): (0.3, 0.2)
#
#   Clock-face quadrants and center:
#       - Q1 (NE): (1, 1)
#       - Q2 (SE): (1, -1)
#       - Q3 (SW): (-1, -1)
#       - Q4 (NW): (-1, 1)
#       - CTR (Center): (0, 0)
#
# DEPENDENCIES:
#   - qgraph (qgraph network visualization function)
#   - dcdt_to_adj_matrix.R (adjacency matrix computation)
#   - read_csk.R (data loading)
#
# NOTES:
#   - This script is intended for qgraph visualization.
#   - The layout is fixed and reusable across multiple dCDT drawings.
#   - Node labels are shortened for clarity, e.g., "CTR" instead of
#     "CF CENTER".
#   - Edge visibility and styling can be controlled using qgraph()
#     arguments such as cut, minimum, edge.color, and edge.width.
#   - The adjacency matrix must follow the expected 19-node ordering.
# =====================================================================

plot_clocknet <- function(adj_mat, ...) {
  
  # Layout matrix for network nodes
  radius <- 0.75
  angles_radians <- seq(0, 330, by = 30) * pi / 180
  x_coords <- radius * cos(angles_radians)
  y_coords <- radius * sin(angles_radians)
  
  clock_coords <- cbind(x_coords, y_coords)
  clock_coords <- clock_coords[nrow(clock_coords):1, ]
  clock_coords <- clock_coords[c(10:12, 1:9), ]
  
  other_coords <- matrix(c(
    -0.3,  0.2,  # Node HH
    0.3,  0.2,   # Node MH
    1,    1,     # Q1
    1,   -1,     # Q2
    -1,   -1,    # Q3
    -1,    1,    # Q4
    0,    0      # CF Center
  ), ncol = 2, byrow = TRUE)
  
  layout_matrix <- rbind(clock_coords, other_coords)
  
  # Shortened label names for visualization purposes
  node_names <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
                  "HH", "MH", "CF Q1", "CF Q2", "CF Q3", "CF Q4", "CTR")
  
  colnames(adj_mat) <- node_names 
  rownames(adj_mat) <- node_names
  
  qgraph::qgraph(adj_mat, layout = layout_matrix, ...)
}

# Example Usage
# file <- "/path/to/dcdt_example.csk"
# out_mat <- dcdt_to_adj_matrix(read_csk(file))
# plot_clocknet(out_mat, details = TRUE, cut = 0.8)
