# clocknet
Build weighted spatial networks from digital clock drawing test (dCDT) pen data.

## Pre-processing functions

dcdt_example.csk: example XML file of clocksketch data

read_csk.R: reads one XML file of clocksketch data and outputs one data frame
of pen points

plot_dcdt.R: plots dCDT pen points into a spatial point pattern using the 
data frame output by read_csk.R. 

## Network functions

symbol_centroid.R: finds the geometric centroid of one dCDT clock component 
(i.e. column symbollabel in the data frame output by read_csk.R). Used in
dcdt_to_adj_matrix.R.

cf_source.R: determines which clock face to use for network construction in the
event that a clock face is drawn twice. Used in dcdt_to_adj_matrix.R.

dcdt_to_adj_matrix.R: converts the data frame of pen points output from
read_csk.R into a network adjacency matrix, transforming pairwise Euclidean
distances of clock components into edge weights. 

plot_clocknet.R: visualizes the dCDT network by plotting a qgraph
of the adjacency matrix output by dcdt_to_adj_matrix.R







