# CLASSIFICATION OF ELECTRIC TOWERS AND DETECTION OF LINE’S

## The program uses LSD to detect Line segments
## The detected line segments are then merged to extract the Tower


## the implementation of line segment detection is lsd_call_example.c
## the implementation of LSD detector is in lsd.c

# To run

$	make
$	./lsd_call_example 2.pgm 5.eps 700 150

# After running following output files are generated
# parallels contains all parallel lines detected for each line
# singles contains all lines that can be merged into one
# final.eps contains all the merged lines and reconstructed tower

