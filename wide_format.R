library(dplyr)
library(tidyr)
#turn new into wide format (new is the name of my long dataset,
#idvar is the name of rows (you can have your sample names here) and timevar will be the name of columns- 
#(you can add the variable here- I think this is for NO3?)
wide<-reshape(new, idvar="Order", timevar="Sample", direction="wide")
#clean up colnames- you might notice that in the new wide dataset, R has added "count." 
#in front of your site name or variable name
#in the columns- clean up using:
for ( col in 1:ncol(wide)){
  colnames(wide)[col] <-  sub("Count.", "", colnames(wide)[col])
}