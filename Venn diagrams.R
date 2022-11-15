###VENN DIAGRAMS## F 
#####################
# Make Venn Diagrams Functions
#note that your data must be TRANSPOSED
# Source Niels's scripts
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram2.r")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")
#join metaG and MetaT table with metadata descriptors
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")

#metaG<-read.csv("MetaG_tranposed.csv")
metaT<-read.csv("MetaT_net-trap_transposed.csv")
eExt<-read.csv("Extraction_test_ASV_avg.csv")
ext<-read.csv("ext_test_global_Orders_SILVA.csv")

asv<-t(asv)
ext<-ext[-4,]
Ext<-t(eExt)
Ext<-as.data.frame(Ext)
#for asvs:

#metaT
long<-gather(metaT, Pathway, RPKM, X12DICHLORETHDEG.PWY:XYLCAT.PWY )
long<-as.data.table(ext_test_long)
long<-gather(ext, Order, Count, Unknown:Xanthomonadales )

#make into data.table
long<-as.data.table(long)
some<-subset(long, Count >0)
a<-unique(some$Order)
a
#MetaG:
long<-gather(metaG, Pathway, RPKM, X12DICHLORETHDEG.PWY:VALSYN.PWY )

metaT2<-as.data.table(long)
#groups<-read.csv("metaG_groups_isa.csv")
groups2<-read.csv("isa_MetaT_groups.csv")
groups2<-read.csv("Ext_test_groups.csv")
groups<-as.data.table(groups2)
#in case you want to make row names the column names
#tmetaG<-t(metaG)
#df <- tibble::rownames_to_column(tmetaG, "SAMPLE_ID")
#library(data.table)
#setDT(tmetaG, keep.rownames = TRUE)[]
setkey(metaT2, SAMPLE_ID)
setkey(groups, SAMPLE_ID)
setkey(long, Type)
setkey(groups, Type)

#join- 
data_metaG <- metaT2[groups]
data_test<-long[groups]
#remove data not part of the venn:
#data_metaG1<-data_metaG[, -c(4,5,8,9)]

# Make my own Venn diagram function for 2 sets as a check 
MakeVennDiagram2 <- function(set1, set2, categories, colors){
  shared <- intersect(set1,set2)
  only1 <- setdiff(set1,set2)
  only2 <- setdiff(set2,set1)
  # make plot
  plot.new()
  venn.plot <- draw.pairwise.venn(length(only1)+length(shared), length(only2)+length(shared), length(shared), category = categories,
                                  fill = colors)
  # return plot object (need to print to plot)
  return(venn.plot)
}



# Do using Niels function as a check ... 
MakeVennDiagramLPSI2 <- function(dt){
  venn_diagram2(a = unique(dt[Zone=="photic" & (RPKM>0)]$Pathway),
                b = unique(dt[Zone=="bathypelagic" & (RPKM>0)]$Pathway),
                name_a = "photic", 
                name_b = "bathypelagic", 
                colors=c("lightblue", "navy"), 
                euler =TRUE)
}
MakeVennDiagramLPSI2(dt = data_metaG)
dev.off()

# Set Difference Analysis - Depth - LP only
MakeVennDiagramLP_depths <- function(dt){
  venn_diagram3(a = unique(dt[Type=="M" & (Count>0)]$Order),
                b = unique(dt[Type=="S" & (Count>0)]$Order),
                c = unique(dt[Type=="K" & (Count>0)]$Order), 
                name_a = "M", 
                name_b = "S", 
                name_c = "K", 
                colors=c("lightblue","dodgerblue", "navy"))
}

MakeVennDiagramLP_depths(dt = long) #note that your data frame needs to be in data.table format!!!!
dev.off()
sum(long$Count)
count(long$M)

install.packages("BiocManager")
BiocManager::install("limma")
library(VennDiagram)
# Chart
long %>%                    # take the data.frame "data"   # Using "data", filter out all rows with NAs in aa 
  group_by(Type) %>%          # Then, with the filtered data, group it by "bb"
  summarise(Unique_Elements = n_distinct(Order)) 
#number of unique orders per extraction type
num_orders<-long %>%
  group_by(Type, Order) %>%
  dplyr::mutate(count = n()) %>% 
  unique()
write.csv(num_orders, "num_orders_unique.csv")

#using a presence absence

#convert matrix into presence/absemce
asv<-read.csv("Greenland_2018__ASV_tableb.csv", row.names=1)
asv[asv>0] <-1
#export
write.csv(asv, "presence_absence_asv.csv")
#sum across
