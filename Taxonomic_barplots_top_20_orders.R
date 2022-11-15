library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(vegan)
library(tidyr)
#####hellinger transformation V6V8 hellinger transform
funct<-read.csv("Greenland_2018_metacyc_pwys_top_20_base.csv", row.names=1)
funct<-decostand(funct, method="hellinger")
comb<-read.csv("combined_V4V5_V6V8_asv_SILVA_low_count_removed.csv", row.names=1) # in v4v5 folder
read<-read.csv("QS_2021__ASV_table.csv", row.names=1)
#for Knorr
taxa<-read.csv("Phyla_knorr_16S.csv",row.names=1)
funct<-decostand(taxa, method="hellinger")
#read in extraction file
ext<-read.csv("Top_20_greengene_orders_ext_test.csv", row.names=1)
#select only for v6v8:
comb_v6v8<-comb[,c(28:49)]
#transpose
comb_v6v8<-t(comb_v6v8)
#hellinger transform:
comb_v6v8<-decostand(comb_v6v8, method="hellinger")
#write it out into csv file
write.csv(comb_v6v8, "hellinger_transform_v6v8_SILVA_low_removed.csv") #saved to v6v8 folder
#select top 20 in excel file
###calc rel abund for top 20####
v6v8<-read.csv("hellinger_top_20_SILVA_v6v8_low_removed.csv", row.names=1)#saved to v6v8 folder
#transpose
funct<-t(funct)
x<-funct/rowSums(funct)
x<-na.omit(x)
tx<-t(x)
write.csv(tx, "relative_abund_amalgamated_knorr_phyla.csv") #in v6v8 folder
#rel abund for extraction test
ext<-t(ext)
x<-ext/rowSums(ext)
x<-na.omit(x)
tx<-t(x)

funct<-t(funct)
x<-funct/rowSums(funct)
x<-na.omit(x)
tx<-t(x)
funct<-read.csv("relative_abundance_Top_20_base_pwys_metacyc.csv")

write.csv(x, "relative_abundance_Top_20_base_pwys_metacyc.csv") #in v6v8 folder
#calculate relative abundance

##use this for hellinger transformed rel abund:
v6v8_r<-read.csv("relative_abundance_hellinger_Top_20_orders_SILVA_v6v8_low_removed.csv")#transpose here or in excel (in v6v8 folder)
long<-read.csv("relative_abundance_Top_20_orders_greengenes_ext.csv")
long<-t(long)
long<-as.data.frame(long)
long<-long[-1,]
long<-read.csv("relative_abund_QS_2021.csv")
#for knorr----------------------------------------------------------------------------------------------------------------------------------------
taxa<-read.csv("relative_abund_amalgamated_knorr_phyla.csv",row.names=1)
taxa2<-t(taxa)
#make long
long<-gather(taxa, Sample_ID, rel_abund,KN_S15_1000m:KN_S7_Surface)
#get grouping info
grouping_info<-data.frame(row.names=rownames(taxa2),t(as.data.frame(strsplit(rownames(taxa2),"_"))))
#write.csv
write.csv(grouping_info, "grouping_info.csv")
write.csv(long,"long_format_phyla.csv")
#re-read into R
long<-read.csv("long_format_phyla.csv")
#order according to depth
long$Depth<-factor(long$Depth, levels=c("5", "25", "50", "75","100", "130", "150", "250", "350", "500", "850", "1000", "1500", "2000", "2500", "3000", "4000"))
#set colours
nb.cols<-45
mycolours<-colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)#<- this will need to be manually edited for your orders to be
#plot
long %>%
  ggplot (aes(x=Depth, y=rel_abund, fill=Phylum))+
  geom_bar(stat="identity")+scale_fill_manual(values=mycolours)+theme_bw()+ facet_wrap(Station ~., drop=T, scales="free")+                     # Change font size
  theme(strip.text.x = element_text(size = 14))+theme(axis.text.y= element_text(size=12))+ theme(axis.title=element_text(size=12))+
  theme(axis.text.x= element_text(angle = 45, hjust = 1, size=12)) +ylab(label="Relative abundance")+ xlab(label="Depth (m)")+ theme(legend.position="bottom") 
#check
long2<-subset(long, Station=="23")
long2 %>%
  ggplot (aes(x=Depth, y=rel_abund, fill=Phylum)) +
  geom_bar(stat="identity")+scale_fill_manual(values=mycolours)+theme_bw()+
  theme(axis.text.x= element_text(angle = 45, hjust = 1, size=10)) +ylab(label="Relative abundance")+ xlab(label="Depth (m)")



#-------------------------------------------------------------------------------------------------------------------------------------------------------

#make long
#gather to make long:
long<-gather(long, Site, rel, NSR.practice.QS2:QS.test.Z)
write.csv(long,"long_rel_abund.csv")
#read
long<-read.csv("long_rel_abund_only_thermus.csv")
long2<-read.csv("long_rel_abund.csv")
t<-subset(long2, Family=="Thermaceae")
long<-na.omit(long)
#order
long$Site<-factor(long$Site, levels=c("SGP_1_a", "SGP.1..b.", "LIM.MS..a.", "LIM.MS..b.", "LIM_MS_c", "LIM.North..a.", "LIM.North..b.", "BIM"))
#set colours
nb.cols<-300
mycolours<-colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)#<- this will need to be manually edited for your orders to be
col=c("Thermaceae"="magenta", "Other"="gray35")
long$Added<-factor(long$Added, levels=c("no_additives", "before_extraction", "after_extraction"))
#shared amongst all graphs
#plot
long %>%
  ggplot (aes(x=Type, y=rel, fill=Family)) +
  geom_bar(stat="identity")+scale_fill_manual(values=col)+theme_bw()+ #facet_wrap(~Added)+
  theme(axis.text.x= element_text(angle = 45, hjust = 1, size=10)) +guides(fill=guide_legend(title="Pathway"))+ylab(label="Relative abundance")+ xlab(label="Additives")# theme(legend.position="none") +

#for greenland data
#first amalgamate orders and phyla
asv<-read.csv("Greenland_2018__ASV_table.csv")

#remove first column
asv<-asv[,c(-1,-3)]
#make long
long_v4<-gather(asv, Site, Count, B:L)
#merge all duplicate orders by summing up counts belonging to orders: 
new<-aggregate(Count~Phylum +Site,data=long_v4,FUN=sum)
#turn new into wide format:
wide<-reshape(new, idvar="Phylum", timevar="Site", direction="wide")
#clean up colnames
for ( col in 1:ncol(wide)){
  colnames(wide)[col] <-  sub("Count.", "", colnames(wide)[col])
}
#wide<-wide[,-1]
#write.csv
write.csv(wide, "Amalgamated_Phyla_SILVA.csv")
#cal rel abund for top 20
##calc rel abund for top 20####
v6v8<-read.csv("Top20_phyla_SILVA.csv", row.names=1)#saved to v6v8 folder
v6v8<-read.csv("ind_for_rel.csv", row.names=1)
#transpose
v6v8<-t(v6v8)
x<-v6v8/rowSums(v6v8)
x<-na.omit(x)
tx<-t(x)
write.csv(x, "relative_abundance_for_ind.csv") #in v6v8 folder
#plot
phyla<-read.csv("relative_abundance_Top_20_phyla_SILVA.csv")
order<-read.csv("relative_abundance_Top_20_orders_SILVA.csv")
#gather to make long:
long<-gather(phyla, Phyla, rel, Proteobacteria:Myxococcota)
long$Phyla<-factor(long$Phyla, levels=c("Proteobacteria", "Bacteroidota", "Actinobacteriota", "Acidobacteriota","Cyanobacteria", "Verrucomicrobiota", "Planctomycetota", "Firmicutes", "Crenarchaeota",
                                        "Chloroflexi", "Gemmatimonadota", "Bdellovibrionota", "Armatimonadota", "WPS.2", "Nitrospirota", "Elusimicrobiota","Desulfobacterota", "Patescibacteria", "Unknown", "Myxococcota"))
long$Type<-factor(long$Type, levels=c("subglacial", "ice_marginal_stream", "ice_marginal_lake"))

#for orders
long<-gather(order, Order, rel, Burkholderiales:Propionibacteriales)
long$Order<-factor(long$Order, levels=c("Burkholderiales", "Flavobacteriales", "Cytophagales", "Chitinophagales","Frankiales", "Sphingobacteriales", "Chloroplast", "Rhodobacterales", "Vicinamibacterales",
                                        "Bacteroidales", "Nitrosopumilales", "Acetobacterales", "Xanthomonadales", "Acidobacteriales", "Microtrichales", "uncultured","Gemmatales", "Oceanospirillales", "Gemmatimonadales", "Propionibacteriales"))
long$Type<-factor(long$Type, levels=c("subglacial", "ice_marginal_stream", "ice_marginal_lake"))
#set colours
nb.cols<-20
mycolours<-colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)#<- this will need to be manually edited for your orders to be
#for orders
mycolours<-colorRampPalette(brewer.pal(9, "Set2"))(nb.cols)
#shared amongst all graphs
#plot
long %>%
  ggplot (aes(x=ID, y=rel, fill=Order)) +
  geom_bar(stat="identity")+ scale_fill_manual(values=mycolours)+theme_bw()+facet_grid(. ~ Type, drop=TRUE,scale="free",space="free_x")+
  theme(axis.text.x= element_text(angle = 45, hjust = 1, size=7)) +guides(fill=guide_legend(title="Order"))+ylab(label="Relative abundance")+ xlab(label="Site ID")# theme(legend.position="none") +
