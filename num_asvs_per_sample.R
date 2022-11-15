library(ggplot2)
library(tidyr)
library(dyplr)
#read in file
asv<-read.csv("Greenland_2018__ASV_table.csv")
#make long
long<-gather(asv, Sample, count, B:L)

col=c("B"="grey", "C"="grey", "G"="lightblue", "H"="lightblue", "I"="lightblue", "J"="tan", "K"="tan", "L"="tan")
#plot
long %>%
  ggplot (aes(x=Sample, y=count, fill=Sample)) +
  geom_bar(stat="identity")+ scale_fill_manual(values=col)+theme_bw()+
  theme(axis.text.x= element_text(angle = 45, hjust = 1, size=7)) +guides(fill=guide_legend(title="Sample"))+ylab(label="Number of unique sequences (ASVs) obtained")+ xlab(label="Sample ID")# theme(legend.position="none") +
##rarefaction curve######
data(BCI) #example data
asv<-read.csv("Greenland_2018__ASV_table.csv", row.names=1)
asv<-t(asv)
S <- specnumber(asv) # observed number of species
(raremax <- min(rowSums(asv))) #rarefy down to lowest number of asvs in a sample,  in this case 7000
Srare <- rarefy(asv, raremax) #rarefies the dataset down to lowest number of asvs
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")# plots total observed species
abline(0, 1)
rarecurve(asv, step = 20, sample = raremax, col = "blue", cex = 0.6)
