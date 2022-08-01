
#ploting clone distribution change
library(ggplot2)
library(dplyr)
#library(benthos) #for shannon calculation
setwd("/home/chengchen/data/NMO_project/11PD6/")
barcodes = read.delim("barcodes_real_collapsed.txt",header = FALSE, sep = " ")
barcodes1 = arrange(barcodes, desc(V2))
plot(1:100,barcodes1[1:100,]$V2/2298770)
D6_barcode = barcodes1

setwd("/home/chengchen/data/NMO_project/LYXbarcodes/")
barcodes = read.delim("barcodes_real_collapsed.txt",header = FALSE, sep = " ")
barcodes1 = arrange(barcodes, desc(V2))
sumBC = sum(barcodes$V2)
plot(1:100,barcodes1[1:100,]$V2/sumBC)
D0_barcode = barcodes1

setwd("/home/chengchen/data/NMO_project/11PD4/")
barcodes = read.delim("barcodes_real_collapsed.txt",header = FALSE, sep = " ")
barcodes1 = arrange(barcodes, desc(V2))
D4_barcode = barcodes1

#head 100 reads
Infection = D0_barcode$V2[1:100]/sum(D0_barcode$V2)
NMO_D4 = D4_barcode$V2[1:100]/sum(D4_barcode$V2)
NMO_D6 = D6_barcode$V2[1:100]/sum(D6_barcode$V2)
order = 1:100

Barcode = rep(1:100,3)
frequency = c(Infection,NMO_D4,NMO_D6)
length(frequency)
Sample = c(rep("Infection",100),rep("NMO_D4",100),rep("NMO_D6",100))
df1 = data.frame(Barcode, frequency, Sample)
ggplot(df1, aes(x=Barcode,y=frequency*100,color=Sample))+geom_point(size = 2)+theme_classic()+
  theme(legend.position = c(0.7,0.7))+ylab("Relative frequency (%)")+ylim(0,4)+
  theme(axis.text = element_text(size = 38,color = "black"),axis.title=element_text(size = 44, color = "black"))+
  theme(legend.title = element_text(colour="black", size=44),legend.text = element_text(colour="black", size=38),legend.key.size = unit(3,"line"))+
  scale_color_manual(values = c("#365A87", "#E64B35CC", "#00A087CC"))
 
#integration barcode diversity and indel diversity
#Sample = c("NMO_D4","NMO_D6")
#Barcode = c(250,230)
#Barcode = c(3.40,4.42)
#Indel = c(150,500)
#Indel = c(4.52,5.24)
#BC_indel = c(600,5500)
#BC_indel = c(9.37,11.18)
#df = data.frame(Sample,Barcode,Indel,BC_indel)

Diversity = c(3.40,4.52,9.37,4.42,5.24,11.18)
Sample = c(rep("NMO_D4",3),rep("NMO_D6",3))
Type = c("Barcode","Indel","BC_indel","Barcode","Indel","BC_indel")
df = data.frame(Diversity,Sample,Type)
df$Type = factor(df$Type,levels = c("Barcode","Indel","BC_indel"))
ggplot(df, aes(x=Sample,y=Diversity,fill=Type))+geom_bar(stat="identity",position = position_dodge(0.75))+
  scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF"))+theme_classic()+
  theme(axis.text = element_text(size = 38,color = "black"),axis.title=element_text(size = 44, color = "black"))+
  theme(legend.title = element_text(colour="black", size=44),legend.text = element_text(colour="black", size=38),legend.position = "top")+
  scale_y_continuous(expand = c(0,0),limits = c(0,12.5))+ylab("Diversity index")

#
setwd("/home/chengchen/data/NMO_project/11PD6/")
indel = read.delim("b11p_d6_indels.txt",header = FALSE, sep = " ")
MajorIndel = filter(indel, cumsum(V1)<sum(indel$V1)*0.95)
dim(MajorIndel)
#write MajorIndel into a new file, processed with python to extract site1 indels
indels = read.table("majorIndel11PD6_screen.tsv", header = TRUE,sep = "\t")
indels$percent = indels$number/sum(indels$number)
indels = indels[!duplicated(indels$indels),]
indel = head(indels,21) #remove the highest indel
indel$indels = factor(indel$indels,levels = indel$indels)
ggplot(indel[-1,],aes(x=indels,y=percent*100))+geom_bar(stat = "identity",fill="#0072B2")+theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 30,color = "black"),axis.title.x=element_blank())+
  ylab("Relative frequency (%)")+
  theme(axis.text.y = element_text(size = 38,color = "black"), axis.title.y = element_text(size = 42))+
  scale_y_continuous(expand = c(0,0),limits = c(0,10))+
  theme(plot.margin = unit(rep(1,4),'lines'))


###check the barcode distribution in NMOD31 datasets
#check bulk distribution
nmoD30BCs = read.csv("/home/chengchen/data/NMO_project/LYX_NMO/nmoD30_mcherry/nmoD30_intBC_distribution.csv", header = F)
totalBCs= sum(nmoD30BCs$V1)
nmoD30BCs$log10cellsperBC = log10(nmoD30BCs$V1)
ggplot(nmoD30BCs, aes(x=log10cellsperBC))+geom_histogram(fill="#00A087FF")+theme_classic()+
  theme(axis.text = element_text(size = 38,color = "black"),axis.title=element_text(size = 44, color = "black"))+
  theme(legend.title = element_text(colour="black", size=44),legend.text = element_text(colour="black", size=38),legend.key.size = unit(3,"line"))+
  ylab("numberOfBCs")+coord_flip()



#head 100 reads
top100 = nmoD30BCs$V1[1:100]/totalBCs
Barcodes = 1:100
frequency = top100
df1 = data.frame(Barcodes, frequency)
ggplot(df1, aes(x=Barcodes,y=frequency*100))+geom_point(size = 3, colour="#00A087FF")+theme_classic()+
  theme(legend.position = c(0.7,0.7))+ylab("Relative frequency (%)")+ylim(0,3)+
  theme(axis.text = element_text(size = 38,color = "black"),axis.title=element_text(size = 44, color = "black"))+
  theme(legend.title = element_text(colour="black", size=44),legend.text = element_text(colour="black", size=38),legend.key.size = unit(3,"line"))

#check the cell barcode distribution across cells
cellBC = data.frame(table(nmoD30BCs$V1))
cellBC$numcells = cellBC$Var1
cellBC$numBCs = cellBC$Freq
#ggplot(cellBC, aes(x=numcells,y=numBCs))+geom_point(size = 3, colour="#00A087FF")+theme_classic()+
#  theme(axis.text.x = element_text(size = 20,color = "black",hjust = 1,angle = 90),axis.text.y = element_text(size = 20,color = "black"),axis.title=element_text(size = 44, color = "black"))+
#  theme(legend.title = element_text(colour="black", size=44),legend.text = element_text(colour="black", size=38),legend.key.size = unit(3,"line"))

meanCellPerBC = sum(nmoD30BCs$V1)/nrow(nmoD30BCs) #24.2
newcellBC = rbind(head(cellBC,10),tail(cellBC,10))
newcellBC$group = c(rep("topCells",10), rep("topBCs",10))
#select first 10 and the last 10 rows
ggplot(newcellBC, aes(x=numBCs,y=numcells,fill=group))+geom_bar(stat = "identity")+theme_classic()+
  theme(axis.text = element_text(size = 25,color = "black"),axis.title=element_text(size = 44, color = "black"))+
  theme(legend.title = element_text(colour="black", size=44),legend.text = element_text(colour="black", size=38),legend.key.size = unit(3,"line"))
  


#intBC per cell
nmoD30BCs = read.csv("/home/chengchen/data/NMO_project/LYX_NMO/nmoD30_mcherry/nmoD30_intBC_per_cell.csv", header = F)
BCpercell = data.frame(table(nmoD30BCs$V1))
meanBCpercell =  sum(nmoD30BCs$V1)/nrow(nmoD30BCs)
#meanBCpercell=2.6

ggplot(BCpercell, aes(x=Var1,y=Freq))+geom_bar(stat="identity",fill="#E64B35FF")+theme_classic()+
  theme(axis.text = element_text(size = 33,color = "black"),axis.title=element_text(size = 35, color = "black"))+ylim(0,2000)+
  theme(legend.title = element_text(colour="black", size=44),legend.text = element_text(colour="black", size=38),legend.key.size = unit(3,"line"))+
  coord_flip()+geom_vline(xintercept = meanBCpercell,linetype=3,size=1)


###nmoD30 CRISPR informations
crispr = read.csv("nmoD30indelStat.csv")
crispr$events = factor(crispr$events, levels = c("edited","AmpliconInSc","ScInAmplicon"))
ggplot(crispr,aes(x=events,y=number,fill=states))+geom_bar(stat = "identity")+theme_classic()+
  theme(axis.text = element_text(size = 30,color = "black"),axis.title=element_text(size = 35, color = "black"))+
  theme(legend.title = element_text(colour="black", size=60),legend.text = element_text(colour="black", size=38),legend.key.size = unit(3,"line"))+
  scale_y_continuous(expand = c(0,0),limits = c(0,10000))

