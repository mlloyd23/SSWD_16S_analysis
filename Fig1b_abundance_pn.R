library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)

map<-read.table("map.txt", header=TRUE)

Day0<-dplyr::filter(map, map$Day=="0")

Day3<-dplyr::filter(map, map$Day=="3")

Day6<-dplyr::filter(map, map$Day=="6")

Day9<-dplyr::filter(map, map$Day=="9")

Day12<-dplyr::filter(map, map$Day=="12")

Day15<-dplyr::filter(map, map$Day=="15")

days<-c("Day 0", "Day 3", "Day 6", "Day 9", "Day 12", "Day 15")
num0<-c(37, sum(Day3$Pheno_num==0), sum(Day6$Pheno_num==0), 
sum(Day9$Pheno_num==0), sum(Day12$Pheno_num==0), sum(Day15$Pheno_num==0))

num1<-c(0, sum(Day3$Pheno_num==1), sum(Day6$Pheno_num==1), 
sum(Day9$Pheno_num==1), sum(Day12$Pheno_num==1), sum(Day15$Pheno_num==1))

num2<-c(sum(Day0$Pheno_num==2), sum(Day3$Pheno_num==2), sum(Day6$Pheno_num==2), 
sum(Day9$Pheno_num==2), sum(Day12$Pheno_num==2), sum(Day15$Pheno_num==2))

num3<-c(sum(Day0$Pheno_num==3), sum(Day3$Pheno_num==3), sum(Day6$Pheno_num==3), 
sum(Day9$Pheno_num==3), sum(Day12$Pheno_num==3), sum(Day15$Pheno_num==3))

num4<-c(sum(Day0$Pheno_num==4), sum(Day3$Pheno_num==4), sum(Day6$Pheno_num==4), 
sum(Day9$Pheno_num==4), sum(Day12$Pheno_num==4), sum(Day15$Pheno_num==4))

num5<-c(sum(Day0$Pheno_num==5), sum(Day3$Pheno_num==5), (sum(Day6$Pheno_num==5)+3), 
(sum(Day9$Pheno_num==5)+8), (sum(Day12$Pheno_num==5)+16), (sum(Day15$Pheno_num==5)+19))

class(num0)
data<-as.data.frame(t(data.frame(num0, num1, num2, num3, num4, num5)))
colnames(data) <- c("0","3", "6", "9", "12", "15")
data

my_datm <- melt(cbind(data, ind = rownames(data)), id.vars = c('ind'))
my_datm
colnames(my_datm)<-c("Symptom_numbers", "Day", "Value")

theme_set(theme_bw())

getPalette = colorRampPalette(brewer.pal(6, "YlGnBu"))
orderPalette = getPalette(6)
orderPalette

labels=c(num0="0",num1="1",num2="2",num3="3",num4="4",num5="5")

pdf("abundance_by_pheno_num.pdf")
ggplot(my_datm,aes(x = Day, y = Value, fill = Symptom_numbers)) + 
    geom_bar(position = "fill",stat = "identity") + 
    xlab("Day") + ylab("Proportion") + scale_color_manual(values=orderPalette) + 
	scale_fill_manual(values= orderPalette, labels=c("0","1","2","3","4","5"))+ 
	theme(axis.title = element_text(color="#666666", face="bold", size=16)) +
	theme(legend.text=element_text(size=12)) + 
	theme(legend.title = element_text(colour="#666666", size=12, face="bold"))+
	guides(fill=guide_legend(title="Symptom\nNumber"))
dev.off()




