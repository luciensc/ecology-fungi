# author: Lucien Schläpfer
# date: 01.11.2017 / updated 22.08.2020

library(ggplot2)
library(reshape2)
library(dplyr) # a.o. for select fxn

#setwd("/Users/lucien_schlaepfer/PycharmProjects/BK2_data_analysis") # TODO: SET RELATIVE PATH??

dat <- read.table("pilze.csv", na.strings="NA", sep = ";", header = T, stringsAsFactors=T)


### create summary columns
# combine by forest type
dat$F <- dat$exk1F+dat$exk2F
dat$B <- dat$exk1B+dat$exk2B
# combine by day (~1 week delta t)
dat$d1 <- dat$exk1F+dat$exk1B
dat$d2 <- dat$exk2F+dat$exk2B

### specify level order of the categorical data features -> right order for plotting
#str(dat)
x <- levels(dat$einteilung)
x <- x[c(3,5,6,2,7,4,1)] # specify order of taxa (loosely by phylogeny)
dat$einteilung <- ordered(dat$einteilung, x) # order by x
dat$ecol <- ordered(dat$ecol, c("M","S","P/H","P"))
#str(dat)

freq.bin <- as.data.frame(+(select(dat, "F", "B", "d1", "d2") >= 1))

#################################################################
######## DATA ANALYSIS ##########################################
#################################################################

# obtain counts of different properties
dat$standortunabh <- (freq.bin$F==1)&(freq.bin$B==1)
dat$fichteonly <- (freq.bin$F==1)&(freq.bin$B==0)
dat$bucheonly <- (freq.bin$F==0)&(freq.bin$B==1)
dat$tagunabh <- (freq.bin$d1==1)&(freq.bin$d2==1)
dat$fichte <- (freq.bin$F==1)
dat$buche <- (freq.bin$B==1)

l <- length(dat$name)
cat("present at pine site only:", sum(dat$fichteonly)/l)
cat("present at beech site only:", sum(dat$bucheonly)/l)
cat("present at both sites:", sum(dat$standortunabh)/l)
cat("present on both days:", sum(dat$tagunabh)/l)

# barplot abundance by taxon
ggplot(dat, aes(x=einteilung))+
  geom_bar()+
  coord_flip()+
  theme_light()+
  scale_fill_grey()+
  labs(x="",y="Anzahl Arten")+
  stat_count(aes(label=..count..), vjust=+0.5, hjust=-0.5, geom="text", position="identity")
ggsave("barplot1 einteilung.png", width=8,height=5)

# detailed ecological function of species by site
ecol.dat <- select(dat, "ecol", "mbaum", "ssubstrat", "fichte", "buche")
ecol.dat <- ecol.dat[complete.cases(ecol.dat$ecol, ecol.dat$mbaum, ecol.dat$ssubstrat),]  #rm rows with NAs in ecol info
# get ecological function information with relevant detail into one column
ecol.dat$finer <- as.character(ecol.dat$ecol)  # factors: encoded as numeric -> need to switch to char to combine vec.s
ecol.dat$finer[ecol.dat$ecol=="M"] <- as.character(ecol.dat$mbaum[ecol.dat$ecol=="M"])
ecol.dat$finer[ecol.dat$ecol=="S"] <- as.character(ecol.dat$ssubstrat[ecol.dat$ecol=="S"])
ecol.dat$finer <- as.factor(ecol.dat$finer)  # generate levels again (for table building)

# obtain counts of ecological function per site
Fichte <- table(ecol.dat$finer[ecol.dat$fichte])
Buche <- table(ecol.dat$finer[ecol.dat$buche])
ecol.tbl <- as.data.frame(rbind(Fichte, Buche))

# distribute hybrid category counts into pure categories
ecol.tbl$B <- ecol.tbl$B+0.5*ecol.tbl$`B/F`
ecol.tbl$F <- ecol.tbl$F+0.5*ecol.tbl$`B/F`
ecol.tbl$H <- ecol.tbl$H+0.5*ecol.tbl$`P/H`+0.5*ecol.tbl$`S/H`
ecol.tbl$P <- ecol.tbl$P+0.5*ecol.tbl$`P/H`
ecol.tbl$S <- ecol.tbl$S+0.5*ecol.tbl$`S/H`
# remove unneeded categories
ecol.tbl$`B/F` <- NULL
ecol.tbl$`P/H` <- NULL
ecol.tbl$`S/H` <- NULL
rownames(ecol.tbl) <- c("Fichtenwald", "Buchenwald")
colnames(ecol.tbl) <- c("Buche (M.)","Fichte (M.)", "Totholz (S.)", "Div. (P.)", "Streu (S.)")
ecol.tbl <- ecol.tbl[,c(1,2,5,3,4)]
#rowSums(ecol.tbl)  # pine: 41 species, beech: 35 species
ecol.f <- ecol.tbl/rowSums(ecol.tbl)  # normalise
ecol.f$wald <- rownames(ecol.f)
ecol.f.m <- melt(ecol.f, id.vars = "wald")

# generate stacked barplot
ggplot(ecol.f.m, aes(x = wald, y = value, fill = variable, label = round(value,2))) + #label rounded to 2 decimal places
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x="",y="Anteil an Artenzahl")+
  scale_fill_manual("oekol. Assoziation", values = c("palegreen3", "palegreen4", "peachpuff3","peachpuff4","red4"))
ggsave("w5h8 barplot2 function stacked.png", width=5,height=8)

