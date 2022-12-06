############################################################################################    
#       S1: R codes for phylogenetic analyses by using FASTA format file

#	A workflow with RStudio: phylogenetic analyses and 
#	      visualizations using mtDNA-Cytb sequences
#		
#	      Emine TOPARSLAN, Kemal KARABAG, Ugur BILGE
#			       June 2020

#Device configration
#Processor        : Intel(R) Core(TM) i5-4200 CPU @2.50 GHz 
#Installed memmory: 4,00 GB
#System type      : 64 bit Operating system, x64-based processor

#R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
#Copyright (C) 2020 The R Foundation for Statistical Computing
#Platform: x86_64-w64-mingw32/x64 (64-bit)
###########################################################################33

#This workflow is planned according to data sets without missing data and gaps.
#Before starting, please check your sequences that have any gaps or missing and fix them.


#The following packages must be installed to perform all the analyses outlined in the article. 
#Use these commands if you have not already installed them. 
#If you have already them, you can skip this part (InstallPackages is set to FALSE by default)

#If your samples belonging to populations or groups, 
#you must modify the name of samples using  "_" between population/group name and number.  
#This naming method allows extracting unique names as population names from sample names with the help of a short command. 
#Thus, the name of the population in all the analyzes do not need to be entered again. 

#If you can't change the name of the samples, please run the commands following. 
#After a few steps, you will see the example of name change codes. Please, follow the commands below.
############################################################################################  

#set working directory and input file name below
#setwd("C:/Users/ugur/desktop/ileribiyoinformatik")
setwd("C:/Users/usr/Documents")

fname = "S2_Appendix.fas"#You must enter the name of the data set you want to read.

InstallPackages = FALSE

if (InstallPackages) {
	if (!requireNamespace("BiocManager", quietly=TRUE)) 
		install.packages("BiocManager")
	BiocManager::install("msa")

	install.packages("adegenet")
	install.packages("ape")
	install.packages("ggtree")
	install.packages("ggplot2")
	install.packages("ips")
	install.packages("bios2mds")
	install.packages("haplotypes")
	install.packages("pegas")
	install.packages("phytools")
	install.packages("stats")
	install.packages("treeio")
}

library(adegenet)
library(ape)
library(ggtree)
library(ggplot2)
library(stats)
library(ips)
library(msa)

##########################   WARNING   ##################################

#We would like to warn those who will use the Align command that the alignment process can take 
#only a minute or a few hours depending on the length of your sequence and the number of samples or strong of your devices. 
#For example, 740 samples that length is 1226 bp may take approximately 4 hours for alignment.
#We aligned FASTA formatted file using msa() command and ordered the names in the "nbin" object, 
#when coloring haplotype networks and renaming samples for grouping.
#Note: "S2_Appendix.fas" file have been non-aligned before.

#Before starting the analysis, it is necessary to align the FASTA format sequences. 
#If your sequences are already aligned, please proceed from "READING AND PLOTTING OF ALIGNMENT-NJ-MSAPLOT" section.
#Even so, if you want to alignment again and export your FASTA file, you can follow instructions as below.
#If your sequences are not aligned, please run the commands following.
#If your samples already aligned, you can also set "AlignNeeded = FALSE" below.

AlignNeeded = TRUE # The program reads fasta file and aligns it

if (AlignNeeded) {
	
	file <- readDNAStringSet(fname)#for reading multiple DNA sequences from msa package
	file

	cb<- msa(file) # multiple sequence alignment from msa package   
	cb

##########   CONVERTING ALIGN FILE TO FASTA FILE    ###############

	cv<-msaConvert(cb, type=c("bios2mds::align"))

##########   EXPORTING  ALIGNED FASTA FILE    ###############           

	library(bios2mds) # for exporting fasta file

	export.fasta(cv, outfile = "outfile.fas", ncol(cb), open = "w")
}
 
#After align and export, please check your sequence that has whether gaps or missing in the "outfile.fas".
#If sequences have gaps or are missing, please fix them and call again using "fasta2DNAbin(fname/"outfile.fas")" as below.
#If your sequences have not gaps or are not missing, you can directly contuniue to "nbin<-as.DNAbin(cb)" command as below.

nbin<-as.DNAbin(cb) #read aligned data from cb above
 

##########  READING AND PLOTTING OF ALIGNMENT-NJ-MSAPLOT   ############### 

#If you don't need to align your file or you get repaired "outfile.fas", 
#you can use "readFASTA2DNAbin(fname)" command to call them to RStudio console as below.

#nbin<-fasta2DNAbin(fname/"outfile.fas")# read your previously aligned FASTA file name.


#After align you must use trimEnds() command as below. 
#If you don't need to trim you can follow the next command.

TRIM = FALSE       # Already trimmed sequence is assumed

if (TRIM) {
	nbin<-trimEnds(nbin)#trimming of sequences ends
}

#############    WARNING    ####################### 
#Before starting the analysis, you should be sure that the names of samples like "name_1; name_2, ...".
#If you couldn't change it before, you can fastly change them as below.
#Firstly, you can write names of aligned sequences to see the name list using the "labels(nbin)" command, 
#then you can rename them as below;
#A, B, and C letters represent names of populations or groups.
#rownames(nbin)[1:10]=paste("A_",1:10, sep="")
#rownames(nbin)[11:20]=paste("B_",1:10, sep="")
#rownames(nbin)[21:30]=paste("C_",1:10, sep="")....
#After changing the name of samples, you can follow codes as below.

an<-as.alignment(nbin)  #converting DNAbin to alignment format
nm<-as.matrix(an)       #converting alignment to matrix
nbinmat<-as.matrix(labels(nbin)) #extraction of the sample names
nbin
class(nbin)
dnbin<-dist.dna(nbin, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree<-nj(dnbin)
ggt<-ggtree(tree, cex = 0.8, aes(color=branch.length))+scale_color_continuous(high='lightskyblue1',low='coral4')+geom_tiplab(align=TRUE, size=2)+geom_treescale(y = - 5, color = "coral4", fontsize = 4)

#The shared S2_Appendix.fas dataset is not contained missing value, different letter (N, R, K etc.) or gaps.
#If your sequences have gaps (-) or N, Y, K, R letter (IUPAC nucleotide code), 
#you have to add additional color in "color argument" below.
#For example, "rep("green",1,)" for gaps, "rep("pink",1,)" for N letter etc.
#However, you have to delete or fix the gaps to be able to do other analysis.
#If you have big data set, you can modified "width" and "height" arguments below for better images.

njmsaplot<-msaplot(ggt, nbin, offset = 0.009, width=1, height = 0.5, color = c(rep("rosybrown", 1), rep("sienna1", 1), rep("lightgoldenrod1", 1), rep("lightskyblue1", 1)))
njmsaplot

dev.new()
njmsaplot


pdf("njmsaplot.pdf", width = 11, height = 9)#save as pdf file
njmsaplot
dev.off()

#########   EXTRACTION SEQUENCE AND HAPLOTYPE INFORMATION    ###############

nrow(nm)#confirmation of number of samples

ncol(nm)#confirmation of sequences size

sat2 <- NULL
for (i in 1:nrow(nm)) {
    sat2[i] <- paste(nm[i, ], collapse="")
}

sat2 <- toupper(sat2) #converts all letters to uppercase
sat3 <- unique(sat2) #gives only unique sequences from all sequences
sat3#that is, it gives complete sequences of haplotypes (20x373).
hfreq <- NULL
for (i in 1:length(sat3)) {
    hcount = 0
    s3 <- sat3[i]
    for (j in 1:length(sat2)) {
        s2 <- sat2[j]
        if (s3 == s2) {
            hcount <- (hcount + 1) #counts the number of individuals with the same haplotype sequence. 
            #print(paste(i, "yes", hcount))
        }
        #print(s2)
    }
    hname<-(paste("H",i, sep =""))
    hfreq[i] <- hcount
    #print(paste(hname, hcount, collapse = ""))
}   #haplotype frequency in the all samples

len <- nchar(sat3[1]) #assume all have same length!!!
cnt <- 1
sat4 = list()
for (j in 1:len) {
    same <- TRUE
    first <- substr(sat3[1], j, j)
    for (i in 2:length(sat3)) {
        ch1 <- substr(sat3[i], j, j)
        if (first != ch1) {
            str <- paste(j, first, ch1)
            print(str)
            same <- FALSE
            break
        }
    }
    if (!same) {
        ss <- NULL
        for (i in 1:length(sat3)) {
            ss <- paste(ss, substr(sat3[i], j, j), sep="")
        }
        sat4[cnt] <- ss
        cnt <- cnt + 1
    }
}#it gives the mutation points and the nucleotide substitutions

len <- nchar(sat3[1]) #assume all have same length!!!
cnt <- 1
sat5 = list() 
for (j in 1:len) { #scan all columnns and if all elements are the same do not copy
    same <- TRUE
    first <- substr(sat3[1], j, j)
    scol <- first
    for (i in 2:length(sat3)) {
        ch1 <- substr(sat3[i], j, j)
        scol <- paste(scol, ch1, sep="")
        if (first != ch1) {
            str <- paste(j, first, ch1)
            #print(str)
            same <- FALSE
            #break
        }
    }
    if (!same) {
        scol <- paste("V_", cnt, " ", scol, sep="")
        ss <- NULL
        for (i in 1:length(sat3)) {
            ss <- paste(ss, substr(sat3[i], j, j), sep="")
        } 
        sat5[cnt] <- ss
        cnt <- cnt + 1
    }
}

sat6 <- as.matrix(sat5)
mat6 = matrix(nrow=nrow(sat6), ncol=nchar(sat6[1]))
for (i in 1:nrow(mat6)) {
    s <- as.vector(strsplit(as.character(sat5[i]), ""))
    for (j in 1:ncol(mat6)) {
        mat6[i, j] <- as.character(s[[1]][j])
    }
}
mat7 <- t(mat6) #sequences of haplotypes and variable sites matrix (20x41)
write.table(mat7,file="mat7.txt", quote=FALSE, sep="\t")
hname<-paste("H", 1:nrow(mat7), sep = "")
rownames(mat7)=hname
write.table(mat7,file="mat7.txt", quote=FALSE, sep="\t") 

str4 <- NULL
str4[1] <- paste(mat7[1, ], collapse="")
for (i in 2:nrow(mat7)) {
    tmp <- NULL
    for (j in 1:ncol(mat7)) {
        chr = "."
        if(mat7[i, j] != mat7[1, j]) chr = mat7[i, j]
        tmp <- paste(tmp, chr, sep="")
    }
    str4[i] <- paste(tmp, collapse="")
}
nchar(str4[1]) #confirmation of number of variable sites
mstr4<-as.matrix(str4)
rownames(mstr4)<-hname
colnames(mstr4)<-paste("sequences length","(", ncol(mat7), "base pairs", ")")
pct<-round((as.matrix(hfreq)*100/colSums(as.matrix(hfreq))), 2)
colnames(pct)<-c("pct")
cmstr4<-as.data.frame(cbind(mstr4, hfreq, pct))
cmstr4
write.table(cmstr4,file="cmstr4.txt", quote=FALSE, sep="\t") 

#############  HAPLOTYPES FREQUENCY  ################
library(haplotypes)

kn<-as.dna(nbin)
kh<-haplotypes::haplotype(kn)

ncb <- as.matrix(labels(nbin))
n2 <- NULL
for (i in 1:nrow(ncb)) {
    n2[i] <- strsplit(ncb[i], '_')[[1]][1] #to get the names of the examples where the name and number are separated by an underscore
}
n2

for (i in 1:nrow(ncb)) {
    n2[i] <- strsplit(n2[i], ' ')[[1]][1] #to get the names of the examples where the name and number are separated by a space
}
n2

hf<-grouping(kh, factors=n2)
hf[["hapvec"]] <- NULL
dhf<-as.data.frame(hf$hapmat) #haplotype frequency matrix per populaion
rownames(dhf)<-paste("H", 1:nrow(mat7), sep = "")
dhf
write.table(dhf,file="dhf.txt", quote=FALSE, sep="\t")

############## DISTANCE MATRIX OF HAPLOTYPES  ############################
library(pegas)
dh<-dist.hamming(mat7) #from pegas package
dh
dhm<-as.matrix(dh)
write.table(dhm,file="dhm.txt", quote=FALSE, sep="\t")

############### HEAT MAP ###########################
sat2 <- NULL
for (i in 1:nrow(nm)) {
    sat2[i] <- paste(nm[i, ], collapse="")
}
sat2

sat2 <- toupper(sat2)
sat3 <- unique(sat2)
comat = matrix(nrow=length(sat3), ncol=length(sat3))
for (i in 1:length(sat3)) { 
    si <- sat3[i]
    for (j in 1:length(sat3)) { 
        sj <- sat3[j]
        difcnt = 0
        s1 = as.vector(strsplit(as.character(si), ""))
        s2 = as.vector(strsplit(as.character(sj), ""))
        for (k in 1:length(s1[[1]])) {
            if (s1[[1]][k] != s2[[1]][k]) {
                difcnt = difcnt + 1
            }
            comat[i, j] = difcnt
            #print(paste(i, " ", j, " ", difcnt))
        }
    }
}
comat	#is Hamming distance matrix
colnames(comat)<-paste("H", 1:nrow(comat), sep = "")
rownames(comat)<-paste("H", 1:nrow(comat), sep = "")
heatmap(comat,scale="none",col=heat.colors(100),keep.dendro=TRUE, symm=TRUE) #stats package

dev.new()
heatmap(comat,scale="none",col=heat.colors(100),keep.dendro=TRUE, symm=TRUE)


pdf("heatmap.pdf", width = 7, height = 7)
heatmap(comat,scale="none",col=heat.colors(100),keep.dendro=TRUE, symm=TRUE)
dev.off()

############## HAPLOTYPE NETWORK-INDIVIDUALS ##########################
library(pegas)
h<-pegas::haplotype(nbin, strict = FALSE, trailingGapsAsN = TRUE)#extracts haplotypes from DNAbin object
hname<-paste("H", 1:nrow(h), sep = "")
rownames(h)= paste(hname)
net<-haploNet(h, d = NULL, getProb = TRUE)#constructs the haplotype network
net
ind.hap<-with(
    utils::stack(setNames(attr(h, "index"), rownames(h))),
    table(hap=ind, individuals=rownames(nbin))
)

par(mar=c(0.01,0.01,0.01,15))
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=15, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.52, ncol=6, x.intersp=0.2, text.width=11)

dev.new()
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=30, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.5, ncol=6, x.intersp=0,003, text.width=9)

pdf("hapind.pdf", width = 11, height = 5)
par(mar=c(0.01,0.01,0.01,15))
plot(net, size=attr(net, "freq"), scale.ratio = 2, cex = 0.6, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x= 57,y=15, colnames(ind.hap), fill=rainbow(ncol(ind.hap)), cex=0.52, ncol=6, x.intersp=0.2, text.width=11)
dev.off()

#For a better position of the name table, you can optimize x and y coordinates in the "legend()" argument.

############## HAPLOTYPE NETWORK-POPULATIONS ##########################

h<-pegas::haplotype(nbin, strict = FALSE, trailingGapsAsN = TRUE)#extracts haplotypes from DNAbin object
hname<-paste("H", 1:nrow(h), sep = "")
rownames(h)= paste(hname)
net<-haploNet(h, d = NULL, getProb = TRUE) #constructs the haplotype network
net
ind.hap<-with(
    utils::stack(setNames(attr(h, "index"), rownames(h))),
    table(hap=ind, individuals=rownames(nbin)[values])
)

#We followed the name order in the ind.hap object for all coloring of samples in the haplotype network circles.
#You can see the order of sample names using the command "colnames(ind.hap)".
#The color of the samples has been assigned as the "bg" list for use in the plot () command, as below.
#The number of colors is determined by the number of samples in the same population (see "colnames (ind.hap)").
#Using the colors() command, you can choose different colors for your samples according to the populations.

bg<-c(rep("dodgerblue4", 15), rep("olivedrab4",15), rep("royalblue2", 15), rep("red",15), rep("olivedrab3",15), 
 rep("skyblue1", 15), rep("olivedrab1", 15),  rep("darkseagreen1", 15))

#If you are not going to use a special order for the population names to be listed in the legend command (population names listed in the legend ("hapcol"), 
#you can automatically select and sort the population names using the following command ("sn2"). 
#This order is fully consistent with the order in the ind.hap object.

#un2<-unique(n2)
#sn2<-sort(un2)

#Then you must the use "sn2" list instead of "hapcol" list in the legend argument.

#But, we wanted to order names of the population names according to groups. 
#Therefore, we reordered the population names as "hapcol" list.
#So, you can see color transitions of the populations by the groups. 
#The population names were assigned using "hapcol" list in the legend () command.

hapcol<-c("Aksu", "Demre", "Kumluca", "Firm", "Bayatbadem", "Geyikbayir", "Phaselis", "Termessos")

#At this point, the colors to be used in the list where the populations will be shown (legend) must match the colors you choose for the samples in the circles.
#Therefore, if you used the "sn2" list, you can directly use the "bg" list for the populations' color symbols. 
#Note, now the number of colors must be changed to "1" for each population.
#For example;
#bgp<-c(rep("dodgerblue4", 1), rep("olivedrab4",1),...)

#Then you must the use "bgp" list instead of "ubg" list in the legend argument (fill).

#But, we used different order for the population names (hapcol list).
#Thus, colors of the population were assigned as "ubg" list and used in the fill argument of the legend () command.
#The colors of the samples assigned in the bg list and the population colors (ubg) belonging samples must match each other.

ubg<-c(rep("dodgerblue4",1), rep("royalblue2",1), rep("skyblue1",1), rep("red",1), rep("olivedrab4",1), 
rep("olivedrab3",1), rep("olivedrab1",1), rep("darkseagreen1",1))

par(mar=c(0.001,0.001,0.001,0.001))
plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x=-36,y=53, hapcol, fill=ubg, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)

########## USING "sn2" NAME LIST AND "bgp" COLOR LIST  ###############

#par(mar=c(0.001,0.001,0.001,0.001))
#plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
#legend(x=-36,y=53, sn2, fill=bgp, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)

dev.new()
par(mar=c(0.001,0.001,0.001,0.001))
plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x=-36,y=53, hapcol, fill=ubg, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)

pdf("happop.pdf", width = 7, height = 4)
par(mar=c(0.001,0.001,0.001,0.001))
plot(net, size=attr(net, "freq"), bg = bg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hap, show.mutation=1, font=2, fast=TRUE)
legend(x=-36,y=53, hapcol, fill=ubg, cex=0.8, ncol=1, bty="n", x.intersp = 0.2)
dev.off()

############## HAPLOTYPE NETWORK-GROUPS ##########################

#To construct a haplotype network according to the groups, we assigned the group names of the samples according to the nbin object. 
#At this stage, we followed the name order in the nbin object using "labels(nbin)".

nmname<-c(rep("Nature", 3), rep("Nature", 6), rep("Nature", 6), rep("Nature", 4), rep("Firm", 4), 
rep("Nature", 11), rep("Greenhouse", 3), rep("Nature", 3), rep("Greenhouse", 2), rep("Greenhouse", 12), 
rep("Nature", 6), rep("Firm", 4),  rep("Nature", 7), rep("Greenhouse", 13), rep("Nature", 6), rep("Greenhouse", 15), 
rep("Firm", 3), rep("Nature", 6), rep("Nature", 2), rep("Firm", 4))

ng<-nbin
rownames(ng)<-nmname
hg<-pegas::haplotype(ng, strict = FALSE, labels = hname, trailingGapsAsN = TRUE) #extracts haplotypes from DNAbin object
hg
netg<-haploNet(hg, d = NULL, getProb = TRUE) #constructs the haplotype network
netg
ind.hapg<-with(
    utils::stack(setNames(attr(hg, "index"), rownames(hg))),
    table(hap=ind, individuals=rownames(ng)[values])
)

#Group colors in the legend command which gbg argument below, were assigned according to the order of names in the ind.hapg object.
#To see name order, you can use the below command;
#as.data.frame(ind.hap)
#Color of samples was assigned using "gbg" argument in the plot () command.

gbg<-c(rep("red"), rep("blue"), rep("green"))

par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
legend(x=-35,y=45, colnames(ind.hapg), fill=c("red","blue", "green"), cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)

dev.new()
par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
legend(x=-35,y=45, colnames(ind.hapg), fill=c("red","blue", "green"), cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)


pdf("hapgrp.pdf", width = 8, height = 5) #save as pdf file
par(mar=c(0.001,0.001,0.001, 0.001))
plot(netg, size=attr(netg, "freq"), bg = gbg, scale.ratio = 2, cex = 0.7, labels=TRUE, pie = ind.hapg, show.mutation=1, font=2, fast=TRUE)
legend(x=-35,y=45, colnames(ind.hapg), fill=c("red","blue", "green"), cex=0.8, ncol=1, bty = "n", x.intersp = 0.2)
dev.off()

############# CIRCULAR TREE-NJ-POPULATIONS #######################

library(pegas)
library(ggtree)
library(ggplot2)
library(phytools)
library(stats)

nname<-as.matrix(sort(labels(nbin)))#to construct population names 
krp<-list(Aksu = nname[1:15,],
Bayatbadem = nname[16:30,],
Demre = nname[31:45,],
Firm = nname[46:60,],
Geyikbayir = nname[61:75,],
Phaselis = nname[76:90,],
Kumluca = nname[91:105,],
Termessos = nname[106:120,])

nbin
class(nbin)
dnbin<-dist.dna(nbin, model = "K80") #computing distance by ape package with K80 model derived by Kimura (1980)
tree<-nj(dnbin) #neighbor-joining tree estimation of Saitou and Nei (1987).

emos<-ggtree(tree, layout = 'circular', branch.length='branch.length', lwd = 0.5)+xlim(-0.1, NA)
groupOTU(emos, krp, 'Populations') + aes(color=Populations)+theme(legend.position="right")+geom_tiplab(names(nbin), cex = 1.7, offset=0.002)+guides(color = guide_legend(override.aes = list(size = 2.5)))+geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9)

dev.new()
emos<-ggtree(tree, layout = 'circular', branch.length='branch.length', lwd = 0.5)+xlim(-0.1, NA)
groupOTU(emos, krp, 'Populations') + aes(color=Populations)+theme(legend.position="right")+geom_tiplab(names(nbin), cex = 1.7, offset=0.002)+guides(color = guide_legend(override.aes = list(size = 2.5)))+geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9)


pdf("njtreepop.pdf", width = 6.5, height = 5)
emos<-ggtree(tree, layout = 'circular', branch.length='branch.length', lwd = 0.5)+xlim(-0.1, NA)
groupOTU(emos, krp, 'Populations') + aes(color=Populations)+theme(legend.position="right")+geom_tiplab(names(nbin), cex = 1.7, offset=0.002)+guides(color = guide_legend(override.aes = list(size = 2.5)))+geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9)
dev.off()

##################### CIRCULAR TREE-NJ-DISTANCE ##############################

nbin
class(nbin)
dnbin<-dist.dna(nbin, model = "K80")#computing distance by ape package with K80 model derived by Kimura (1980)
tree<-nj(dnbin)#neighbor-joining tree estimation of Saitou and Nei (1987).

njdistree<-ggtree(tree,layout = 'circular', branch.length='branch.length', aes(color=branch.length), lwd = 0.5)+xlim(-0.1, NA)+
           geom_tiplab(names(nbin), size = 1.7, offset=0.002)+scale_color_continuous(high='lightskyblue1',low='coral4')+
           geom_treescale(x=-0.1, color = "coral4", fontsize = 3, offset = 9) 

njdistree

dev.new()
njdistree


pdf("njdistree.pdf", width = 6.5, height = 5)
njdistree
dev.off()

############## HAPLOTYPE TREE ##########################

library(treeio)
D <- dist.hamming(mat7)#pegas package
class(D)
htre<-nj(D)#This function performs the neighbor-joining tree estimation of Saitou and Nei (1987).
bp <- boot.phylo(htre, mat7, B=100, function(x) nj(dist.hamming(x)))
bp2 <- data.frame(node=1:Nnode(htre) + Ntip(htre), bootstrap = bp)
htree <- full_join(htre, bp2, by="node")
boothap<-ggtree(htree, size=1, branch.length='branch.length')+geom_tiplab(size=4)+
geom_nodepoint(aes(fill=cut(bootstrap, c(0,50,70,85,100))), shape=21, size=4)+
theme_tree(legend.position=c(0.85, 0.2))+ 
scale_fill_manual(values=c("black", "red", "pink1", "white"), guide='legend', name='Bootstrap Percentage (BP)',breaks=c('(85,100]', '(70,85]', '(50,70]', '(0,50]'), labels=expression(BP>=85, 70<=BP*"<85",50<=BP*"<70", BP<50))

boothap

dev.new()
boothap


pdf("boothap.pdf", width = 11, height = 5)
boothap
dev.off()