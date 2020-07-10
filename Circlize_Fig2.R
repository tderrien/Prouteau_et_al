## R program for circular representation of SVs
# May 2019
#		* add log2 ratio from dnacopy
# Aim 	: Visualize Strucutral Variants (from Delly)
#
# need 3 input files from the same individuals:
#		 * delly2 .vcf somatic (which ones p a m) => see Stephanie emails
#		 * segments files from Benoit computed using dnacopy?
#		 * numreads to filter SVs


###########################################################################################
library(circlize)
library(tools) # file_path_sans_ext(basename(filepath))
library(vcfR) # dedicated R package for vcf manipulation : https://knausb.github.io/vcfR_documentation/vcfR_object.html
library(dplyr)
library(tidyr)
library(tidyverse)



###########################################################################################
## import the whole variant file using vcfR
###########################################################################################
args <- commandArgs(TRUE)

cat ("* Importing data...\n");

vcffile  		<-      args[1]   # vcf file from Delly 
segfile  		<-      args[2]   # segmentation file from dnacopy
numreads  		<-      args[3]   # number of reads to SR + PR validating the SVs

vcf 	<- read.vcfR( file=vcffile , verbose = FALSE )
dat 	<- read.table(file=segfile, h=T)


cat ("* Transform into tidyR object...\n");
## Transform into tidyR object
vcft 	<- vcfR2tidy(vcf)

cat ("* Parse vcf to filter SVs...\n");
minRead 		<- numreads
minRead_bnd		<- numreads

svtable	<- vcft$fix %>% filter (FILTER == "PASS" & 
                                IMPRECISE== FALSE & 
                                CHROM %in% c(seq(1:38), "X") 
                                & (SR +PE > minRead) ) %>% select (CHROM,POS, END, CHR2,  POS2, SVTYPE, SR, PE) %>% as.data.frame()

# add chr to CHROM and CHR2 fields to be used by circlize which uses chr (e.g UCSC) convention
svtable$CHROM	<- paste('chr', svtable$CHROM, sep='')
svtable$CHR2	<- paste('chr', svtable$CHR2 , sep='')


####################################
#### Get log2ratio data from dnacopy
####################################
cat ("* Parse Segment files ...\n");

dat.clean = dat %>% select (chrom, loc.start, loc.end  , seg.mean ) %>% rename (chr=chrom, start = loc.start, end = loc.end , value= seg.mean)
## add chr
dat.clean$chr	<- paste('chr', dat.clean$chr, sep='')


###########################################################################################
## END Parsing => Draw circles :-)
###########################################################################################
cat ("* Draw SVs...\n");


pdf(file=paste(file_path_sans_ext(basename(vcffile)), minRead,"_bnd",minRead_bnd,"_INV_DUP.pdf", sep=""), height=8, width=8, compress=TRUE)

circos.clear()

## draw canFam3 reference ideograms:
## Warnings : UCSC definition is canFam3 (case sensitive)
circos.initializeWithIdeogram(species = "canFam3")

## Draw log2ratio
circos.genomicTrackPlotRegion(dat.clean,
	panel.fun = function(region, value, ...) { 
		circos.genomicRect(region, value, ytop.column = 1, ybottom = 0,
	col = ifelse(value[[1]] > 0, "red", "green"), ...)
	cell.xlim = get.cell.meta.data("cell.xlim")
	circos.lines(cell.xlim, c(0, 0), lty = 2, col = "#00000040")
})

## Draw SVs
typeE=c("DEL","DUP","INS","INV")
colE=c("blue","red","orange","green")
for (i in 1:4) { 
        bed_list=svtable[svtable$SVTYPE==typeE[i],]
        circos.genomicTrackPlotRegion(bed_list, track.height = uh(5, "mm"), stack=TRUE, panel.fun = function(region, value, ...) {
                circos.genomicPoints(region, value, cex = 0.5, pch = 16, col = colE[i], ...)
        })
}

# BND links
bnd			= svtable %>% filter (SVTYPE == "BND") %>% mutate (SRPE = SR + PE) %>% arrange(SRPE) %>%
filter (SRPE > minRead_bnd)

# group color
nb_grp = 3
bnd$SRPE_grp = cut(bnd$SRPE, nb_grp, labels =seq(1:nb_grp))
colfunc<-colorRampPalette(c("grey","black"))

for (i in 1:nrow(bnd)) {

		# get last element of the color => i.e corresponding to the group
		color_link= tail(colfunc(as.numeric(bnd[i,c("SRPE_grp")])), n=1)
		
		# draw link
        circos.link(bnd[i,c("CHROM")],bnd[i,c("POS")],bnd[i,c("CHR2")],bnd[i,c("POS2")], col=color_link)
}

# Legend
title(paste(file_path_sans_ext(basename(vcffile)), " => Delly SVs calls (filter SR&PR=", minRead, " BND=", minRead_bnd," )", sep=""))


dev.off()

circos.clear()
