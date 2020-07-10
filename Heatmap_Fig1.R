# July 2020
# This R program is meant to be used to reproduce Fig1 of the paper.
#####################################################################

# Load library
library(ComplexHeatmap)
library(circlize)
library(scales) 
library(tidyverse)

# import gene signatures for k = 2 : as given by the Supplementary table 3
genes=rbind(
   data.frame(NMFgeneSign = rep("1", length(rownames(datmad[k2l1,])) ), 
                    RLOC_id     = rownames(datmad[k2l1,])
                    ),
data.frame(NMFgeneSign = rep("2", length(rownames(datmad[k2l2,])) ), 
                    RLOC_id     = rownames(datmad[k2l2,])
                    )
)

# Extraction genes from expression matrice 
dat_sign = datmad [match( genes$RLOC_id, rownames(datmad)), ]

# Important diganostic data from Supp table "Table32RNAseq-30-05"
############################
diag = read.table("Table32RNAseq-30-05.txt", h=T, stringsAsFactors=F)
diag = diag %>% mutate(temp=sample_id) %>% separate(temp, c("dog_id", "breed")) # add breed and dog_id column
diag = diag[match( colnames(dat), diag$sample_id), ] #  reorder rows of diag according to colnames (dat)


#######################################################################################
# Annotation
# colorblind patette
# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#000000", "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
col_SV        <- cbPalette[1:2]
names(col_SV) <- levels(as.factor(diag$max_densite_alteration_cfa_filtre5sondes_manuel))

col_MDM2        <- cbPalette[c(6,8)];names(col_MDM2) <- levels(as.factor(diag$MDM2))

col_breed        <- cbPalette[c(5,7,9)]
names(col_breed) <- levels(as.factor(diag$breed))

# HeatmapAnnotation column
column_ha = HeatmapAnnotation(SV_Content      = diag$max_densite_alteration_cfa_filtre5sondes_manuel,
                              MDM2            = as.factor(diag$MDM2),
                              CDK4            = as.factor(diag$CDK4),
                              CFA10_30        = as.factor(diag$CFA10_CFA30),
                              breed           = diag$breed,
                              col = list( SV_Content      = col_SV,MDM2 = col_MDM2, CDK4 = col_MDM2, CFA10_30= col_MDM2,breed = col_breed
                              ),
                              annotation_name_side = "left"
                            )

# Row Annotation
col_row        <- cbPalette[3:4]
names(col_row) <- levels(as.factor(genes$NMFgeneSign))
row_ha    = rowAnnotation( NMFgeneSign = as.factor(genes$NMFgeneSign), col=list(NMFgeneSign = col_row ))

set.seed(123456)

Heatmap (t(scale(t(dat_sign))), 
  top_annotation = column_ha,
  cluster_columns = res2consmap$Colv,       # from clustering samples from NMF
  show_column_dend = T,
  column_split = 2,
  column_names_gp = gpar(fontsize = 9), 
  
  cluster_rows = FALSE,      # need to be reordered
  show_row_dend    = FALSE,
  show_row_names   = FALSE,
heatmap_legend_param = list(title = "NormExpr"),

left_annotation = row_ha
)
