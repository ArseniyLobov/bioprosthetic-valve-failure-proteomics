###This is a code to reproduce statistical analysis of proteomics data from Proteolytic degradation is the culprit behind a bioprosthetic heart valve failure

##Openinig_the_data
#We recommend you to set the working directory for easy reproduction  of the code
#setwd("your directory")

#exctracting the data
dat <- data.frame(read.csv("protein_data.csv"))

dat$Accession <- sub("\\|.*", "", dat$Accession)                   # Extract first three characters

dat1 <- dat[,c(3,21:35)]
head(dat1)

str(dat1)
rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
head(dat1)

## Qualititative analysis
#Extraction group-specific proteins
native <- dat1[which(rowMeans(!is.na(dat1[,c(1:5)])) >= 3/5), ]
bio_1 <- dat1[which(rowMeans(!is.na(dat1[,c(6:10)])) > 3/5), ]
bio_2 <- dat1[which(rowMeans(!is.na(dat1[,c(11:15)])) >= 3/5), ]

#Venn diagram
library(VennDiagram)
library(RColorBrewer)
myCol1 <- brewer.pal(3, "Pastel2")

venn.diagram( #This code will draw diagram to working directory
  x = list(rownames(native), rownames(bio_1), rownames(bio_2)), 
  category.names = c("native" , "bio_1" , "bio_2"),
  filename = '#14_venn_diagramm.png',
  resolution = 600,
  fill = myCol1,
  output=TRUE
)

#Creting the list of group-specific proteins
library(gplots)
v.table <- venn(list(rownames(native), rownames(bio_1), rownames(bio_2)))
print(v.table)

native1 <- c(rownames(native), rep(0, 971 - length(rownames(native))))
bio_11 <- c(rownames(bio_1), rep(0, 971 - length(rownames(bio_1))))
bio_21 <- c(rownames(bio_2), rep(0, 971 - length(rownames(bio_2))))

venn_diff_dat <- data.frame(native1,bio_11,bio_21)
write.csv(venn_diff_dat, "venn_diff_dat.csv")


## Quantitative analysis
#Removing rows with a lot of missing values
dat2 <- dat1[which(rowMeans(!is.na(dat1)) >= 0.85), ]
mean(complete.cases(dat2))
colSums(is.na(dat2))

#knn imputation of missng values
library(impute)
tdat <- t(dat2)
dat_knn1 <- impute.knn(tdat, k = 5)
dat_knn <- t(dat_knn1$data)
mean(complete.cases(dat_knn))

#Opening the sample info
library(readxl)
fact <- data.frame(read_excel("sample_info.xlsx"))

rownames(fact) <- fact[,1]
fact <- fact[,-1]
fact$Type <- as.factor(fact$Type)
fact$Type
fact$Type <- relevel(fact$Type, ref = "Native")
fact$Type

fact$typ2 <- as.factor(fact$typ2)
fact$typ2
fact$typ2 <- relevel(fact$typ2, ref = "Native")
fact$typ2


#Normalization and data QC
library(RColorBrewer)
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$Type]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)
colSums(dat_knn)
#log transformation (погуглите зачем это делать). В чем разница с тем, что было? почему так удобно? зачем прибавлять 1 к dat_knn?
dat_log <- log2(dat_knn+1)
head(dat_log)
mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)

#Quantile normalization
library(limma)
dat_norm <- normalizeQuantiles(dat_log)
head(dat_norm)
boxplot(dat_norm, col = cols, main = "Normalized data")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)
mean(complete.cases(dat_norm))
colSums(is.na(dat_norm))


##sPLS-DA clusterization
t_dat1 <- t(dat_norm)

library(mixOmics)
#PLS-DA 
ordination.optimum.splsda <- splsda(t_dat1, fact$Type, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('PLSDA_cell_contr1.tiff', units="in", width=11, height=6, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination of different valve types", size.title = 1.5, legend=TRUE)
#dev.off()
layout(1,1)

#diff. expression - Naiv_vs_bio

X <- model.matrix(~ fact$typ2)
X

fit <- lmFit(dat_norm, design = X, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit <- eBayes(fit)

# Dif_expr_table
topTable(efit, coef = 2)
full_list_efit <- topTable(efit, number = length(dat_norm))
#write.csv(full_list_efit,'Dif_expr_Native_vs_bio_all.csv')
head(full_list_efit)
#Vulcano-plot
library(EnhancedVolcano)

#tiff('Vulcano_bioall_N.tiff', units="in", width=8, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_efit,
                lab = rownames(full_list_efit),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 2,
              # xlim = c(-8, 10),
              #  ylim = c(0, 10),
                title ="Native versus Bio",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()

#diff expression Bio_1_bio2
fact_bio <- fact[c(2, 4:8, 12:15),]
dat_bio <- dat_norm[, c(2, 4:8, 12:15)]

dim(dat_bio)
dim(fact_bio)

fact_bio$Type <- as.factor(as.character(fact_bio$Type))
fact_bio$Type
X_N2 <- model.matrix(~ fact_bio$Type)
X_N2

fit_N2 <- lmFit(dat_bio, design = X_N2, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit_N2 <- eBayes(fit_N2)

# Dif_expr_table
topTable(efit_N2, coef = 2)
full_list_efit_N2 <- topTable(efit_N2, number = length(dat_bio))
#write.csv(full_list_efit_N2,'Dif_expr_Native_vs_bio2.csv')
head(full_list_efit_N2)
#Vulcano-plot

#tiff('Vulcano_bio1_bio_2.tiff', units="in", width=8, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_efit_N2,
                lab = rownames(full_list_efit_N2),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
             #   ylim = c(0, 7.5),
                FCcutoff = 2,
                title ="Bio1 versus Bio_2",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()