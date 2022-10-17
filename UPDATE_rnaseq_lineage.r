### RERUN TRANSCRIPTOMIC PSEUDOTIME ANALYSIS WITH UPDATED AMP-AD 2.0 DATA. USE RNASEQ HARMONIZATION PROJECT COUNTS AND DE GENES

#function to run monocle analysis
RunMonocleTobit <- function(Dat, Labels, max_components=2, meth = 'DDRTree',C_by = NULL, 
                            gene_short_name = NULL){ 
  
  library(monocle)
  
  HSMM_expr_matrix <- Dat
  names(HSMM_expr_matrix)<-seq(1,dim(Dat)[2])
  
  if(is.null(gene_short_name)){
    gene_short_name <- c(1:dim(Dat)[1])
  }
  
  
  gene_short_name <- data.frame(gene_short_name)
  Labels <- data.frame(Labels)
  rownames(Labels) <- seq(1,dim(Dat)[2])
  
  pd <- new("AnnotatedDataFrame", data = Labels)
  fd <- new("AnnotatedDataFrame", data = gene_short_name)
  
  
  
  HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                         phenoData = pd,
                         featureData = fd,
                         expressionFamily=tobit())
  
  #HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth, residualModelFormulaStr = ~pmi+educ)
  HSMM <- reduceDimension(HSMM, max_components=max_components, reduction_method = meth)
  #HSMM <- orderCells(HSMM, reverse=TRUE)
  HSMM <- orderCells(HSMM)
  if(is.null(C_by)){
    plot_cell_trajectory(HSMM, color_by="Labels")
  }
  else{
    plot_cell_trajectory(HSMM, color_by=C_by)
  }
  
  
  return(HSMM)
  
}







#try residualized counts matrix (diagnosis-sex)
dlpfcCPMObj <- synapser::synGet('syn26967455')
#Dat <- read.delim(dlpfcCPMObj$path)
Dat <- readr::read_tsv(dlpfcCPMObj$path)

#load rosmap clinical metadata
dlpfcCovObj <- synapser::synGet('syn25808375')
metadata <- read.csv(dlpfcCovObj$path,stringsAsFactors = F)
#limit to DLPFC samples
metadata <- dplyr::filter(metadata, tissue=='DLPFC')

#load differential expression results (AD case vs control, by sex)
de_file <- synapser::synGet('syn26967458')
de1 <- read.delim(de_file$path)
de_male <- dplyr::filter(de1,Comparison=='AD_male_DLPFC - CT_male_DLPFC')
#make short gene names unique
de_male$hgnc_symbol<-make.unique(de_male$hgnc_symbol)
de_female <- dplyr::filter(de1,Comparison=='AD_female_DLPFC - CT_female_DLPFC')
de_female$hgnc_symbol<-make.unique(de_female$hgnc_symbol)

#match short gene names in de file to ensembl features in Dat and use as row names in Dat
genekey <- subset(de_female, select=c(ensembl_gene_id, hgnc_symbol))
names(genekey)[names(genekey) == 'ensembl_gene_id'] <- 'feature'
Dat2 <- dplyr::left_join(Dat, genekey)
Dat2 <- Dat2 %>%
  relocate(hgnc_symbol)
Dat2 <- as.data.frame(Dat2)
rownames(Dat2) <- Dat2$hgnc_symbol

#subset genes by sex
InM <- which(de_male$adj.P.Val<0.1)
MaleGenes <- de_male[InM,]
InF <- which(de_female$adj.P.Val<0.1)
FemaleGenes <- de_female[InF,]


#Limit rnaseq matrix by DE genes, by sex
GeneNames <- Dat2$hgnc_symbol
GeneNamesMale <- MaleGenes$hgnc_symbol
GeneNamesFemale <- FemaleGenes$hgnc_symbol

Names <- colnames(Dat2)
cNames <- metadata$specimenID
l <- length(Names)

#deleting columns not in the covariate list
temp <- rep(T,l)
for (i in 1:l){
  if (!(Names[i] %in% cNames)){
    temp[i] <- F
  }
}

In <- which(temp)
#print(temp)
Dat3 <- Dat2[,In]

#deleting extra rows in covariate list
Names <- Names[In]
l <- length(cNames)
temp <- rep(T,l)
for (i in 1:l){
  if (!(cNames[i] %in% Names)){
    temp[i] <- F
  }
}
In <- which(temp)
metadata <- metadata[In,]

# ColNorm <- function(Dat3){
#   
#   M = max(colSums(Dat3))
#   l <- length(colnames(Dat3))
#   
#   for( i in 1:l){
#     
#     Dat3[,i] = Dat3[,i]*(M/sum(Dat3[,i]))
#     
#   }
#   
#   return(Dat3)
# }

#DatNorm <- ColNorm(Dat3)
DatNorm <- Dat3


#male dataset
In_S <- which(metadata$sex == 'male')
DatNorm_male <- DatNorm[,In_S]
#save matrix and metadata for DE statistics:
saveRDS(DatNorm_male, file="~/prot-lineage/data_objects/Male_fulldatamatrix.RDS")
metadata_male <- metadata[In_S,]
In_genes <- which(GeneNames %in% GeneNamesMale)
DatNorm_male <- DatNorm_male[In_genes,]
GeneNamesMale <- GeneNames[In_genes]



#female dataset
In_S <- which(metadata$sex == 'female')
DatNorm_female <- DatNorm[,In_S]
saveRDS(DatNorm_female, file="~/prot-lineage/data_objects/Female_fulldatamatrix.RDS")
metadata_female <- metadata[In_S,]
In_genes <- which(GeneNames %in% GeneNamesFemale)
DatNorm_female <- DatNorm_female[In_genes,]
GeneNamesFemale <- GeneNames[In_genes]




#prepare for monocle (rerun for each sex)
#males
temp <- DatNorm_male
temp2 <- metadata_male
rownames(temp)<-NULL
rownames(temp2)<-NULL
gene_short_name <- rownames(DatNorm_male)
#save gene list for later use:
rnaseq_genesM <- as.data.frame(gene_short_name)
write.csv(rnaseq_genesM, file="~/prot-lineage/data_objects/rnaseq_genesM.csv", row.names=FALSE)


#females
temp <- DatNorm_female
temp2 <- metadata_female
rownames(temp)<-NULL
rownames(temp2)<-NULL
gene_short_name <- rownames(DatNorm_female)
#save gene list for later use:
rnaseq_genesF <- as.data.frame(gene_short_name)
write.csv(rnaseq_genesF, file="~/prot-lineage/data_objects/rnaseq_genesF.csv", row.names=FALSE)



#Run Monocle2: (ignore warning messages that occur)
#male rnaseq dataset needs to be run with order Reversed (HSMM <- orderCells(HSMM, reverse=TRUE)
MonRun <- RunMonocleTobit(temp, temp2, C_by = 'Pseudotime',gene_short_name = gene_short_name)
g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g

g<- plot_cell_trajectory(MonRun,color_by = "State",show_branch_points=F,use_color_gradient = F,cell_size = 1)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="diagnosis")
g


plot_cell_trajectory(MonRun,color_by = "final_batch",show_branch_points=F,use_color_gradient = F,cell_size = 1)
plot_cell_trajectory(MonRun,color_by = "RIN",show_branch_points=F,use_color_gradient = F,cell_size = 1)
plot_cell_trajectory(MonRun,color_by = "pmi",show_branch_points=F,use_color_gradient = F,cell_size = 1)
table(MonRun$State)


#male samples
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 3] <- 2
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 5] <- 3
MonRun$State2[MonRun$State == 6] <- 5
MonRun$State2[MonRun$State == 7] <- 4

#female samples
MonRun$State2 <- MonRun$State
MonRun$State2[MonRun$State == 7] <- 2
MonRun$State2[MonRun$State == 2] <- 3
MonRun$State2[MonRun$State == 4] <- 3
MonRun$State2[MonRun$State == 3] <- 4
MonRun$State2[MonRun$State == 6] <- 5
MonRun$State2[MonRun$State == 5] <- 6

MonRun$State2 <- as.numeric(MonRun$State2)
MonRun$State2 <- as.factor(MonRun$State2)
table(MonRun$State2)

#save Monocle object for later
saveRDS(MonRun, file='~/prot-lineage/data_objects/MonRun_female.RDS')
saveRDS(MonRun, file='~/prot-lineage/data_objects/MonRun_male.RDS')

tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/rnaseq_figures/MALE_tree_state.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "State2",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="State")
g
dev.off()

tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/rnaseq_figures/MALE_tree_diagnosis.tiff',height=85,width=100,units='mm',res=300)
g<- plot_cell_trajectory(MonRun,color_by = "diagnosis",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
g <- g + ggplot2::scale_color_viridis_d()
g <- g + ggplot2::labs(color="Diagnosis")
g
dev.off()

# #tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
# tiff(file='~/prot-lineage/rnaseq_figures/MALE_tree_braak.tiff',height=85,width=100,units='mm',res=300)
# g<- plot_cell_trajectory(MonRun,color_by = "braaksc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
# g <- g + ggplot2::scale_color_viridis_d()
# g <- g + ggplot2::labs(color="Braak Score")
# g
# dev.off()
# 
# #MonRun$APOE <- factor(MonRun$APOE,levels=c(0,1,2))
# MonRun$apoe_genotype <- as.character(MonRun$apoe_genotype)
# #tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_tree_apoe.tiff',height=85,width=100,units='mm',res=300)
# tiff(file='~/prot-lineage/rnaseq_figures/figures/MALE_tree_apoe.tiff',height=85,width=100,units='mm',res=300)
# g<- plot_cell_trajectory(MonRun,color_by = "apoe_genotype",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
# g <- g + ggplot2::scale_color_viridis_d()
# g <- g + ggplot2::labs(color="APOE genotype")
# g
# dev.off()
# 
# 
# #tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
# tiff(file='~/prot-lineage/rnaseq_figures/figures/MALE_tree_cerad.tiff',height=85,width=100,units='mm',res=300)
# g<- plot_cell_trajectory(MonRun,color_by = "ceradsc",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
# g <- g + ggplot2::scale_color_viridis_d()
# g <- g + ggplot2::labs(color="CERAD Score")
# g
# dev.off()
# 
# #tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_tree_cogdx.tiff',height=85,width=100,units='mm',res=300)
# tiff(file='~/prot-lineage/rnaseq_figures/figures/MALE_tree_cogdx.tiff',height=85,width=100,units='mm',res=300)
# g<- plot_cell_trajectory(MonRun,color_by = "cogdx",show_branch_points=F,use_color_gradient = F,cell_size = 0.5)
# g <- g + ggplot2::scale_color_viridis_d()
# g <- g + ggplot2::labs(color="Cognitive Diagnosis")
# g
# dev.off()

MonRun$braaksc <- as.factor(MonRun$braaksc)
tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/rnaseq_figures/MALE_bargraph_braak.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=braaksc, y=scale(Pseudotime,center=F),fill=braaksc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Braak\nScore",y="Pseudotime",x="Braak Score")
g
dev.off()

MonRun$ceradsc <- as.factor(MonRun$ceradsc)
MonRun$ceradsc <- fct_rev(MonRun$cerads)
tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/rnaseq_figures/MALE_bargraph_cerad.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=ceradsc, y=scale(Pseudotime,center=F),fill=ceradsc)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="CERAD\nScore",y="Pseudotime",x="CERAD Score")
g
dev.off()

MonRun$cogdx <- as.factor(MonRun$cogdx)
tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
#tiff(file='~/prot-lineage/rnaseq_figures/MALE_bargraph_cogdx.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(MonRun@phenoData@data, aes(x=cogdx, y=scale(Pseudotime,center=F),fill=cogdx)) 
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::stat_summary(fun.y=mean, geom="point", shape=23, size=2)
g <- g + ggplot2::theme(axis.text=element_text(size=15), axis.title=element_text(size=15,face="bold"),
                        legend.text=element_text(size=15)) 
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(fill="Cognitive\nDiagnosis",y="Pseudotime",x="Cognitive Diagnosis")
g
dev.off()





######## for stats & other figures, create a dataframe with all relevant covariates & pseudotime & state

x <- list()
x$msex <- MonRun$sex
x$SampleID <- MonRun$specimenID
x$individualID <- MonRun$individualID
x$tissue <- MonRun$tissue
x$State <- MonRun$State2
x$Pseudotime <- MonRun$Pseudotime
x$diagnosis <- MonRun$diagnosis
x$braaksc <- MonRun$braaksc
x$ceradsc <- MonRun$ceradsc
x$cogdx <- MonRun$cogdx
x$apoe4_allele <- MonRun$apoe4_allele
x$pmi <- MonRun$pmi
x$RIN <- MonRun$RIN
x$batch <- MonRun$final_batch
#x$mmse <- MonRun$cts_mmse30_lv
#x$age_death <- MonRun$age_death
x$SampleID <- as.character(x$SampleID)

#rename and create a scaled pseudotime variable
pseudo <- as.data.frame(x)
pseudo$pseudotime_sc <- scale(pseudo$Pseudotime, center=F)

#save variables file for later
write.csv(pseudo, file="~/prot-lineage/data_objects/female_pseudotimes_states.csv", row.names=FALSE)
file <- synapser::File(path='~/prot-lineage/data_objects/female_pseudotimes_states.csv', parentId='syn38349639')
file <- synapser::synStore(file)

#write.csv(pseudo, file="~/prot-lineage/data_objects/male_pseudotimes_state.csv", row.names=FALSE)
#file <- synapser::File(path='~/prot-lineage/data_objects/male_pseudotimes_state.csv', parentId='syn38349639')
#file <- synapser::synStore(file)


#run logistic regression comparing pseudotime between cases and controls only
casecontrol <- subset(pseudo, pseudo$diagnosis=='AD'|pseudo$diagnosis=='CT')
casecontrol$diag2 <- ifelse(casecontrol$diagnosis=='AD', 1, 0)

summary(glm(diag2 ~ pseudotime_sc,casecontrol,family='binomial'))

tiff(file='~/prot-lineage/rnaseq_figures/FEMALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)
#tiff(file='~/prot-lineage/rnaseq_figures/MALE_bargraph_diagnosis.tiff',height=170,width=200,units='mm',res=300)

g <- ggplot(casecontrol,aes(x=diagnosis,
                             y=pseudotime_sc,
                             color=diagnosis)) + geom_boxplot()
g <- g + ggplot2::geom_boxplot() + theme(text = element_text(size = 22)) + theme(legend.position = "none")
g <- g + ggplot2::geom_point(size=2.5, position=ggplot2::position_jitterdodge())
g <- g + ggplot2::scale_color_manual(values=viridis::viridis(3)[1:2])
g
dev.off()

#run proportional odds logistic regression for neuropath/cognitive endpoints:
braakfit <- MASS::polr(braaksc ~ pseudotime_sc,pseudo)
ceradfit <- MASS::polr(ceradsc ~ pseudotime_sc,pseudo)
cogdxfit <- MASS::polr(cogdx ~ pseudotime_sc,pseudo)

cat('braak p-value: ',pt(abs(summary(braakfit)$coef[1,3]),braakfit$df.residual,lower.tail=F)*2,'\n')
cat('cerad p-value: ',pt(abs(summary(ceradfit)$coef[1,3]),ceradfit$df.residual,lower.tail=F)*2,'\n')
cat('cogdx p-value: ',pt(abs(summary(cogdxfit)$coef[1,3]),cogdxfit$df.residual,lower.tail=F)*2,'\n')



#look for correlations with GWAS LOAD genes
ad_gwas <- c("CR1",
             "BIN1",
             "INPP5D",
             "HLA-DRB1",
             "TREM2",
             "MEF2C",
             "NME8",
             "CD2AP",
             "NYAP1",
             "EPHA1",
             "PTK2B",
             "CLU",
             "SPI1",
             "MS4A2",
             "PICALM",
             "SORL1",
             "FERMT2",
             "SLC24A4",
             "ABCA7",
             "APOE",
             "CASS4",
             "ECHDC3",
             "ACE",
             "NDUFAF6",
             "ECHDC3",
             "ADAMTS20",
             "SPPL2A",
             "ADAM10",
             "IQCK",
             "MIR142",
             "ACE",
             "ADAMTS1",
             "SUCLG2P4",
             "FST",
             "OARD1",
             "WWOX",
             "MAF",
             "CD55",
             "YOD1",
             "HLA-DRB1",
             "PSMB8",
             "C4A",
             "GPSM3",
             "HLA-DPA1",
             "HLA-DQA1",
             "HLA-DRA",
             "HLA-DRB5",
             "PSMB9",
             "CD2AP",
             "AGFG2",
             "PILRA",
             "EPHB4",
             "C7orf43",
             "GAL3ST4",
             "ZKSCCAN1",
             "FAM131B",
             "PSMC3",
             "ACP2",
             "C1QTNF4",
             "CELF1",
             "MTCH2",
             "NDUFS3",
             "NUP160",
             "MS4A6A",
             "MS4A7",
             "MS4A4A",
             "EED",
             "PICALM",
             "STYX",
             "RIN3",
             "HMHA1",
             "CNN2",
             "WDR18",
             "CASS4")

Dat3 <- Dat2
Dat3$gene_names <- rownames(Dat3)
Dat3$gene_short_name <- gsub("\\|.*", "", Dat3$gene_names)


# dlpfcCPMObj <- synapser::synGet('syn8456638')
# Dat <- data.table::fread(dlpfcCPMObj$path,data.table=F)
sampleIds <- colnames(Dat2)#[-120]
#sampleIds <- sampleIds[-120]
sampleIds
geneIds <- Dat3$gene_short_name
Dat3$gene_short_name<-NULL
Dat3$gene_names<-NULL
Dat3 <- t(Dat3)
colnames(Dat3) <- geneIds
Dat3 <- data.frame(Dat3,stringsAsFactors=F)
Dat3$sampleId <- sampleIds
dlpfc <- dplyr::left_join(Fvariables,Dat3,by=c('SampleID'='sampleId'))
dlpfc2 <- dlpfc[,14:1850]

corvec <- cor(dlpfc2,dlpfc$Pseudotime,method='spearman')
corDfdlpfc <- data.frame(geneid=colnames(dlpfc2),cor=corvec,stringsAsFactors=F)
corDfdlpfc2 <- corDfdlpfc
corDfdlpfc2$cor <- NULL
corDfdlpfc2 <- dplyr::left_join(corDfdlpfc,corDfdlpfc2,by=c('geneid'))

corDfdlpfc <- corDfdlpfc2
mean(abs(corDfdlpfc$cor))
mean(abs(corDfdlpfc[corDfdlpfc$geneid %in% ad_gwas,]$cor))

mean(corDfdlpfc$cor)
mean(corDfdlpfc[corDfdlpfc$geneid %in% ad_gwas,]$cor)


corDfdlpfc$adGwas <- corDfdlpfc$geneid %in% ad_gwas
colnames(corDfdlpfc)[3] <- 'LOADGWASGene'
corDfdlpfc$LOADGWASGene2 <- ifelse(corDfdlpfc$LOADGWASGene==FALSE, "NOT GWAS GENE", "GWAS GENE")

#tiff(file='~/prot-lineage/figures/FEMALE_loadgwas_cor.tiff',height=85,width=100,units='mm',res=300)
tiff(file='~/prot-lineage/figures/MALE_loadgwas_cor.tiff',height=85,width=100,units='mm',res=300)
g <- ggplot2::ggplot(corDfdlpfc,ggplot2::aes(x=LOADGWASGene2,y=cor,fill=LOADGWASGene2))
g <- g + ggplot2::geom_boxplot() + theme(legend.position="none")
g <- g + ggplot2::scale_fill_viridis_d()
g <- g + ggplot2::labs(x = 'LOAD GWAS STATUS',y='Correlation with pseudotime')
g
dev.off()


