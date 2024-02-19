#Compare the DE in the healthy control cells in HSC versus more mature cells
#I'll try with early GMPs first (50 cells in total)
#Load libraries
library(ExperimentHub)
install.packages("Seurat")
install.packages("Matrix") #loaded separately, MATRIX 1.6-1.5 required
install.packages("SeuratObject")
library(Seurat)
library(DESeq2)
library(tidyverse)
library(RColorBrewer)

#Read in the file including CNAs and annotations from reference map
df<-readRDS("/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/240110_Mapping_to_reference/Seurat_obj_CNA_referenceannotations.rds")

#Create a subset with only control cells 
df<-subset(x=df, subset=Patient=='control')
df<-subset(x=df, subset=predicted_CellType %in% c("HSC", "Early GMP"))
#Check whether it worked
table(df$Patient, df$predicted_CellType) #shows 356 HSC and 50 early GMPs
table(df$Sample_well)
#Sample_well is not a unique identifier

df$samples <- paste0(df$predicted_CellType, df$Sample_well, df$mapping_error_score)

table(df$samples)

DefaultAssay(df)

#Create dataset with counts and genes
cts <- AggregateExpression(df, 
                           group.by = c("samples"),
                           assays = 'RNA',
                           slot = "counts",
                           return.seurat = FALSE)
cts <- cts$RNA
View(cts)

# transpose
cts.t <- t(cts)

# convert to data.frame
cts.t <- as.data.frame(cts.t)

# get values where to split cell type and sample
splitRows <- gsub('_.*', '', rownames(cts.t))


# split data.frame according row cell types which have been generated above
#cts.split <- split.data.frame(cts.t,
#f = factor(splitRows)) #only necessary when looking at a specific cell type

# fix colnames and transpose
cts.t.modified <- t(cts.t)
cts.t.modified[(1:10), (1:10)]


# Create Count matrix for DE
View(cts.t.modified)
#Analyse for all cells
# 1. generate sample level metadata
colData_all <- data.frame(samples = colnames(cts.t.modified))
View(colData_all)

#3. generate a variable for the condition to compare (e.g. HSC vs GMP)
colData_all <- colData_all %>%
  mutate(condition = ifelse(grepl('HSC', samples), 'HSC', 'earlyGMP')) 



# perform DESeq2 
# 1.Create DESeq2 object   
?DESeqDataSetFromMatrix
dds_all <- DESeqDataSetFromMatrix(countData = cts.t.modified,
                                  colData = colData_all,
                                  design = ~ condition)

# 2.filter -> optional to exclude genes with a low number of reads
keep <- rowSums(counts(dds_all)) >=10
dds_all <- dds_all[keep,]

# 3. run DESeq2
dds_all <- DESeq(dds_all)

# 4. Check the coefficients for the comparison
resultsNames(dds_all)


# 5. Generate results object
res_all <- results(dds_all)  #at this point you could determine the direction of the comparison, in this case resultsNames gets shows inly 1 condition 
res_all
res_all<-as.data.frame(res_all)

#Visualize the results
#Create a table
res_tbl_all <-res_all%>% 
  data.frame()%>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
res_tbl_all
#Save the table
write.csv(res_tbl_all,"/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/DE_healthy_cells_HSCvsearlyGMP.csv" )

#Filter for only significant Genes
#Set threshold
padj_cutoff <- 0.05
#Subset for signifcant results 
sig_res_all<- subset(res_tbl_all, res_tbl_all$padj<0.05)
write.csv(sig_res_all, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/DE_healthy_cells_HSCvsearlyGMP_pvalue_adj.csv" )


#Volcano Plot
#Create a vector where TRUE values denot padj ,0.05 and fold change >1.5
res_table_thres_all<- res_tbl_all %>% mutate(threshold=padj <0.05 & abs(log2FoldChange)>0.58)
ggplot(res_table_thres_all)+
  geom_point(aes(x=log2FoldChange, y= -log10(padj), color = threshold))+
  ggtitle("Volcano plot HSC vs earlyGMP")+
  xlab("log2 fold Change")+
  scale_y_continuous (limits = c(0,200)) +
  theme(legend.position = "none",
        plot.title = element_text (size =rel(1.5), hjust = 0.5),
        axis.title=element_text(size = rel(1.25)))
View(res_table_thres_all)
write.csv(res_table_thres_all, "/Users/rabeamecklenbrauck/Library/CloudStorage/OneDrive-Nexus365/Ivo-Ven-Project/Transcriptome analysis/Genes_GMP_HSC_healthycells_Volcano_Plot.csv")
