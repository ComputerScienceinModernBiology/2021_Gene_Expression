library(tximport) #exporting the GFF file for processing by kallisto 
library(GenomicFeatures) # basic package needed for lot of genomic analysis
library(AnnotationDbi) # basic package for all annotation packages 
library(DESeq2) #DEG processing package 
library(pheatmap) #making heatmap

#Concatenating all the count file in same variable named samples
samples <- c("L1.tsv","L2.tsv", "L3.tsv","P1.tsv","P2.tsv","P3.tsv")

#### Assigning the heading to the count file 
names(samples) <- c("Control_1","Control_2","Control_3","Pten_KO_1","Pten_KO_2","Pten_KO_3")
######## Assigning the condition 
sampleTable <- data.frame(condition = factor(c("Control","Control","Control","Pten_KO","Pten_KO","Pten_KO")))
TxDB <- makeTxDbFromGFF("./gencode.vM25.annotation.gff3.gz") #extracting genome feature format file 
pk <- keys(TxDB, keytype = "TXNAME") # took all the transcript name from the table filtering all the exons, introns and other junk  
tx2gene1<- select(TxDB, pk, "GENEID", "TXNAME") # taking pk list for transcript, selecting corresponding genes from TxDb and creating the table with all the transcript and relating the genes associated with it.
txi <- tximport(samples, type = "kallisto", tx2gene =  tx2gene1, 
                ignoreAfterBar = TRUE)## using samples(text file with location) go through this and find all the kasllisto file. If you are using salmon u have to chnage the argument .We are gonna summarize the information on gene level using tx2gene argument. 
#The ignoreAfterbar command is essential for taking only the transcript id before the first bar or first separator. 
## specific data i.e. transcript in this case. This is essential for mapping 1:1 gene to trnascript id.
row.names(sampleTable) <- colnames(txi$counts) # take the column name from object txi and make it row name 

head(txi$counts) # head commands just shows the first few lines of the result 
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition) # starts the differential equation test 
#changing the object type from txi object to desq object since it will be enriched with much more information by linking 
# all the data from two tables 
# ~condition guides the deseq to do differential expression testing between the condition 
#Define factor levels here
dds$condition <- factor(dds$condition, levels = c("Control","Pten_KO")) #### set factor so the orientation or order is constant 
row.names(dds) <- substr(row.names(dds), 1,18) ### Extracts the ensemble id with 18 character and discards remaining 
#character after period  
as.data.frame(colData(dds)) #helps to identify the dds matrix of comparison 
dds <- DESeq(dds) # actual differential expression testing 

########### heatmap 
counttable<-counts(dds,normalized=TRUE)
y<-counttable[c("ENSMUSG00000013663", 
                "ENSMUSG00000005087", "ENSMUSG00000004655"),]
pheatmap(y,cluster_cols = FALSE,scale = "row")

######################
Epithelial_Cell_Proliferation <- df[c("Cdk6",
"Gata3",
"Pten",
"Ptn",
"Thbs1",
"Pbld1",
"Aqp1",
"Slit2",
"Bmper",
"Cenpv",
"Fzd7",
"Cd44",
"Inhba",
"Nrk",
"Ano1",
"Piezo2",
"Adam12",
"Adgrb2",
"Igf1r",
"Ltbp3"),]

pheatmap(Epithelial_Cell_Proliferation,cluster_cols = F,scale = "row",
         main = "Negative Regulation of Epithelial Cell Proliferation ", 
         cutree_rows = 1,
         gaps_col = c(3),
         cellheight = 10)


############## End of the code 