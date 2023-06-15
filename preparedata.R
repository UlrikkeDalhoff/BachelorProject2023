# Preparing of data
library(ape)
library(BiocGenerics)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(GenomicRanges)
library(tidyverse)
library(ggplot2)
library(Repitools)
library(rtracklayer)
library(data.table)
library(kmer)

#NB requires BiocManager 
# install.packages("BiocManager")
list.of.packages <- c("ape", "BSgenome", "BSgenome.Athaliana.TAIR.TAIR9", "GenomicRanges", "tidyverse", "ggplot2", "rtracklayer", "BiocParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

#Importing the iCLIP data
my_data_iCLIP <- rtracklayer::import("D:/Data_bachelor/iCLIP.bed")


#Importing the m6a data
my_data_m6a <- rtracklayer::import("D:/Data_bachelor/m6a_nanopore.bed")

#Make the datasets comparable 
seqlevels(my_data_m6a, pruning.mode = "coarse") <- c("Chr1","Chr2","Chr3","Chr4","Chr5")


#Importing the araport data 
my_data_Araport <- rtracklayer::import("D:/Data_bachelor/Araport11.gtf") 

#Finding the unique genes for iclip by using the araport data
araport_iclip <- my_data_Araport[subjectHits(findOverlaps(my_data_iCLIP,my_data_Araport))]

length(unique(araport_iclip$gene_id))



#Finding unique genes for m6a by using araport data
araport_m6a <- my_data_Araport[subjectHits(findOverlaps(my_data_m6a, my_data_Araport))]

length(unique(araport_m6a$gene_id))

#Import the background data
background <- rtracklayer::import("D:/Data_bachelor/background.bed")

#make the background comparable 
seqlevels(background, pruning.mode= "coarse") <- c("Chr1","Chr2","Chr3","Chr4","Chr5")


#making the groups for my 5 fold cross validation
genes_iclip <- my_data_Araport[nearest(my_data_iCLIP, my_data_Araport)]$gene_id

genes_background <- my_data_Araport[nearest(background, my_data_Araport)]$gene_id


genes_all <- c(genes_iclip, genes_background)

sample_vec <- sample(1:5, length(unique(genes_all)), replace = TRUE)

unique_genes <- unique(genes_all)
names(sample_vec) <- unique_genes
my_group_vector <- sample_vec[genes_all]
my_group_vector

write.csv(my_group_vector, file ="my_groups.csv")




#Making a + - 50 nt window for iCLIP
extended_data_iCLIP <- my_data_iCLIP + 50

extended_data_iCLIP

#Getting the sequence
extended_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_iCLIP)

#Converting the sequence to a matrix
matrix_seq <- as.matrix(extended_seq)

#counting the 5mers
five_kmer <- kcount(matrix_seq, k = 5)
five_kmer












# Preparing the hyperTRIBE data
hypertribe_roots <- readxl::read_excel("D:/hypertribe.xlsx", sheet = 3)
hypertribe_roots

hypertribe_aerial <- readxl::read_excel("D:/hypertribe.xlsx", sheet = 5)

#Combine the two datasets so we have one set
hypertribe <- rbind(hypertribe_roots, hypertribe_aerial)
length(hypertribe$GENE)


#Sort the hypertribe, so it does not contain any genes from the iclip data
no_iclip_hypertribe <- hypertribe[!hypertribe$GENE %in% genes_iclip, ]

no_iclip_hypertribe

hyper_araport <- my_data_Araport[my_data_Araport$gene_id %in% no_iclip_hypertribe$GENE]
hyper_araport

#Subsetting the Araport11 to only contain the hits in 3' UTR
three_utr_araport <- my_data_Araport[my_data_Araport$type == "three_prime_UTR", ]
three_utr_araport



#Extracting the values in three_utr_araport by using the genes from no_iclip_hypertribe 
hypertribe_araport <- three_utr_araport[three_utr_araport$gene_id %in% no_iclip_hypertribe$GENE, ]
hypertribe_araport

unique_genes <- unique(hypertribe_araport$gene_id)
unique_genes

#Dataen splittes i gener
gene_split <- split(hypertribe_araport, hypertribe_araport$gene_id)

#Vi gemmer kun de unikke (1 af hver)
unique_genes_hyper <- lapply(gene_split, function(x) x[1])

#making the dataset so it only contain unique genes. 
unique_hypertribe <- do.call(c, unique_genes_hyper)


width(unique_hypertribe[["AT1G01020"]])


unique_hypertribe[["AT1G01020"]]

#creating an empty grangeobject for my results



#GRanges(
#  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
#  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#  score = 1:10,
#  GC = seq(1, 0, length=10))

start(unique_hypertribe[["AT1G01020"]])

#Create the position loop 
for (i in 1:width(unique_hypertribe[["AT1G01020"]])){
  position <- unique_hypertribe[start(unique_hypertribe[["AT1G01020"]]) = start(unique_hypertribe[["AT1G01020"]]),
                                end(unique_hypertribe[["AT1G01020"]]) = (start(unique_hypertribe[["AT1G01020"]]) + i),
                                width((unique_hypertribe[["AT1G01020"]]) = i)]
  print(position)
  
}
gene_1








TPM_data <- read.csv("D:new_data_bachelor/arabidopsis_roots_shoots_average_TPM_genes.csv")

TPM_data

TPM_data_filter <- filter(TPM_data, value >= 1)

hyper_df_araport <- as.data.frame(hypertribe_araport$gene_id)

hyper_tpm1 <- unique(hyper_df_araport)

colnames(hyper_tpm1)[1] <- "gene"

hyper_tpm1

TPM_data_filter

hyper_tpm <- merge(TPM_data_filter, hyper_tpm1)

hyper_tpm

non_hyper_tpm <- anti_join(TPM_data_filter, hyper_tpm1)

non_hyper_tpm

average_df <- data.frame(gene = c("target_genes", "non_target_genes"),
                         average = c(mean(hyper_tpm$value), mean(non_hyper_tpm$value)))
average_df

ggplot(data = average_df, aes(x = gene, y = average, fill = gene)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set3") +
  ggtitle("TPM of target vs non_target of ECT2, HyperTRIBE") +
  theme_classic() 
  



















gene_1 <- unique_hypertribe[["AT1G01020"]]


#See if there is an overlap in the m6A data 
length(subjectHits(GenomicRanges::findOverlaps(gene_1, araport_m6a)))
# result 0

width(gene_1)

start_pos <- start(gene_1)
end_pos <- end(gene_1)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_1),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_1))
  results <- c(results, subinterval)
}

print(results)




results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf10_gene1.csv" )











#All of the above for gene AT1G01080

width(unique_hypertribe[["AT1G01080"]])


unique_hypertribe[["AT1G01080"]]

#creating an empty grangeobject for my results



start(unique_hypertribe[["AT1G01080"]])


gene_2 <- unique_hypertribe[["AT1G01080"]]

gene_2

#See if there is an overlap in the m6A data 
#length(subjectHits(GenomicRanges::findOverlaps(gene_2, araport_m6a)))
# result 20


width(gene_2)

start_pos <- start(gene_2)
end_pos <- end(gene_2)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_2),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_2))
  results <- c(results, subinterval)
}

print(results)




results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 50


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf50_gene2cut.csv" )













#All of the above for gene AT2G01060



width(unique_hypertribe[["AT2G01060"]])


(unique_hypertribe[["AT2G01060"]])



#creating an empty grangeobject for my results



start(unique_hypertribe[["AT2G01060"]])


gene_3 <- unique_hypertribe[["AT2G01060"]]

gene_3

#See if there is an overlap in the m6A data 
length(subjectHits(GenomicRanges::findOverlaps(gene_3, araport_m6a)))
# result 0

width(gene_3)

start_pos <- start(gene_3)
end_pos <- end(gene_3)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_3),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_3))
  results <- c(results, subinterval)
}

#results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

#extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

#hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf10_gene3.csv" )












#All of the above for gene AT2G40400



width(unique_hypertribe[["AT2G40400"]])


(unique_hypertribe[["AT2G40400"]])



#creating an empty grangeobject for my results



start(unique_hypertribe[["AT2G40400"]])


gene_4 <- unique_hypertribe[["AT2G40400"]]

gene_4

#See if there is an overlap in the m6A data 
length(subjectHits(GenomicRanges::findOverlaps(gene_4, araport_m6a)))
# result 5


width(gene_4)

start_pos <- start(gene_4)
end_pos <- end(gene_4)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_4),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_4))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf10_gene4.csv" )










#All of the above for gene AT3G02460



width(unique_hypertribe[["AT3G02460"]])


(unique_hypertribe[["AT3G02460"]])



#creating an empty grangeobject for my results



start(unique_hypertribe[["AT3G02460"]])


gene_5 <- unique_hypertribe[["AT3G02460"]]

gene_5


#See if there is an overlap in the m6A data 
length(subjectHits(GenomicRanges::findOverlaps(gene_5, araport_m6a)))
# result 16


width(gene_5)

start_pos <- start(gene_5)
end_pos <- end(gene_5)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_5),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_5))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf10_gene5.csv")







#All of the above for gene AT3G22220



width(unique_hypertribe[["AT3G22220"]])


(unique_hypertribe[["AT3G22220"]])



#creating an empty grangeobject for my results



start(unique_hypertribe[["AT3G22220"]])


gene_6 <- unique_hypertribe[["AT3G22220"]]

gene_6

#See if there is an overlap in the m6A data 
length(subjectHits(GenomicRanges::findOverlaps(gene_6, araport_m6a)))
# result 0


width(gene_6)

start_pos <- start(gene_6)
end_pos <- end(gene_6)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_6),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_6))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf10_gene6.csv")








#All of the above for gene AT4G11850



width(unique_hypertribe[["AT4G11850"]])


(unique_hypertribe[["AT4G11850"]])



#creating an empty grangeobject for my results



start(unique_hypertribe[["AT4G11850"]])


gene_7 <- unique_hypertribe[["AT4G11850"]]

gene_7

#See if there is an overlap in the m6A data 
length(subjectHits(GenomicRanges::findOverlaps(gene_7, araport_m6a)))
# result 0

width(gene_7)

start_pos <- start(gene_7)
end_pos <- end(gene_7)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_7),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_7))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf10_gene7.csv")







#All of the above for gene AT4G39850



width(unique_hypertribe[["AT4G39850"]])


(unique_hypertribe[["AT4G39850"]])



#creating an empty grangeobject for my results



start(unique_hypertribe[["AT4G39850"]])


gene_8 <- unique_hypertribe[["AT4G39850"]]

gene_8

#See if there is an overlap in the m6A data 
length(subjectHits(GenomicRanges::findOverlaps(gene_8, araport_m6a)))
# result 39

width(gene_8)

start_pos <- start(gene_8)
end_pos <- end(gene_8)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_8),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_8))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf10_gene8.csv")






#All of the above for gene AT5G50310



width(unique_hypertribe[["AT5G50310"]])


(unique_hypertribe[["AT5G50310"]])



#creating an empty grangeobject for my results



start(unique_hypertribe[["AT5G50310"]])


gene_9 <- unique_hypertribe[["AT5G50310"]]

gene_9

#See if there is an overlap in the m6A data 
length(subjectHits(GenomicRanges::findOverlaps(gene_9, araport_m6a)))
# result 0 


width(gene_9)

start_pos <- start(gene_9)
end_pos <- end(gene_9)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_9),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_9))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf10_gene9.csv")




#All of the above for gene AT5G67640



width(unique_hypertribe[["AT5G67640"]])


(unique_hypertribe[["AT5G67640"]])



#creating an empty grangeobject for my results



start(unique_hypertribe[["AT5G67640"]])


gene_10 <- unique_hypertribe[["AT5G67640"]]

gene_10

#See if there is an overlap in the m6A data 
length(subjectHits(GenomicRanges::findOverlaps(gene_10, araport_m6a)))
# result 0 


width(gene_10)

start_pos <- start(gene_10)
end_pos <- end(gene_10)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_10),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_10))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_hypertribe_tester <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_hypertribe_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_hypertribe_tester)

extended_hypertribe_seq



#Converting the sequence into a matrix 
hypertribe_seq_matrix <- as.matrix(extended_hypertribe_seq)

hypertribe_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
hypertribe_kmer <- kcount(hypertribe_seq_matrix, k=5)


write.csv(hypertribe_kmer, file = "Hypertribe_for_rf10_gene10.csv")











#Creating the data that are non-target genes
no_iclip_araport <- my_data_Araport[!my_data_Araport$gene_id %in% araport_iclip$gene_id, ]

no_hyper_araport <- no_iclip_araport[!no_iclip_araport$gene_id %in% hyper_araport$gene_id, ]

no_m6a_araport <- no_hyper_araport[!no_hyper_araport$gene_id %in% araport_m6a$gene_id, ]


no_m6a_araport




#Extracting the values in three_utr_araport
araport_three <- no_m6a_araport[no_m6a_araport$type == "three_prime_UTR", ]
araport_three

unique_genes_a <- unique(araport_three$gene_id)
unique_genes_a

#Dataen splittes i gener
gene_split_araport <- split(araport_three, araport_three$gene_id)

#Vi gemmer kun de unikke (1 af hver)
unique_genes_araport <- lapply(gene_split_araport, function(x) x[1])

#making the dataset so it only contain unique genes. 
unique_araport <- do.call(c, unique_genes_araport)

unique_araport






#Creating first test gene from araport [AT1G12070] with model 10
width(unique_araport[["AT1G12070"]])


(unique_araport[["AT1G12070"]])


start(unique_araport[["AT1G12070"]])


gene_test1 <- unique_araport[["AT1G12070"]]

gene_test1

width(gene_test1)

start_pos <- start(gene_test1)
end_pos <- end(gene_test1)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_test1),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_test1))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_araport <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_araport_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_araport)

extended_araport_seq



#Converting the sequence into a matrix 
araport_seq_matrix <- as.matrix(extended_araport_seq)

araport_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
araport_kmer <- kcount(araport_seq_matrix, k=5)


write.csv(araport_kmer, file = "araport_for_rf_gene_test1_10.csv")



#Creating first test gene from araport [AT1G80810] with model 10

width(unique_araport[["AT1G80810"]])


(unique_araport[["AT1G80810"]])


start(unique_araport[["AT1G80810"]])


gene_test2 <- unique_araport[["AT1G80810"]]

gene_test2

width(gene_test2)

start_pos <- start(gene_test2)
end_pos <- end(gene_test2)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_test2),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_test2))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_araport <- results + 10


#Then I obtain the sequence from the BSgenome tool 
extended_araport_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_araport)

extended_araport_seq



#Converting the sequence into a matrix 
araport_seq_matrix <- as.matrix(extended_araport_seq)

araport_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
araport_kmer <- kcount(araport_seq_matrix, k=5)


write.csv(araport_kmer, file = "araport_for_rf_gene_test2_10.csv")






# using the 50 model
#Creating first test gene from araport [AT1G12070] with model 10
width(unique_araport[["AT1G12070"]])


(unique_araport[["AT1G12070"]])


start(unique_araport[["AT1G12070"]])


gene_test50_1 <- unique_araport[["AT1G12070"]]

gene_test50_1

width(gene_test50_1)

start_pos <- start(gene_test50_1)
end_pos <- end(gene_test50_1)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_test50_1),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_test50_1))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_araport <- results + 50


#Then I obtain the sequence from the BSgenome tool 
extended_araport_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_araport)

extended_araport_seq



#Converting the sequence into a matrix 
araport_seq_matrix <- as.matrix(extended_araport_seq)

araport_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
araport_kmer <- kcount(araport_seq_matrix, k=5)


write.csv(araport_kmer, file = "araport_for_rf_gene_test1_50.csv")








#Creating first test gene from araport [AT1G80810] with model 50

width(unique_araport[["AT1G80810"]])


(unique_araport[["AT1G80810"]])


start(unique_araport[["AT1G80810"]])


gene_test2_50 <- unique_araport[["AT1G80810"]]

gene_test2_50

width(gene_test2_50)

start_pos <- start(gene_test2_50)
end_pos <- end(gene_test2_50)

results <- GRanges(
  seqnames = character(),
  ranges = IRanges(start = integer(), end = integer()),
  strand = character()
)

for (i in start_pos:end_pos){
  subinterval <- GRanges(seqnames = seqnames(gene_test2_50),
                         ranges = IRanges(start = i, end = i),
                         strand = strand(gene_test2_50))
  results <- c(results, subinterval)
}

results

#Preparing the hypertribe data so it can be used in my rf model
#Choosing a window there is + - 50 nt, since it gave the best results in the 
#model training
extended_data_araport <- results + 50


#Then I obtain the sequence from the BSgenome tool 
extended_araport_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_araport)

extended_araport_seq



#Converting the sequence into a matrix 
araport_seq_matrix <- as.matrix(extended_araport_seq)

araport_seq_matrix

#Create 5mers of the sequence by using the kmer package
#This is our final data I will use to predict new binding sites by using my RF
#model. 
araport_kmer <- kcount(araport_seq_matrix, k=5)


write.csv(araport_kmer, file = "araport_for_rf_gene_test2_50.csv")

