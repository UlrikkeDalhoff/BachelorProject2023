library(ape)
library(AnnotationDbi)
library(annotate)
library(base64enc)
library(BiocFileCache)
library(BiocIO)
library(BiocManager)
library(BiocParallel)
library(BiocVersion)
library(Biostrings)
library(BSgenome)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(carData)
library(cellranger)
library(cluster)
library(data.table)
library(DNAcopy)
library(edgeR)
library(factoextra)
library(fansi)
library(genefilter)
library(generics)
library(GenomicRanges)
library(GenomeInfoDbData)
library(GenomicAlignments)
library(GenomicFeatures)
library(GGally)
library(gtable)
library(jsonlite)
library(KEGGREST)
library(kmer)
library(locfit)
library(Matrix)
library(readr)
library(readxl)
library(Repitools)
library(rjson)
library(RSQLite)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(gridExtra)

library(rtracklayer)



#Data

#Import of my data:
my_data_iCLIP <- rtracklayer::import("/Users/ulrikkedalhoff/Desktop/BIOKEMI/år\ 3\ xP/Bachelor/data/E2mCh_110kDa_1.merged.unique.peaks.bed")

#Count up and plot how iCLIP sites there are on each chromosome. 

GenomicRanges::show(my_data_iCLIP)

length(seqnames(my_data_iCLIP))

#If you look at the data you can see that there are 6 chromosomes where we find iCLIP sites.

df <- annoGR2DF(my_data_iCLIP)

ggplot(df) + 
  geom_bar(colour = "black", aes(x=chr, fill = chr)) + theme_classic() + scale_fill_brewer(palette = "Set3") + xlab("Chromosomes")

#Count the number of iCLIP sites with overlap m6A sites and find the average distance between ECT2 iCLIP sites and m6A sites.

my_data_m6a <- rtracklayer::import("/Users/ulrikkedalhoff/Desktop/BIOKEMI/år\ 3\ xP/Bachelor/data/elife-49658-fig5-data3-v1_nanopore_m6A_sites.bed")

seqlevels(my_data_m6a, pruning.mode = "coarse") <- c("Chr1","Chr2","Chr3","Chr4","Chr5")

seqlevels(my_data_iCLIP)

GenomicRanges::findOverlaps(my_data_iCLIP, my_data_m6a)

overlaps <- GenomicRanges::countOverlaps(my_data_m6a, my_data_iCLIP, type = "any")

overlaps

#We find 883 overlaps between the datasets. 

Overlap_data <- subsetByOverlaps(my_data_m6a, my_data_iCLIP)

mean(nearest(Overlap_data))
# 442.1257 

#The average length between the iCLIP and m6A sites is 442.1257

#Make a plot centered on the m6A sites of the relative numbers of ECT2 iCLIP 
#peaks over distance on the x axis (that is how many iCLIP peaks are 1nt away 
#from an m6A site, 2nt ect, for -100 nt to +100 nt.)

GenomicRanges::findOverlaps(my_data_m6a, my_data_iCLIP)

overlap_results <- c()

for (i in -100:100){
  hits <- GenomicRanges::countOverlaps(shift(my_data_m6a, i), my_data_iCLIP)
  hits_plus <- hits[as.vector(strand(my_data_m6a)) == "+"]
  overlap_results <- c(overlap_results, sum(hits_plus))
}
overlap_results

xes <- -100:100

df_results <- data.frame(xes, overlap_results)
df_results

ggplot(df_results) +
  geom_line(aes(xes, overlap_results)) + theme_classic() + xlab("distance") + ggtitle("Distance for the forwarded strand")

overlap_results_2 <- c()

for (i in -100:100){
  hits_2 <- GenomicRanges::countOverlaps(shift(my_data_m6a, i), my_data_iCLIP)
  hits_minus <- hits_2[as.vector(strand(my_data_m6a)) == "-"]
  overlap_results_2 <- c(overlap_results_2, sum(hits_minus))
}
overlap_results_2

df_results_2 <- data.frame(xes, overlap_results_2)
df_results_2

ggplot(df_results_2) +
  geom_line(aes(xes, overlap_results_2)) + theme_classic() + xlab("distance") + ggtitle("Distance for the reverse strand")

#Using the Araport11 GTF annotations in the folder Data/annotations: 
  #How many unique genes are represented by the iCLIP and/or m6A sites and what is the average number of peaks/m6A sites per gene?
  #Make a barplot of the counts for the annotation type of the peaks - CDS, 3/5UTR.

my_data_Araport <- rtracklayer::import("/Users/ulrikkedalhoff/Desktop/BIOKEMI/år\ 3\ xP/Bachelor/data/Araport11_GTF_genes_transposons.current.gtf")

#How many unique genes are represented by the iCLIP and/or m6A sites and what is the average number of peaks/m6A sites per gene?

araport_iclip <- my_data_Araport[subjectHits(findOverlaps(my_data_iCLIP,my_data_Araport))]

my_data_iCLIP
my_data_Araport

length(unique(araport_iclip$gene_id))

#There are 2289 unique genes represented by iCLIP 

araport_m6a <- my_data_Araport[subjectHits(findOverlaps(my_data_m6a, my_data_Araport))]

length(unique(araport_m6a$gene_id))
araport_m6a

# there are 5231 unique genes represented by m6A.


#Make a barplot of the counts for the annotation type of the peaks - CDS, 3/5UTR.

sort_type_m6a <- c()
sort_type_m6a <- (araport_m6a[araport_m6a$type == "CDS" | araport_m6a$type == "three_prime_UTR" | araport_m6a$type == "five_prime_UTR"])

Df_peaks_m6a <- as.data.frame(sort_type_m6a)
Df_peaks_m6a

ggplot(Df_peaks_m6a) +
  geom_bar(colour = "black", aes(x=type, fill = type)) +
  theme_classic() + scale_fill_brewer(palette = "Set3") + xlab("Annotation types m6A")


#Using BSGenome::getSeq() and supplying the Arabidopsis genome and GenomicRanges 
#objects use the package to find the nucleotide (A,U,C,G) for all of the 
#individual m6A and ECT2 iCLIP positions and produce a plot of the counts for each nucleotide.


?BSgenome::getSeq

m6a_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), my_data_m6a)

iCLIP_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), my_data_iCLIP)


m6a_seq_df <- as.data.frame(m6a_seq)

iCLIP_seq_df <- as.data.frame(iCLIP_seq)

m6a_seq_df["x"][m6a_seq_df["x"] == "T"] <- "U"

ggplot(m6a_seq_df) +
  geom_bar(colour = "black", aes(x=x, fill = x)) + theme_classic() + xlab("m6A nucleotides") + 
  scale_fill_brewer(palette = "Set3") 

iCLIP_seq_df["x"][iCLIP_seq_df["x"] == "T"] <- "U"

ggplot(iCLIP_seq_df) +
  geom_bar(colour = "black", aes(x=x, fill = x)) + theme_classic() + xlab("iCLIP nucleotides") +
  scale_fill_brewer(palette = "Set3")
  


#Using BSGenome::getSeq() extract all the sequences in windows around ECT2 iCLIP 
#sites. Plot the average number of A,C,U,G as a function of increasing distance 
#from the iCLIP site.
#Here I adjust the window by changing the + n
extended_data_iCLIP <- my_data_iCLIP + 10

extended_data_iCLIP

length(unique(araport_iclip$gene_id))

extended_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_iCLIP)

extended_seq
matrix_seq <- as.matrix(extended_seq)


matrix_seq
columns <- c("A", "U", "G", "C")
df_nucleotides <- data.frame(matrix(nrow = 101, ncol = length(columns)))
colnames(df_nucleotides) <- columns



for (i in 1 : 101){
  hits_A <- mean(matrix_seq[, i] == "A")
  df_nucleotides$A[i] <- hits_A
  hits_U <- mean(matrix_seq[, i] == "T")
  df_nucleotides$U[i] <- hits_U
  hits_G <- mean(matrix_seq[, i] == "G")
  df_nucleotides$G[i] <- hits_G
  hits_C <- mean(matrix_seq[, i] == "C")
  df_nucleotides$C[i] <- hits_C
}

df_nucleotides$x <- -50 : 50

df_nucleotides

df_nucleotides %>%
  pivot_longer(!x)

df_nucleotides %>%
  pivot_longer(cols = !x) %>%
  ggplot(aes(x = x, y = value)) + 
  geom_line() +
  facet_wrap(~name, scale = "free") +
  ggtitle("Average of nucleotide as a function of increasing distance from the iCLIP site") +
  theme_linedraw() + 
  xlab("distance") + 
  ylab("Average")

ggplot(df_nucleotides) +
  geom_line(colour = "green",  aes(x = x, y = A)) +
  ylim(0, 0.8) + ylab("Average") + xlab("distance") + theme_linedraw() +
  geom_line(colour = "black", aes(x = x, y = U))  + 
  geom_line(colour = "blue", aes(x = x, y = G)) + 
  geom_line(colour = "orange", aes(x = x, y = C))

 #Plot with A = green, U = black, G = blue and C = orange  


#Find frequencies of 5-mers (use R-package kmer). Are there any which occur 
#more often than random chance?

#use the background dataset

background <- rtracklayer::import("/Users/ulrikkedalhoff/Desktop/BIOKEMI/år\ 3\ xP/Bachelor/data/iCLIP_ECT2_background_coordinates.bed")

#We have to change the seqnames to use the different methods
seqlevels(background, pruning.mode = "coarse") <- c("Chr1","Chr2","Chr3","Chr4","Chr5")

#Extend or else we cannot use the kmer count, since there will not be more than one nucleotide per position
background_extended <- background + 10

background_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), background_extended)

background_seq

background_seq_matrix <- as.matrix(background_seq)

background_seq_matrix

five_kmer_background <- kcount(background_seq_matrix, k = 5)

five_kmer_background

# The matrix_seq is the DNAStringset of our data converted into a matrix

five_kmer <- kcount(matrix_seq, k = 5)


five_kmer

?kmer

five_sub <- subset(five_kmer, colnames(five_kmer) %in% colnames(five_kmer_background))

sum(five_kmer[,2])

total_five <- colSums(five_kmer[,-1]) 

total_back <- colSums(five_kmer_background[,-1]) 

rest_size <- c(sum(total_five) - as.numeric(total_five)) 

rest_size_back <- c(sum(total_back) - as.numeric(total_back)) 


data <- data.frame(total = c( as.numeric(total_five)),
                   total_back = c(as.numeric(total_back)),
                   rest = c(rest_size),
                   rest_back = c(rest_size_back),
                   row.names = c(names(total_five)))

for (i in 1:nrow(data)){
  df <- data.frame(total = c(as.numeric(total_five[i]), as.numeric(total_back[i])),
                   rest = c(rest_size[i], rest_size_back[i]))
  data$test[i] <- fisher.test(df)$p.value
  }
colnames(data)[5] <- "p_value"

data

data$q_value <- p.adjust(data$p_value ,method = "fdr")

data_full <- subset(data, q_value <= 0.05)

data_full

data_ascending <- data_full[order(data_full$q_value),]

data_ascending_50 <- data_ascending[1:50,]

data_ascending_50

data_ascending_10 <- data_ascending[1:10,]

matrix_10_50 <- five_kmer[, (rownames(data_ascending_50))]

colnames(five_kmer)
rownames(data_full)
length(matrix_10_50)

matrix_10_50_background <- five_kmer_background[, (rownames(data_ascending_50))]


matrix_25_full_background


write.csv(matrix_10_50_background, file = "matrix_10_50_background.csv")


# Which criteria is there for choosing the most enriched 

# Using the expression (tags per million, TPM) file - for all expressed genes 
# in the file " arabidopsis_roots_shoots_average_TPM_genes.csv" (TPM>1) plot the 
# average expression for those genes which are vs are not targets of ECT2. 
# Ditto with m6A vs non-m6A genes.

TPM_data <- read.csv("/Users/ulrikkedalhoff/Desktop/BIOKEMI/år\ 3\ xP/Bachelor/data/arabidopsis_roots_shoots_average_TPM_genes.csv")

TPM_data

TPM_data_filter <- filter(TPM_data, value >= 1)

iclip_df_araport <- as.data.frame(araport_iclip$gene_id)

iclip_tpm1 <- unique(iclip_df_araport)

colnames(iclip_tpm1)[1] <- "gene"

iclip_tpm1

TPM_data_filter

iclip_tpm <- merge(TPM_data_filter, iclip_tpm1)

iclip_tpm

non_iclip_tpm <- anti_join(TPM_data_filter, iclip_tpm1)

non_iclip_tpm

average_df <- data.frame(gene = c("target_genes", "non_target_genes"),
                         average = c(mean(iclip_tpm$value), mean(non_iclip_tpm$value)))
average_df

ggplot(data = average_df, aes(x = gene, y = average, fill = gene)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ggtitle("TPM of target vs non_target of ECT2, iCLIP") +
  theme_classic() + 
  scale_fill_brewer(palette = "Set3")


m6a_df <- as.data.frame(araport_m6a$gene_id)

m6a_df_unique <- unique(m6a_df)

colnames(m6a_df_unique)[1] <- "gene"

m6a_tpm <- merge(TPM_data_filter, m6a_df_unique)

non_m6a_tpm <- anti_join(TPM_data_filter, m6a_df_unique)

average_m6a_df <- data.frame(gene = c("genes", "non_m6A_genes"),
                         average = c(mean(m6a_tpm$value), mean(non_m6a_tpm$value)))

ggplot(data = average_m6a_df, aes(x = gene, y = average, fill = gene)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  ggtitle("TPM of m6A vs non-m6A genes") +
  theme_classic() +
  scale_fill_brewer(palette = "Set3")

# to show the expression bias 










# Creating the training data for my model:
#adjust the window in the iclip data to 10, 25 and 50 and run all of the code below
extended_data_iCLIP <- my_data_iCLIP + 10

# Get the sequence
extended_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), extended_data_iCLIP)

#convert it into a matrix 
matrix_seq <- as.matrix(extended_seq)

# Get the 5mers
five_kmer <- kcount(matrix_seq, k = 5)

# Do the same for the background data
#Adjust the window to 10, 25, 50 and run all of the code below
background_extended <- background + 10

background_seq <- BSgenome::getSeq((BSgenome.Athaliana.TAIR.TAIR9), background_extended)

background_seq_matrix <- as.matrix(background_seq)

five_kmer_background <- kcount(background_seq_matrix, k = 5)


# Use the data and background to get the q_values. 

five_sub <- subset(five_kmer, colnames(five_kmer) %in% colnames(five_kmer_background))

total_five <- colSums(five_kmer[,-1]) 

total_back <- colSums(five_kmer_background[,-1]) 

rest_size <- c(sum(total_five) - as.numeric(total_five)) 

rest_size_back <- c(sum(total_back) - as.numeric(total_back)) 


data <- data.frame(total = c( as.numeric(total_five)),
                   total_back = c(as.numeric(total_back)),
                   rest = c(rest_size),
                   rest_back = c(rest_size_back),
                   row.names = c(names(total_five)))

for (i in 1:nrow(data)){
  df <- data.frame(total = c(as.numeric(total_five[i]), as.numeric(total_back[i])),
                   rest = c(rest_size[i], rest_size_back[i]))
  data$test[i] <- fisher.test(df)$p.value
}
colnames(data)[5] <- "p_value"

data

data$q_value <- p.adjust(data$p_value ,method = "fdr")

data_full <- subset(data, q_value <= 0.05)


#Create the datasets for the randomforest:
data_ascending <- data_full[order(data_full$q_value),]

#The 50 most significant hits
data_ascending_50 <- data_ascending[1:50,]

data_ascending_50

#The 10 most significant hits 
data_ascending_10 <- data_ascending[1:10,]

data_ascending_10


#Creating the data matrix by using data_ascending_10 if we want 10 most 
#significant hits, data_ascending_50 if we want 50 most significant hits or 
#data_full for all the significant hits
matrix_10_50 <- five_kmer[, (rownames(data_ascending_50))]


colnames(five_kmer)
rownames(data_full)
length(matrix_50_10)
ncol(matrix_10_10)

#Creating the background data matrix
matrix_10_50_background <- five_kmer_background[, (rownames(data_ascending_50))]



#save the datasets
write.csv(matrix_10_50_background, file = "matrix_10_50_background_new.csv")




