library(msa)
library(seqinr)
library(ggplot2)
library(adegenet)
library(ape)
library(ggtree)
library(DECIPHER)
library(viridis)
library(ggmsa)
library(rentrez)
library(Biostrings)


email <- "a00837960@tec.mx"
#wuhan, variante lambda, alpha, gamma, delta
acc_numbers <- c("NC_045512","OQ564806","OQ551287","OQ551274","OQ892331")
sequences <- entrez_fetch(db = "nucleotide", id = acc_numbers, rettype = "fasta", retmode = "text", email = email)
writeLines(sequences, "neucleotide.fasta")
DNA_seqs <- readDNAStringSet("neucleotide.fasta",format = "fasta")
DNA_seqs

cat("La longitud de la variante wuhan es", nchar(DNA_seqs[1]), "\n")
cat("La longitud de la variante lambda es", nchar(DNA_seqs[2]), "\n")
cat("La longitud de la variente alpha es", nchar(DNA_seqs[3]), "\n")
cat("La longitud de la variante gamma es", nchar(DNA_seqs[4]), "\n")
cat("La longitud de la variante  delta es", nchar(DNA_seqs[5]), "\n")

#SARS-COV-2 WUHAN________________________
sars_cov2 <- entrez_fetch(db="nucleotide", id="NC_045512.2", rettype="fasta", retmode="text", email = email)
writeLines(sars_cov2, "sars_cov2.fasta")
sars_cov2_seq <- read.fasta("sars_cov2.fasta")

#COVID
COV2_seq <- read.fasta(file = "sars_cov2.fasta")
COV2seq <- COV2_seq[[1]]
cuentas<-table(COV2seq)
#barplot(table(COV2seq))

barplot(table(COV2seq), main='Genoma del SARS-COV2',xlab='bases', ylab='frecuencia', col = c('blue','lightblue','#02F5B5','lightgreen'))

#COVID VARIANTE LAMBDA  ______________________
covid_lambda <- entrez_fetch(db="nucleotide", id="OQ564806", rettype="fasta", retmode="text", email = email)
writeLines(covid_lambda, "covid_lambda.fasta")
covid_lambda_seq <- read.fasta("covid_lambda.fasta")

#COVID
LAMBDA_seq <- read.fasta(file = "covid_lambda.fasta")
LAMBDAseq <- LAMBDA_seq[[1]]
cuentas<-table(LAMBDAseq)
#barplot(table(LAMBDAseq))

barplot(table(LAMBDAseq), main='Genoma del covid variante lambda',xlab='bases', ylab='frecuencia', col = c('blue','lightblue','#02F5B5','lightgreen'))

#COVID VARIANTE ALPHA  ______________________
covid_alpha <- entrez_fetch(db="nucleotide", id="OQ551287", rettype="fasta", retmode="text", email = email)
writeLines(covid_alpha, "covid_alpha.fasta")
covid_alpa_seq <- read.fasta("covid_alpha.fasta")

#COVID
ALPHA_seq <- read.fasta(file = "covid_alpha.fasta")
ALPHAseq <- ALPHA_seq[[1]]
cuentas<-table(ALPHAseq)
#barplot(table(ALPHAseq))

barplot(table(ALPHAseq), main='Genoma del covid variante alpha',xlab='bases', ylab='frecuencia', col = c('blue','lightblue','#02F5B5','lightgreen'))

#COVID VARIANTE GAMMA  ______________________
covid_gamma <- entrez_fetch(db="nucleotide", id="OQ551274", rettype="fasta", retmode="text", email = email)
writeLines(covid_gamma, "covid_gamma.fasta")
covid_gamma_seq <- read.fasta("covid_gamma.fasta")

#COVID
GAMMA_seq <- read.fasta(file = "covid_gamma.fasta")
GAMMAseq <- GAMMA_seq[[1]]
cuentas<-table(GAMMAseq)
#barplot(table(GAMMAseq))

barplot(table(GAMMAseq), main='Genoma del covid variante gamma',xlab='bases', ylab='frecuencia', col = c('blue','lightblue','#02F5B5','lightgreen'))

#COVID VARIANTE DELTA  ______________________
covid_delta <- entrez_fetch(db="nucleotide", id="OQ892331", rettype="fasta", retmode="text", email = email)
writeLines(covid_delta, "covid_delta.fasta")
covid_delta_seq <- read.fasta("covid_delta.fasta")

#COVID
DELTA_seq <- read.fasta(file = "covid_delta.fasta")
DELTAseq <- DELTA_seq[[1]]
cuentas<-table(DELTAseq)
#barplot(table(DELTAseq))

barplot(table(DELTAseq), main='Genoma del covid variante delta',xlab='bases', ylab='frecuencia', col = c('blue','lightblue','#02F5B5','lightgreen'))

email <- "a00837960@tec.mx"
#wuhan, variante lambda, alpha, gamma, delta
acc_numbers <- c("NC_045512","OQ564806","OQ551287","OQ551274","OQ892331")
sequences <- entrez_fetch(db = "nucleotide", id = acc_numbers, rettype = "fasta", retmode = "text", email = email)
writeLines(sequences, "neucleotide.fasta")
DNA_seqs <- readDNAStringSet("neucleotide.fasta",format = "fasta")
DNA_seqs

se1_1 <- as.DNAbin(DNA_seqs[[1]])
se1_2 <- as.DNAbin(DNA_seqs[[2]])
se1_3 <- as.DNAbin(DNA_seqs[[3]])
se1_4 <- as.DNAbin(DNA_seqs[[4]])
se1_5 <- as.DNAbin(DNA_seqs[[5]])

GC1 <- (GC.content(se1_1)*100)
cat("El porcentaje de GC de variante 1, wuhan: NC_045512, es:", GC1, "%\n")
GC2 <- (GC.content(se1_2)*100)
cat("El porcentaje de GC de variante 2, lambda: OQ564806, es:", GC2, "%\n")
GC3 <- (GC.content(se1_3)*100)
cat("El porcentaje de GC de variante 3, alpha: OQ551287,  es:", GC3, "%\n")
GC4 <- (GC.content(se1_4)*100)
cat("El porcentaje de GC de variante 4, gamma: OQ551274,  es:", GC4, "%\n")
GC5 <- (GC.content(se1_5)*100)
cat("El porcentaje de GC de variante 5, delta: OQ892331,  es:", GC5, "%\n")

contrasentido_wuhan <- reverseComplement(DNAString(DNA_seqs[[1]]))
contrasentido_alpha <- reverseComplement(DNAString(DNA_seqs[[2]]))
contrasentido_gamma <- reverseComplement(DNAString(DNA_seqs[[3]]))
contrasentido_beta <- reverseComplement(DNAString(DNA_seqs[[4]]))
contrasentido_delta <- reverseComplement(DNAString(DNA_seqs[[5]]))
contrasentido_wuhan
contrasentido_alpha
contrasentido_gamma
contrasentido_beta
contrasentido_delta


