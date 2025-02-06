#library("import")
#library("knitr")
#if (!require("BiocManager", quietly = TRUE))
  # install.packages("BiocManager")
#BiocManager::install("BiocStyle")
library("biomformat")
library("here")
library("BiocStyle")
#BiocManager::install("ggplot2")
library("ggplot2")
#BiocManager::install("gridExtra")
library("gridExtra")
#BiocManager::install("Rcpp")
library("Rcpp")
#BiocManager::install("dada2")
library("dada2")
#BiocManager::install("phyloseq")
library("phyloseq")
#BiocManager::install("DECIPHER")
library("DECIPHER")
#BiocManager::install("ape")
library("ape")
#BiocManager::install("phangorn")
library("phangorn")
#BiocManager::install("ShortRead")
library("ShortRead")
#install.packages("remotes")
#remotes::install_github("MadsAlbertsen/ampvis2")
library("ampvis2")
#BiocManager::install("ggpubr")
#library("ggpubr")
#remotes::install_github("vmikk/metagMisc")
library("metagMisc")
#BiocManager::install("vegan")
library("vegan")
#BiocManager::install("plyr")
library("plyr")
#BiocManager::install("dplyr")
library("dplyr")
#BiocManager::install("MASS")
#library("MASS")
#BiocManager::install("tibble")
library("tibble")
#BiocManager::install("microbiome")
#library("microbiome")
#install.packages("devtools")
#devtools::install_github("microsud/microbiomeutilities")
#library("microbiomeutilities")
#BiocManager::install("reshape2")
library("reshape2")
#BiocManager::install("ggplot2")
library("ggplot2")
#BiocManager::install("viridis")
#library("viridis")
#BiocManager::install("hrbrthemes")
#library("hrbrthemes")
#library("devtools")
#install_github("zdk123/SpiecEasi")
#library("SpiecEasi")
#BiocManager::install("igraph")
#library("igraph")
#BiocManager::install("ALDEx2")
#library("ALDEx2")
#BiocManager::install("metagenomeSeq")
#library("metagenomeSeq")
#BiocManager::install("DESeq2")
#library("DESeq2")


# Reading files
fastq_files= "/work/samodha/sadams23/jl/Fastq/" #here you will give the folder name where your fastq files are.
list.files(fastq_files) #listing fastq files


#filtering and trimming fastq files

fnFs <- sort(list.files(fastq_files, pattern="_R1_001.fastq.gz")) #forward reads
fnRs <- sort(list.files(fastq_files, pattern="_R2_001.fastq.gz")) #reverse reads

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sampleNames <- sapply(strsplit(fnFs, "_"), `[`, 1)

# Specify the full path to the fnFs and fnRs


fnFs <- file.path(fastq_files, fnFs)
fnRs <- file.path(fastq_files, fnRs)


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
fnFs[1:3]
fnRs[1:3]

#quality plot (foward) to look at Q score

#The first two forward reads:

plotQualityProfile(fnFs[1:3])
#can change to view more than two plots at a time example: plotQualityProfile(fnFs[1:10]) 
###be sure to know if you are looking at a negative control! Your Q score plots should not look as "good" on a negative vs positive vs sample

#quality plot (reverse) to look at Q score

#The first two reverse reads:

plotQualityProfile(fnRs[1:2])
#overlapped beween forward and reverse reads
#threshold
#can change to view more than two plots at a time example: plotQualityProfile(fnRs[1:10])


#trimming and filtering the F/R reads
#This can take some time depending on your computer speed, memory, and samples

filt_path <- file.path(fastq_files, "filtered") 
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

#can trim the reads according to what fits your dataset by changing the numbers in 'truncLen=c(240,160)'
#additionally you can change your maxEE based off of what you are wanting (based off of expected errors)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,230),
		                         maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
					                   compress=TRUE, multithread=TRUE, matchIDs = T)
#truncLen = C(240,150): truncating forward and reverse reads to desired length
#N: can not decide which nucleotides
#maxEE: maximum expected errors: based on q score, likelihode to encounters sequence
#maxEE = c(2 -> forward,2 -> reverse reads): maximum 2 errors per sequence
#truncQ = 2: in fastq file, 2 indicates unsure (low quality) => removes that
#rm.phix: remove phix contaminants
#multithread (which to False if run in Windows)

out
fnFs
fnRs

#Data Statistics after Trimming

sum(out[,1]) #total reads in---746,124
sum(out[,2]) #total reads out--- 586,255
sum(out[,1]) - sum(out[,2]) #reads lost--- 159,869
sum(out[,2])/sum(out[,1]) # percentage data retained -- 79%--this number is not ideal 

#to avoid error, due to very low reads

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

#Dereplication

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names

names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

errF <- learnErrors(filtFs, multithread=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)


plotErrors(errF)
plotErrors(errR)

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=T)
head(mergers[[1]])
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]
#dim(seqtab2)

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

#seqtab.nochim2 <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
#dim(seqtab.nochim2)

saveRDS(seqtab, "JL.rds")
save(out,dadaFs,dadaRs,errF,errR,mergers,derepFs,derepRs,seqtab, file="JL_Samples_Study.RData")
#load("Selby_Calf_Study_Samples_Selby_Calf_Study.RData")

taxa <- assignTaxonomy(seqtab.nochim, "/work/samodha/sadams23/jl/silva_nr99_v138_wSpecies_train_set.fa.gz", multithread=TRUE)

save(seqtab.nochim,seqtab,taxa, file="JL_Samples_seqtab.nonchim_taxa.RData")


#applying the learned error rates

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]] #summary of first sample

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, file = "Processing_Reads_new.csv")


t(seqtab.nochim) %>%
	  make_biom() %>%
	    write_biom(here::here("table-run1.biom"))

uniquesToFasta(seqtab.nochim, fout = here::here("rep-seqs-run1.fna"), 
	                  ids = colnames(seqtab.nochim))

#giving our seq headers more manageable names (ASV_1, ASV_2...)

#dir.create("Analysis")

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
	  asv_headers[i] <- paste(">ASV", i, sep="_")
}

 # making and writing out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "Analysis/ASVs.fa")

 # count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "Analysis/ASVs_counts.txt", sep="\t", quote=F)

 # tax table:
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "Analysis/ASVs_taxonomy.txt", sep="\t", quote=F)
write.table(taxa, "Analysis/taxonomy.txt", sep="\t", quote=F)


#commands on mothur to make tree


#mothur
#mothur > align.seqs(fasta=ASVs.fa, reference=silva.nr_v138_1.align, processors=8)

#summary.seqs(fasta=ASVs.align, processors=8)

#mothur > summary.seqs(fasta=ASVs.align, processors=8)
#mothur > dist.seqs(fasta= ASVs.align,output=lt, processors=8)
#mothur > clearcut(phylip=ASVs.phylip.dist)

#Command on mothur to make tree

#Chmod +x mothur.sh
#bash mothur.sh

#final phyloseq object
tree= read_tree("Analysis/ASVs.phylip.tre")
library("readxl")
JL_pig_Study_meta= readxl::read_xlsx("Mapping_JL.xlsx")
JL_pig_Study_meta= as.data.frame(JL_pig_Study_meta)
rownames(JL_pig_Study_meta) =JL_pig_Study_meta$SampleID
summary(JL_pig_Study_meta)

ps= phyloseq(otu_table(asv_tab, taxa_are_rows=T),
	                    sample_data(JL_pig_Study_meta),
			                   tax_table(asv_tax),
			                   phy_tree(tree))

save(ps, file= "ps.RData")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#	    install.packages("BiocManager")
 
#BiocManager::install("genefilter")

#library("genefilter")
#set the function parameters, 0.15% abundance w/in a sample, and must have that criteria in at least 2 samples
#flist<- filterfun(kOverA(2, 0.0015))
#Use function on your data:
#this should be performed on your asv table transformed to proportional data

#dd2.logi <- filter_taxa(dd2.ps.prop, flist) #create a list of ASVs that meet flist criteria
#Now filter out ASVs that do not meet the criteria kOverA i.e. dd2.logi list...
#dd2.filt.ps = prune_taxa(dd2.logi, dada_pc_ps)
#taxa_sums(dd2.filt.ps)
