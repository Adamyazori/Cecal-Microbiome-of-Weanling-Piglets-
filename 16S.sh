#!/bin/sh
#SBATCH --licenses=common
#SBATCH --partition=samodha
#SBATCH --ntasks-per-node=4
#SBATCH --mem=100GB
#SBATCH --nodes=1
#SBATCH --time=168:00:00
#SBATCH --job-name=16S_analysis
#SBATCH --error=download.%J.stdout
#SBATCH --output=download.%J.stderr
#SBATCH --mail-user=adamyazori@gmail.com
#SBATCH --mail-type=ALL

module load qiime2/2022.2
module load picrust2/2.4

echo This is Seidu Adams 16S Metagenome analysis pipeline


echo import sequences to qiime format
qiime tools import \
	                 --type 'SampleData[PairedEndSequencesWithQuality]' \
	                             --input-path Fastq \
	                                           --input-format CasavaOneEightSingleLanePerSampleDirFmt \
	                                                           --output-path demux-paired-end.qza

echo visualization of the input data
qiime demux summarize \
	                 --i-data demux-paired-end.qza \
	                             --o-visualization demux-paired-end.qzv

echo Quality control
qiime dada2 denoise-paired \
	                 --i-demultiplexed-seqs demux-paired-end.qza \
	                             --p-trim-left-f 15 \
	                                           --p-trim-left-r 15 \
	                                                           --p-trunc-len-f 250 \
	                                                                             --p-trunc-len-r 230 \
	                                                                                                 --o-table table.qza \
	                                                                                                                       --o-representative-sequences rep-seqs.qza \
	                                                                                                                                               --o-denoising-stats denoising-stats.qza
echo To count statistics
qiime metadata tabulate \
	                 --m-input-file denoising-stats.qza \
	                             --o-visualization stats-dada2.qzv
echo To visualize table Rep seq
qiime feature-table summarize \
	                 --i-table table.qza \

qiime feature-table tabulate-seqs \
                 --i-data rep-seqs.qza \
                             --o-visualization rep-seqs.qzv

echo To export OTU table or ASV
mkdir phyloseq
qiime tools export \
	               --input-path table.qza \
	                       --output-path phyloseq

echo To Convert biom format to tab-separated text format:
biom convert \
	               -i phyloseq/feature-table.biom \
	                       -o phyloseq/otu_table.tsv \
	                              --to-tsv

echo Modify otu_table.txt to make it easier to read into R. Use sed to delete the first line and "#OTU ID" from the. second line
cd phyloseq
sed -i '1d' otu_table.tsv
sed -i 's/#OTU ID//' otu_table.tsv
cd ../

echo Export representative sequences:
qiime tools export \
	               --input-path rep-seqs.qza \
	                       --output-path phyloseq

