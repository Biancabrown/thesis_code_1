##Host Phylogeny 
###Alignment 
~~~
Alignment was build in geneious 
~~~
###Jmodel Test
~~~
J model test was used for model selection TRCAT
~~~
###Sanger RaxML 1000 boot strap 
~~~
#!/bin/bash

# Request an hour of runtime:
#SBATCH --time=1:00:00
# Default resources are 1 core with 2.8GB of memory.
# Use more memory (4GB):
#SBATCH --mem=4G
# Specify a job name:
#SBATCH -J raxml

# Specify an output file
#SBATCH -o raxml.out
#SBATCH -e raxml.out

# Run a command

module load raxml

raxmlHPC -f a -x 12345 -p 12345 -N 10 -m GTRGAMMA -s Nucleotide_alignment_BRPB.fasta -# 1000 -n Boot1000
#raxmlHPC -m GTRCAT -p 12345 -f b -t RAxML_bestTree.Boot100 -z RAxML_bootstrap.Boot100 -n Boot101 

~~~


Load Libraries 

~~~r
library(dada2); packageVersion("dada2")
library(ggplot2)
library(phyloseq)
source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite('phyloseq')
library("devtools")
devtools::install_github("benjjneb/dada2")
install.packages("ggplot2")
~~~

Check Mock Community Accuracy 

~~~r
unqs.mock <- seqtab.nochim["ZymoMock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) 
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path("BioPool_genomes/16S-18S/combine_v1.fasta"))
mock.ref
match.ref <- (sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")

~~~

Carry files over to QIIME II

1. Make Seqnece table

~~~r
write.table(t(seqtab.nochim), "seqtab-nochim.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
~~~

2. Make representative sequence

~~~r
uniquesToFasta(seqtab.nochim, fout='rep-seqs.fna', ids=colnames(seqtab.nochim))
~~~ 

3. Make representative contamination table 

~~~r
mock_contam <- sort(unqs.mock[! match.ref])
uniquesToFasta(mock_contam, fout='mock_contam.fna', ids=colnames(mock_contam))

~~~

#Import to qiime II

~~~Linux

#Rep-set sequneces 
qiime tools import  --input-path dada2-analysis/rep-seqs.fna --type 'FeatureData[Sequence]' --output-path rep-seqs.qza

#Add a special header for BIOM:

echo -n "#OTU Table" | cat - dada2-analysis/seqtab-nochim.txt > dada2-analysis/biom-table.txt

#Convert to BIOM v2.1:

biom convert -i dada2-analysis/biom-table.txt -o dada2-analysis/table.biom --table-type="OTU table" --to-hdf5

#Now we can import that as well:
qiime tools import --input-path table.biom --type 'FeatureTable[Frequency]' --source-format BIOMV210Format --output-path table.qza

#Import mock contaminants 
qiime tools import  --input-path mock_contam.fna --type 'FeatureData[Sequence]' --output-path mock_contam.qza
~~~

#Make Tree

~~~Linux

qiime alignment mafft --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza

iime alignment mask --i-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza

qiime phylogeny fasttree --i-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

~~~

#Evaluate Seqeunces 

~~~
qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file mappingfile_2.txt 
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
~~~

##Assign Taxonomy 
~~~
qiime tools import --type 'FeatureData[Sequence]' --input-path /Users/biancabrown/Desktop/all_small_mammal/gg_13_8_otus/rep_set/97_otus.fasta --output-path 97_otus.qza
  
qiime tools import --type 'FeatureData[Taxonomy]' --source-format HeaderlessTSVTaxonomyFormat --input-path /Users/biancabrown/Desktop/all_small_mammal/gg_13_8_otus/taxonomy/97_otu_taxonomy.txt --output-path ref-taxonomy.qza

GTGCCAGCMGCCGCGGTAA
GGACTACHVGGGTWTCTAAT

qiime feature-classifier extract-reads --i-sequences 97_otus.qza --p-f-primer GTGCCAGCMGCCGCGGTAA --p-r-primer GGACTACHVGGGTWTCTAAT --p-trunc-len 253 --o-reads ref-seqs.qza

qiime feature-classifier fit-classifier-naive-bayes --i-reference-reads ref-seqs.qza --i-reference-taxonomy ref-taxonomy.qza --o-classifier classifier.qza

qiime feature-classifier classify-sklearn --i-classifier classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

qiime taxa barplot --i-table table.qza --i-taxonomy taxonomy.qza --m-metadata-file mappingfile_2.txt --o-visualization taxa-bar-plots.qzv
~~~

##Contamination
~~~
#Filter Zymo experimental mock
  
qiime feature-table filter-samples --i-table table.qza --m-metadata-file mappingfile_2.txt --p-where '"Species" IN ("BLANK")' --o-filtered-table mock.qza

#Summarize experimental mock 
qiime feature-table summarize --i-table mock.qza --o-visualization mock.qzv

#Dada2 on contaminated sequences didn't work. Need to return to evaluate.

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path /Users/biancabrown/Desktop/all_small_mammal/mock_ref --source-format CasavaOneEightSingleLanePerSampleDirFmt --output-path mock_demux-paired-end.qza

qiime dada2 denoise-paired --i-demultiplexed-seqs mock_demux-paired-end.qza --o-table mock_table  --p-trunc-len-f 150 --p-trunc-len-r 150 --o-representative-sequences mock_rep-seq 

~~~
###Check Contamination in Dada2
~~~r
unqs.mock <- seqtab.nochim["ZymoMock",]
unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
mock.ref <- getSequences(file.path("BioPool_genomes/16S-18S/combine_v1.fasta"))
match.ref <- (sapply(names(unqs.mock), function(x) any(grepl(x, mock.ref))))
cat("Of those,", sum(match.ref), "were exact matches to the expected reference sequences.\n")
mock_contam <- sort(unqs.mock[! match.ref])

#Write mock contamination into a txt file
write.table(t(mock_contam), "mock_contam.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

~~~
###Prepare Mock contamiantion for qiime2
~~~
echo -n "#OTU Table" | cat - mock_contam.txt > mock_biom-table.txt

#Transpose mock_biom-table in excel 

biom convert -i mock_biom-table.txt -o mock_biom-table.biom --table-type="OTU table" --to-hdf5


~~~

###Filter contamination in qiime
~~~

qiime tools import --input-path mock_biom-table.biom --type 'FeatureTable[Frequency]' --source-format BIOMV210Format --output-path mock_biom-table.qza

qiime feature-table summarize  --i-table mock.qza --o-visualization mock.qzv


#Download CSV file for contamination sequences and then convert the csv file.

echo 'Feature ID\tFrequency' | cat - sample-frequency-detail.csv | tr "," "\\t" > sample-frequency-detail.tsv

#Filter Contaminants from experimental mock
qiime feature-table filter-features --i-table mock.qza --m-metadata-file sample-frequency-detail.tsv  --p-exclude-ids --o-filtered-table filtered-table.qza

#View expeirmental mock without contamination

qiime taxa barplot --i-table filtered-table.qza --m-metadata-file /Users/biancabrown/Desktop/all_small_mammal/mappingfile_2.txt --i-taxonomy /Users/biancabrown/Desktop/all_small_mammal/taxonomy.qza  --o-visualization mock_without_contam_plots.qzv

#View expeirmental mock with contamination
 qiime taxa barplot --i-table mock.qza --m-metadata-file /Users/biancabrown/Desktop/all_small_mammal/mappingfile_2.txt --i-taxonomy /Users/biancabrown/Desktop/all_small_mammal/taxonomy.qza  --o-visualization mock_with_contam_plots.qzv
 
 #View expeirmental mock contamination only

 qiime taxa barplot --i-table mock_biom-table.qza --m-metadata-file /Users/biancabrown/Desktop/all_small_mammal/mappingfile_2.txt --i-taxonomy /Users/biancabrown/Desktop/all_small_mammal/taxonomy.qza  --o-visualization mock_contam_only_plots.qzv

~~~
###Contamination in entire sequence table 
~~~
qiime feature-table filter-features --i-table table.qza --m-metadata-file /Users/biancabrown/Desktop/all_small_mammal/mock_ref/sample-frequency-detail.tsv  --p-exclude-ids --o-filtered-table full-table_contam_removed.qza

qiime taxa barplot --i-table full-table_contam_removed.qza --m-metadata-file /Users/biancabrown/Desktop/all_small_mammal/mappingfile_2.txt --i-taxonomy /Users/biancabrown/Desktop/all_small_mammal/taxonomy.qza  --o-visualization full-table_contam_removed.qza.qzv

#Remove chloroplast and mitochondria 
qiime taxa filter-table --i-table full-table_contam_removed.qza --i-taxonomy /Users/biancabrown/Desktop/all_small_mammal/taxonomy.qza --p-exclude mitochondria,chloroplast --o-filtered-table table-no-mitochondria-no-chloroplast.qza


##DO I need to remake tree?

##Qiime suggestion 
##You will not need to remove any samples or recreate the tree. The tree is made from only the representative sequences, agnostic of the samples present. As long as the analysis you are performing is based on sequences that are a subset of the sequences used to make the tree (or all of the sequences) you can proceed with whatever analysis. One issue that can arise is when samples are added to an analysis, after the tree was created, which contain sequences that are not present in the tree. In this situation a new tree would need to be created with all of the sequences.
https://forum.qiime2.org/t/rooted-tree-for-filtered-feature-table/658


~~~

Taking data back to R so I can perform phyloseq, vegan, dada2 etc. 

~~~R
qiime tools export table-no-mitochondria-no-chloroplast.qza --output-dir qiime_out/

qiime tools export taxonomy.qza --output-dir qiime_out/
~~~

##Beta Jack Knife Bray Curtis

~~~
qiime feature-table filter-samples --i-table table-no-mitochondria-no-chloroplast.qza --m-metadata-file mappingfile_2.txt --p-where "Species NOT IN ('RHY', 'Porcupine', 'BLANK', 'TAHA', 'ARNI', 'MUSO')" --o-filtered-table table-no-mitochondria-no-chloroplast-sans-samples.qza

qiime diversity beta-rarefaction --i-table table-no-mitochondria-no-chloroplast-sans-samples.qza --m-metadata-file mappingfile_2.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  table-no-mitochondria-no-chloroplast-sans-samples_10000.qzv

qiime tools extract table-no-mitochondria-no-chloroplast-sans-samples_10000.qzv --output-dir qiime_out/
~~~

##Beta Jack Knife Bray Curtis Sum vs. Mean vs. Median

~~~
qiime feature-table group --i-table table-no-mitochondria-no-chloroplast-sans-samples.qza --p-axis sample --m-metadata-column Species --m-metadata-file mappingfile_2.txt --p-mode sum --o-grouped-table sum_table-no-mitochondria-no-chloroplast-sans-samples.qza

qiime diversity beta-rarefaction --i-table sum_table-no-mitochondria-no-chloroplast-sans-samples.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  sum_table-no-mitochondria-no-chloroplast-sans-samples_10000.qzv

qiime feature-table group --i-table table-no-mitochondria-no-chloroplast-sans-samples.qza --p-axis sample --m-metadata-column Species --m-metadata-file mappingfile_2.txt --p-mode mean-ceiling --o-grouped-table mean_table-no-mitochondria-no-chloroplast-sans-samples.qza

qiime diversity beta-rarefaction --i-table mean_table-no-mitochondria-no-chloroplast-sans-samples.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  mean_table-no-mitochondria-no-chloroplast-sans-samples_10000.qzv

qiime feature-table group --i-table table-no-mitochondria-no-chloroplast-sans-samples.qza --p-axis sample --m-metadata-column Species --m-metadata-file mappingfile_2.txt --p-mode median-ceiling --o-grouped-table median_table-no-mitochondria-no-chloroplast-sans-samples.qza

qiime diversity beta-rarefaction --i-table median_table-no-mitochondria-no-chloroplast-sans-samples.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  median_table-no-mitochondria-no-chloroplast-sans-samples_10000.qzv



~~~

#notes on summing,mean, median qiime
~~~
Could you provide an example in which it would be better to use the mean and median ceiling? What benefit would that provide over summing?

To be honest, I’m not certain. Summing seems to make the most sense to me for standard uses, since the data will be normalized downstream anyway (e.g., even rarefaction prior to alpha diversity analysis).

I suppose if you suspected that certain samples within a group are outliers that would skew the mean, you might want to use median-ceiling instead to help control that.

Perhaps mean/median would also make the most sense if the frequencies will be used directly without normalization. For example, if a feature table contained quantitatively informative frequencies (e.g., metabolite concentrations or microbe abundances calculated using an absolutely quantitative method). Then summing those features would not make sense if you want to then look at the mean/median abundance of features in each sample type. But this would not really concern most current uses with microbiome sequencing data.

~~~
In addition to doing Jackknife in Qiime2, I performed a bray curtis analysis in phyloseq 

~~~
library(vegan)
mammals_bray <- phyloseq::distance(physeq_pruned, method = "bray")
mammals_bray
colorScale <- rainbow(length(levels(get_variable(physeq_pruned, "Species"))))
cols <- colorScale[get_variable(physeq_pruned, "Species")]
tip.labels <- as(get_variable(physeq_pruned, "Species"), "character")
tip.labels
meta<-sample_data(physeq_pruned)
meanDist <- meandist(mammals_bray, meta$Speices)
~~~

```{r}
taxa<-read.csv("taxa.csv", sep=",",row.names=1)
taxa<-as.matrix(taxa)
```
```{r}
trefile <- "rootedtree.nwk"
mapfile <- "mappingfile_2.txt"
map <-import_qiime_sample_data(mapfile)
otu<-otu_table(otu_table, taxa_are_rows=TRUE)
tax<-tax_table(taxa)
phy_tree<-read_tree(trefile)
physeq<-phyloseq(tax,otu,map,phy_tree)
physeq
physeq_pruned<- subset_samples(physeq, Species=ARNI)

physeq_pruned = subset_samples(physeq, Species != "ARNI" & Species != "Porcupine" & Species !="BLANK" & Species !="RHY" & Species !="TAHA"& Species !="MUSO")
physeq_pruned
```
#Back to Qiime2. There are outlirs that potential skewing the dataset 
####Remove  SampleID != "LM0035" & SampleID != "TK0000680"
~~~
qiime feature-table filter-samples --i-table table-no-mitochondria-no-chloroplast-sans-samples.qza --m-metadata-file mappingfile_2.txt --p-where "SampleID NOT IN ('LM0035', 'TK0000680')" --o-filtered-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza

qiime feature-table filter-samples --i-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mappingfile_2.txt --p-where "Species NOT IN ('GRMI')" --o-filtered-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza

qiime feature-table summarize  --i-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --o-visualization table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qzv

qiime diversity beta-rarefaction --i-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mappingfile_2.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv
~~~

 SampleID != "LM0035" & SampleID != "TK0000680"
 
 ##Jackknife II
 
 ~~~
 qiime feature-table group --i-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --p-axis sample --m-metadata-column Species --m-metadata-file mappingfile_2.txt --p-mode sum --o-grouped-table sum_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza
 


qiime diversity beta-rarefaction --i-table sum_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  sum_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv

qiime feature-table group --i-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --p-axis sample --m-metadata-column Species --m-metadata-file mappingfile_2.txt --p-mode mean-ceiling --o-grouped-table mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza

qiime diversity beta-rarefaction --i-table mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv

 ~~~
 
 ##Distance Matrix
 ~~~
 qiime diversity core-metrics-phylogenetic --i-phylogeny rooted-tree.qza --i-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --p-sampling-depth 5000 --m-metadata-file mappingfile_2.txt --output-dir core-metrics-results
 ~~~
 
 ~~~
 qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza  --m-metadata-file mappingfile_2.txt --m-metadata-column Plot_nonum --o-visualization core-metrics-results/bray-location-significance.qzv 
 ~~~
 ~~~
 
qiime taxa collapse --i-table mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza  --i-taxonomy taxonomy.qza  --p-level 7 --o-collapsed-table 7_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza

qiime diversity beta-rarefaction --i-table 7_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  7_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv


qiime taxa collapse --i-table mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza  --i-taxonomy taxonomy.qza  --p-level 6 --o-collapsed-table 6_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza
qiime diversity beta-rarefaction --i-table 6_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  6_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv

qiime taxa collapse --i-table mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza  --i-taxonomy taxonomy.qza  --p-level 5 --o-collapsed-table 5_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza
qiime diversity beta-rarefaction --i-table 5_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  5_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv


qiime taxa collapse --i-table mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza  --i-taxonomy taxonomy.qza  --p-level 4 --o-collapsed-table 4_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza
qiime diversity beta-rarefaction --i-table 4_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  4_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv

qiime taxa collapse --i-table mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza  --i-taxonomy taxonomy.qza  --p-level 3 --o-collapsed-table 3_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza
qiime diversity beta-rarefaction --i-table 3_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  3_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv

qiime taxa collapse --i-table mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza  --i-taxonomy taxonomy.qza  --p-level 2 --o-collapsed-table 2_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza
qiime diversity beta-rarefaction --i-table 2_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  2_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv

qiime taxa collapse --i-table mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza  --i-taxonomy taxonomy.qza  --p-level 1 --o-collapsed-table 1_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza

qiime diversity beta-rarefaction --i-table 1_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --m-metadata-file mapping_id_species.txt --p-metric braycurtis --p-clustering-method upgma --i-phylogeny rooted-tree.qza --p-sampling-depth 10000 --o-visualization  1_mean_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers_10000.qzv
~~~

~~~

qiime feature-table rarefy --i-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --p-sampling-depth 10000 --o-rarefied-table rarefied-table_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza

qiime diversity beta --i-table rarefied-table_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --p-metric braycurtis  --o-distance-matrix bray_distance-matrix_rarefied-table_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza

qiime tools export bray_distance-matrix_rarefied-table_table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --output-dir bray_distance-matrix_rarefied
~~~

with GRMI

~~~
qiime feature-table filter-samples --i-table table-no-mitochondria-no-chloroplast-sans-samples.qza --m-metadata-file mappingfile_2.txt --p-where "SampleID NOT IN ('LM0035', 'TK0000680')" --o-filtered-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza

qiime feature-table rarefy --i-table table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --p-sampling-depth 10000 --o-rarefied-table rarefied-table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza


qiime diversity beta --i-table rarefied-table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --p-metric braycurtis  --o-distance-matrix bray_distance-matrix_rarefied-table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza

qiime tools export bray_distance-matrix_rarefied-table-no-mitochondria-no-chloroplast-sans-samples-sansoutliers.qza --output-dir bray_distance-matrix_rarefied_2
~~~








