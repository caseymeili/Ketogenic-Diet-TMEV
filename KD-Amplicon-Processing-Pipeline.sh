#!/bin/bash

# 16S Data Processing (V3-V4)
# mothur v.1.48.3
# adapted from MiSeq SOP (accessed August 2025), https://mothur.org/wiki/miseq_sop/

####################################################################################

# launch mothur, set working directory, and set processors
ml mothur/1.48.3
set.dir(input=/Users/caseymeili/Desktop/28160R/Fastq/KD)
set.current(processors=12)

# make stability file (file containing gz files in directory)
make.file(inputdir=/Users/caseymeili/Desktop/27023R/Fastq/MS, type=gz, prefix=stability) 

# combine two sets of reads for each sample, and then combine the data from all of the samples
make.contigs(file=stability.files)
summary.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table) 

# screen to remove sequences longer than 460 bp, shorter than 400 bp, with ambiguous bases, or more than 8 homopolymers
screen.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table, maxambig=0, minlength=400, maxlength=460, maxhomop=8) 
summary.seqs(fasta=stability.trim.contigs.good.fasta, count=stability.contigs.good.count_table)

# merge duplicated sequences
unique.seqs(fasta=stability.trim.contigs.good.fasta, count=stability.contigs.good.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.fasta, count=stability.trim.contigs.good.count_table)

# create database customized to the region of interest (V3-V4)
# SILVA database (silva.nr_v138_2.align & silva.nr_v138_2.tax) obtained from mothur silva reference files, release 138.2

# determine where primers sit on the alignment using oligos file
# pcr.seqs(fasta=silva.nr_v138_2.align, taxonomy=silva.nr_v138_2.tax, oligos=zymo_primers.oligos, keepdots=F)

# rename output files to friendlier names
# rename.file(input=silva.nr_v138_2.pcr.align, new=silva.V3V4.pcr.align)
# rename.file(input=silva.nr_v138_2.pcr.tax,   new=silva.V3V4.pcr.tax)

# align sequences to new V3–V4 SILVA reference
align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.V3V4.pcr.align)
summary.seqs(fasta=stability.trim.contigs.good.unique.align)

# remove sequences that align outside of expected positions (start=4965, end=21977), usually due to poor alignment or non-specific amplification
screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, start=4965, end=21977)                      

# remove redundancy that could have been created
unique.seqs(fasta=stability.trim.contigs.good.unique.good.align, count=stability.trim.contigs.good.good.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.unique.align, count=stability.trim.contigs.good.unique.good.count_table)

# precluster allowing for up to 2 differences between sequences
pre.cluster(fasta=stability.trim.contigs.good.unique.good.unique.align, count=stability.trim.contigs.good.unique.good.count_table, diffs=2)

# remove chimeras using VSEARCH algorithm
chimera.vsearch(fasta=stability.trim.contigs.good.unique.good.unique.precluster.align, count=stability.trim.contigs.good.unique.good.unique.precluster.count_table, dereplicate=t)

# classify sequences using reference files
# PDS reference files (trainset19_072023.pds.fasta & trainset19_072023.pds.tax) were obtained from mothur RDP reference files version 19 (https://mothur.org/wiki/rdp_reference_files/)
classify.seqs(fasta=stability.trim.contigs.good.unique.good.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.unique.precluster.denovo.vsearch.count_table, reference=trainset19_072023.pds.fasta, taxonomy=trainset19_072023.pds.tax)

# remove unwanted sequences (18S, archaea, chloroplasts, mitochondria, unknown)
remove.lineage(fasta=stability.trim.contigs.good.unique.good.unique.precluster.denovo.vsearch.fasta, count=stability.trim.contigs.good.unique.good.unique.precluster.denovo.vsearch.count_table, taxonomy=stability.trim.contigs.good.unique.good.unique.precluster.denovo.vsearch.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.tax(taxonomy=current, count=current)

# rename files before analysis (so the file names aren't so yucky)
rename.file(fasta=current, count=current, taxonomy=current, prefix=KD)

# clustering into OTUs
# default cutoff used for clustering (0.03)
dist.seqs(fasta=kd.fasta, cutoff=0.03)
cluster(column=kd.dist, count=kd.count_table)

# make shared file
make.shared(list=kd.opti_mcc.list, count=kd.count_table, label=0.03)

# consensus taxonomy for each OTU
# taxonomy file is the output of the classify.seqs command
# generates kd.opti_mcc.0.03.cons.taxonomy which is used for the taxon page of phyloseq object
classify.otu(list=kd.opti_mcc.list, count=kd.count_table, taxonomy=kd.taxonomy, label=0.03)

# get number of sequences in each sample
count.groups(shared=kd.opti_mcc.shared)

# generate subsampled file for analysis
# size of smallest group from count.groups was 4 sequencing, subsampling to 18670 sequencing (second smallest)
# one sample (28160X24) was discarded from subsequent analyses
sub.sample(shared=kd.opti_mcc.shared, size=18670)

# alpha diversity
# generate rarefaction curves for plotting in R
# output file: kd.opti_mcc.0.03.subsample.groups.rarefaction
rarefaction.single(shared=kd.opti_mcc.0.03.subsample.shared, calc=sobs)

# calculate coverage, number of seqs, species observed, and diversity indices (inverse Simpson, Shannon, pielou)
# be sure to subsample to the size of the smallest group for all alpha diversity analysis
# if the smallest group is too small, discard the sample
summary.single(shared=kd.opti_mcc.shared, calc=nseqs-coverage-sobs-invsimpson-shannon-shannoneven, subsample=18670)
