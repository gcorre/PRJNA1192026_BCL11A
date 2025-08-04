# PRJNA1192026_BCL11A
Bioinformatics material for the analysis of ONT LongReads

This repository contains scripts and intermediate files used to analyse the raw data under SRA id PRJNA1192026.
The 4 libraries contains long reads sequences from targeted Cas9 enrichement at BCL11A loci.

# Reads generation
Basecalling was performed using the HAC model (High Accuracy).

# Reads processing 
Reads were analyzed using the snakemake workflow 'ONT.crisprV5.snakemake' located in 01-Reads_Processing'.
For each library, the YAML configuration file was placed in the folder containing the fastq file.

The Tools.YAML file contains the dependencies used to create the conda environment containing the required programs.

The workflow generates bam, bed, vcf files used to perform the post-processing.
mpileup was used to extract nucleotide composition at loci position (5nt window).

~~~
# ================================================================
# Mpileup file generation for the 2 loci from the bam file
# ================================================================

 for lib in *; do 
   samtools mpileup $lib/locus/$lib.locusAligned.sorted.end2end.bam \
   -r BCL11A:8490-8495 \
   --output-QNAME  \
   --no-output-ins \
   --no-output-del \
   --no-output-ins-mods \
   --no-output-ends \
   --min-BQ 0 \
   -aa \
   -o $lib"_ATF4.mpileup" \ 
   echo $lib;\
 done


 for lib in *; do \
   samtools mpileup $lib/locus/$lib.locusAligned.sorted.end2end.bam \
   -r BCL11A:5264-5270 
   --output-QNAME  \
   --no-output-ins \
   --no-output-del \
   --no-output-ins-mods \
   --no-output-ends \
   --min-BQ 0 \
   -aa \
   -o $lib"_GATA1.mpileup"; \
   echo $lib;\
 done
~~~
Folders in 01-Reads_Processing contains VCF, BED, mpileup and composition intermediate files for the post-processing step below.


# Post-Processing
Folder 02-Post_Processing contains 2 R scripts : 
 - allele_from_pileup.R
Perform analysis of the mpileup files to classify loci as modified or not and also perform a co-editing analysis for each single read

 - base_composition_logo.R
Perform the SV analysis from bed and VCF files to quantify INDELs size and frequency. Also generate seqLogo figures.

Sub-folders contains intermediate files used to generate graphs.
