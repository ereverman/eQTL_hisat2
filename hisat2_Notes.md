# Notes on my process for figuring out how to apply hisat to DSPR RNAseq data.

Starting from my hisat_eQTL/dm6_HISAT2_index directory, after running the config.sh script.
I also already downloaded the UCSC genome fasta from http://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/
### File information: dm6.fa.gz (14-Oct-2014,10:55,43M)
This is the only file in this directory to start with

1. For troubleshooting, open an interactive node and create a conda environment:
```
srun --time=4:00:00 --ntasks=1 --nodes=1 --partition=sixhour --pty /bin/bash -l

module load anaconda
conda create -n HISAT2_env
# note that this may not make an environment that the final array script will be able to access. It will have a different name.

conda activate HISAT2_env
```

2. Generate index files that are specific to DSPR that are SNP aware:
```
# code from Stuart

# install hisat2 and samtools:
conda install -c bioconda -c conda-forge hisat2=2.1.0=py27pl5.22.0_0
conda install -c bioconda samtools

# unpack the genome:
gunzip dm6.fa.gz # -k flag doesn't seem to work. this overwrites the original file

# Check the chromosome sizes/contig sizes:
samtools faidx dm6.fa

# get flybase genome fasta file:
wget ftp://ftp.flybase.net/releases/FB2018_06/dmel_r6.25/fasta/dmel-all-chromosome-r6.25.fasta.gz

# unpack the flybase genome:
gunzip dmel-all-chromosome-r6.25.fasta.gz

# Check chromosome/contig sizes:
samtools faidx dmel-all-chromosome-r6.25.fasta

# Get flybase gtf file:
wget  ftp://ftp.flybase.org/genomes/dmel/current/gtf/dmel-all-r6.33.gtf.gz

# Unpack:
gunzip dmel-all-r6.33.gtf.gz



# format the gtf file:

# (1) Gather the 1st column (chr IDs) of the FlyBase GTF via
cut -f1 dmel-all-r6.33.gtf > dmel-all-r6.33_chrArm.txt

# (2) Count number of unique items in the column using
awk '{A[$1]++}END{for(i in A)print i,A[i]}' dmel-all-r6.33_chrArm.txt > dmel-all-r6.33_chrArmCounts.txt

# (3) Copy original "dmel-all-r6.21.gtf" to "dmel-all-r6.21-edit.gtf" using
cp dmel-all-r6.33.gtf dmel-all-r6.33-edit.gtf

# (4) Using information on the names of the unwanted annotations from
#     "dmel-all-r6.21_chrArmCounts.txt" use find/remove rows in BBEdit
#     to eliminate from "dmel-all-r6.21-edit.gtf" those annotations
#     I don't care about. Use the Counts.txt file for search terms. Remember
#     to forward and backward search.

# (5) Now that "dmel-all-r6.21-edit.gtf" only contains annotations from
#     the 8 desired fragments, convert chromosome IDs in this file, again
#     with BBEdit using find/replace with the 'grep' function on
#     (2R\s to chr2R\t, 3R\s to chr3R\t, 2L\s to chr2L\t, 3L\s to chr3L\t, X\s to chrX\t, Y\s to chrY\t,
#      4\sFlyBase to 4\tFlyBase, mitochondrion_genome\s to chrM\t)

# (6) To check the edited "dmel-all-r6.21-edit.gtf" still contains all
#     the relevant annotations as the original file run
cut -f1 dmel-all-r6.33-edit.gtf > dmel-all-r6.33-edit_chrArm.txt
awk '{A[$1]++}END{for(i in A)print i,A[i]}' dmel-all-r6.33-edit_chrArm.txt > dmel-all-r6.33-edit_chrArmCounts.txt

# (7) Compare the pair of *_chrArmCounts.txt files (from the original
#     and the edited FlyBase GTF) and ensure that (a) the edited version
#     only contains annotations from the 8 fragments we care about, and
#     (b) that for those fragments the numbers of elements are the same
#     between the edited and original file


# Extract splice site positions from edited file (.py script already exists with package)
hisat2_extract_splice_sites.py dmel-all-r6.33-edit.gtf > dmel-all-r6.33-splicesites.txt

# Extract exon positions from edited file:
hisat2_extract_exons.py dmel-all-r6.33-edit.gtf > dmel-all-r6.33-exons.txt

# fix additional compatibility issues:
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' dm6.fa > dm6_allcaps.fa

# use custom script from John Kelly to change "." for SNPs to "CHROM.POS"
gunzip DSPR.r6.SNPs.vcf.gz
python ../scripts/fixID.py DSPR.r6.SNPs.vcf # this creates a new file
mv x.DSPR.r6.SNPs.vcf DSPR.r6.SNPs_edit.vcf

# Run the hisat2 vcf conversion script. needs genome file, vcf filenames, and base file name (output)
hisat2_extract_snps_haplotypes_VCF.py --non-rs -v dm6_allcaps.fa DSPR.r6.SNPs_edit.vcf DSPR_dm6_var


# Generate the custom index:
# This got killed on the cluster so make a small script to submit.

nano hisat-build.sh

hisat2-build -f -p 8 --ss dmel-all-r6.33-splicesites.txt --exon dmel-all-r6.33-exons.txt --snp DSPR_dm6_var.snp --haplotype DSPR_dm6_var.haplotype dm6.fa dm6_annot_dsprsnp

```
I loaded a small number of files for testing the code.
