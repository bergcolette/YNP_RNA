# Generating a pseudoreference for aligning RNAseq reads
_originally written by Andrew Demaree in 2020_

## Part 1: Pseudo-reference genome construction and competitive transcriptome alignment
- Need to have a reference to align to that doesn’t bias one species over another 
- The goal is to construct a  diploid pseudoreference genome from two genotypes and competitively map all RNAseq reads against it
- this tutorial assumes you have the packages Trimmomatic, bwa, picard, samtools, GATK4, and STAR installed


### BWA-MEM Align trimmed reads to the reference genome.
https://github.com/lh3/bwa
- BWA is used to align your reads (in fastq format) to the reference genome. These commands requires all index files to be present in the same directory as the reference
- If reference index files are not present, run the following command:
  
```bash
bwa index reference.fa
bwa mem -t 8 <file/path/to/reference.fa> <file/path/to/trimmed_forward.1.paired.fastq> <file/path/to/trimmed_reverse.2.paired.fastq> > <file/path/to/output/aligned.sam>
```

### Samtools Sort sort aligned sam file.
http://www.htslib.org/doc/samtools-sort.html

```bash
samtools sort -T <name_of_file> -o <file/path/to/output/name_of_file.sort.sam> <file/path/to/aligned.sam>
```

###	AddOrReplaceReadgroups adds read groups to sorted sam file.
https://broadinstitute.github.io/picard/

```java
java -jar /file/path/to/java_jar/picard.jar AddOrReplaceReadGroups \ 
INPUT=<file/path/to/output/name_of_file.sort.sam> OUTPUT=<file/path/to/output/name_of_file.RG.sort.sam> \ 
RGSM=<name_of_file> RGLB=<library_type(e.g.TruSeq)> RGPL=<sequencer_type(e.g.Illumina)> RGPU=UNKNOWN VALIDATION_STRINGENCY=LENIENT
```

###	Picard Mark Duplicates Removes optical and  PCR duplicates from thealignment. 
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
- Requires coordinate sorted bam or sam as input 
- Don’t always specify removal of sequencing duplicates, however in this instance the sequencing_duplicates tag is used to remove optical duplicates 
- The remove_duplicates tag uses the Library (LB) read group (@RG) field to remove duplicates from the same library sequenced on different lanes
  
```java 
java -jar /file/path/to/java_jar/picard.jar MarkDuplicates \ 
INPUT=<file/path/name_of_file.RG.sort.sam> OUTPUT=<file/path/name_of_file.rmdup.RG.sort.sam> \ 
METRICS_FILE=<name_of_file.rmdup_metrics_fix \ 
VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE
```

## Samtools View/Index Filter out reads with a quality score below 20 from the SF5 alignment, convert it to a bam file, and create an index file.
http://www.htslib.org/doc/samtools-view.htm
- The bam file format is required for haplotype calling 
  
```bash
samtools view -bS  <name_of_file.sam> -q 20 > <name_of_file.bam>
samtools index <name_of_file.bam>
```

### Generating SNPs
- Goal: generate a set of high-quality SNPs for pseudoreference genome construction.
- GATK HaplotypeCaller (in GVCF mode) Identify phased SNP and insertion/deletion variants from the SF5 alignment.
https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
- Single sample GVCF calling without allele specific annotation was used because biallelic sites will be extracted and used, so multiallelic site calling doesn’t really matter

```java
gatk --java-options "-Xmx4g" HaplotypeCaller  \ 
-R <file/path/to/reference.fasta> \ 
-I  <file/path/to/name_of_file.bam> \ 
-O <file/path/to/output.g.vcf.gz> \ -ERC GVCF
```

###	GATK GenotypeGVCFs 
https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs

```java
gatk --java-options "-Xmx4g" GenotypeGVCFs \ 
-R <file/path/to/reference.fasta> \  
-V <file/path/to/output.g.vcf.gz> \  
-O <file/path/to/output.vcf.gz>
```

###	GATK SelectVariants Extract only the biallelic SNPs from the vcf file containing the genotype calls for the SF5 alignment.
https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants

```java
gatk SelectVariants \ 
-R <file/path/to/reference.fasta> \ 
-V <file/path/to/output.vcf.gz> \ 
--select-type-to-include SNP \ 
--restrict-alleles-to BIALLELIC  
-O <file/path/to/biallelic_snps.output.vcf>
```

###	GATK VariantFiltration Mark variants with “QUAL” in the filter column if the quality score (phred-scaled confidence that a variant is present) is less than 40 and “DEPTH” if the read depth in support of a variant site is less than 2.
https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration
https://commons.apache.org/proper/commons-jexl/reference/syntax.html

```java
gatk VariantFiltration \ 
-R <file/path/to/reference.fasta> \ 
-V <file/path/to/biallelic_snps.output.vcf> \ 
-O <file/path/to/biallelic_snps.filtered.output.vcf> \ 
--filterExpression "QUAL &lt; 40.0" --filterName "QUAL" \ 
--filterExpression "DP &lt; 2" --filterName "DEPTH"
```

###	GATK SelectVariants Remove all variants with anything other than ‘.’ or ‘PASS’ in the Filter field, in effect filtering out sites with mapping quality below 40 or read depth below two.
https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants

```java
gatk SelectVariants \ 
-R <file/path/to/reference.fasta> \ 
-V <file/path/to/biallelic_snps.filtered.output.vcf> \ 
--select-type-to-include SNP \ 
--exclude-filtered=TRUE  -O <file/path/to/biallelic_snps.filtered.dropped.output.vcf>
```

###	GATK FastaAlternateReferenceMaker Use filtered SNPs to generate a M. Nasutus SF5 pseudoreference genome.
https://gatk.broadinstitute.org/hc/en-us/articles/360037594571-FastaAlternateReferenceMaker#--snp-mask
- Don’t want to use masking file because it replaces SNPs with Ns
- Don’t need to provide intervals because we want SNPs across the whole genome

```java
gatk FastaAlternateReferenceMaker \ 
-R <file/path/to/reference.fasta> \ 
-O <reference_pseudo.fasta> \ 
-V <file/path/to/biallelic_snps.filtered.dropped.output.vcf> 
```

###	 Appending Allelic Identifiers: Add allelic identifiers (IM62 and SF5) to the chromosome names in the M. guttatus v2.0 reference and M. nasutus SF5 pseudoreference fasta files.
- constructed a diploid M. guttatus IM62-M. nasutus SF5 pseudoreference genome by appending allelic (i.e. IM62 or SF5) identifiers onto the chromosome names in the M. guttatus v2.0 reference and M. nasutus SF5 pseudoreference fasta files
- M. guttatus allelic reference:
```bash
sed '/>/ s/$/_Allelic_ID/' <file/path/to/reference.fasta> > <file/path/to/reference_Allelic_ID.fasta>
```
- M. nasutus SF5 allelic pseudo reference:
```bash
sed '/>/ s/$/_SF5/' <file/path/to/reference_pseudo.fasta> > <file/path/to/reference_pseudo_Allelic_ID.fasta>
```

###	Creating Diploid genome File: Manually merge the M.guttatus allelic reference and the M. nasutus SF5 allelic pseudo reference.
```bash 
cat <file/path/to/reference_Allelic_ID.fasta> <file/path/to/reference_pseudo_Allelic_ID.fasta> > <file/path/to/reference_Allelic_ID_reference_pseudo_Allelic_ID.fasta
```

###	Append Allelic Identifiers in GFF3 Files Create a diploid annotation file by appending allelic identifiers onto chromosome, gene, and transcript names, then manually combine the two files and convert the resulting GFF3 file to a GFT file. 
- The script to append allelic identifiers is called append_identifiers_GFF3.py and its file path on carnation is home/andrew_demaree/scripts/append_identifiers_GFF3.py. A copy has also been added to the RNA2020 Dropbox.The output of this scriptis two GFF3 files (SF5_Mguttatus_v2.0_256_gene_exons.gff3 and IM62_Mguttatus_v2.0_256_gene_exons.gff3)
- To change the input file, output file, allelic identifiers, or matched attributes, open the script using VIM (or any text editor that won’t alter the encoding) and change the appropriate string.
- After running the script, append the two output files together in the same order as the reference files.
```bash
cat <file/path/to/Allelic_ID_1.gff3> <file/path/to/Allelic_ID_2.gff3> > <file/path/to/output/Allelic_ID_1_2.gff3>
```

###	Gffread Convert the merged gff3 file to a gft file using gffread.
https://github.com/gpertea/gffread
`gffread my.gff3 -T -o my.gtf`

## Part 2: RNAseq Alignment Map IM62 and SF5 genotypes to the diploid pseudoreference genome.
- Trimmomatic Trim adapter sequences and low-quality bases from the raw RNAseq reads, then filter out reads shorter than 36 bp.
http://usadellab.org/cms/?page=trimmomatic

```java
java -jar /file/path/to/trimmomatic-0.35.jar SE \ 
-threads 8 -phred33 <file/path/to/forward.1.fastq> <file/path/to/reverse.2.fastq> \ 
<file/path/to/output_dir/trimmed_forward.1.paired.fastq> <file/path/to/output_dir/trimmed_forward.1.unpaired.fastq> \ 
<file/path/to/output_dir/trimmed_reverse.2.paired.fastq> <file/path/to/output_dir/trimmed_reverse.2.unpaired.fastq> \ 
ILLUMINACLIP:file/path/to/adapter/(illumina/Many.TruSeq.PE.fa):2:20:10:4 \ 
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### STAR Align the 36-50-bp Single-end RNAseq reads to the diploid pseudoreference generated in part 1.
https://github.com/alexdobin/STAR
- Generating Genome index files:
- --sjdbOverhang: get the read length from the trimmed read fastqc and subtract 1 from it.
- genomeDir is the directory where the genome indices are kept. Need to make this directory before running and should be empty.
- --sjdbGTFfile is the file path to the diploid GTF3 file generated for the pseudoreference.

```bash
~/STAR/STAR/bin/Linux_x86_64_static/STAR \
--runThreadN 8 --runMode genomeGenerate \
--genomeDir /path/to/genomeDir \
--genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 … \
--sjdbGTFfile /path/to/annotations.gtf --sjdbOverhang ReadLength-1
```

### Mapping Reads:
-	-- readFilesIn can take multiple samples in the form sample1,sample2,sample3,...,sampleN
-	Reads that overlap a SNP are expected to map uniquely to an SF5 or IM62 allele in the pseudoreference genome, while reads that do not overlap a SNP will map equally well to both alleles
-	limit each read to a maximum of two alignments by specifying --outFilterMultimapNmax 2
-	Randomly designate one of the two allowed alignments as primary by specifying --outMultimapperOrder Random

```bash 
	~/STAR/STAR/bin/Linux_x86_64_static/STAR \
    --runThreadN NumberOfThreads \
    --genomeDir /path/to/genomeDir \
    --readFilesIn /path/to/read1 /path/to/read2 … /path/to/readN \
    --genomeDir <file/path/to/genomeDir> --outFilterMultimapNmax 2 \
    --outMultimapperOrder Random
```
```bash
for i in *.P.fastq; 
do ~/STAR/STAR/bin/Linux_x86_64_static/STAR \
--runThreadN 8 \
--genomeDir GenomeDir/ \
--readFilesIn $i \
--outFilterMultimapNmax 2 \
--outMultimapperOrder Random \
--outFileNamePrefix ${i:0:-8}; 
done
```

### SAMtools view Remove secondary alignments and unmapped reads, leaving just unique (i.e. allele-specific) and primary alignments.
-	-F tells samtools to omit sequences with the int value that follows the tag (e.g. 4 = reads that are unmapped, 100 = reads that are not primary)
```bash
samtools view -b -F 4 -F 100 <file/path/to/aligned.bam> > <file/path/to/aligned.filtered.bam>
```
```bash 
for i in *.sam; do samtools view -b -F 4 -F 100 $i > ${i:0:-15}.filtered.bam; done
```

### SAMtools sort Convert the filtered bams to coordinates sorted bams.
```bash
samtools sort <file/path/to/aligned.filtered.bam> >  <file/path/to/aligned.filtered.sort.bam>
```
```bash
for i in *filtered.bam; do samtools sort $i > ${i:0:-4}.sort.bam; done
```

### Picard Mark Duplicates Remove optical and  PCR duplicates from the RNAseq alignments. 
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
- Requires coordinate sorted bam or sam as input 

```java
java -jar <file/path/to/java/jar/picard.jar MarkDuplicates \ 
INPUT=<file/path/name_of_coordinate_sorted.bam> OUTPUT=<file/path/name_of_coordinate_sorted.rmdup.bam> \ METRICS_FILE=<name_of_file.rmdup_metrics_fix \ 
VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE
```

```bash
for i in *.sort.bam; do java -jar /home/thom_nelson/opt/picard.jar MarkDuplicates INPUT=$i OUTPUT=${i:0:-4}.rmdup.filtered.bam METRICS_FILE=SF5.rmdup_metrics_fix VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE; done
```

### Picard Create Sequence Dictionary Creates a .dict file required by future tools.
https://gatk.broadinstitute.org/hc/en-us/articles/360036729911-CreateSequenceDictionary-Picard-
```java 
java -jar picard.jar CreateSequenceDictionary \ R=<file/path/name_of_pseudoreference.fa
```

###	Samtools faidx Creates a .fai index file required by future tools. 
http://www.htslib.org/doc/samtools-faidx.html
```bash
samtools faidx <file/path/to/name_of_pseudoreference.fa>
```

### GATK SplitNCigarReads Parse intron-spanning reads into exon segments and trim off bases that extend into intronic regions.
- Splits reads that contain Ns in their Cigar String into K+1 sequences, where k is the number of Ns.

```java
gatk SplitNCigarReads \ 
-R <file/path/to/Mguttatus_V2refmtcp_Mnasutus_pseudo.fasta \ 
-I <file/path/name_of_coordinate_sorted.rmdup.bam> \ 
-O <file/path/name_of_coordinate_sorted.rmdup.split.bam>
```

```bash
for i in *rmdup.filtered.bam; do ~/gatk-4.1.7.0/gatk SplitNCigarReads -R pseudoreference/MguttatusV2refmtcp_IM62_SF5_pseudo.fa -I $i -O ${i:0:-4}.split.bam; done
```

### SAMtools view filter alignments based on mapping quality.
- Obtain allele-specific and primary mapping (i.e. total) reads using a Q20 threshold.
```bash
samtools view -b <file/path/name_of_coordinate_sorted.rmdup.split.bam> -q 20 > <file/path/name_of_coordinate_sorted.rmdup.split.Q20.bam>
```
```bash
for i in *.split.bam; do samtools view -b -q 20 -o ${i:0:-4}.q20.bam $i;done
```

###	Htseq-count 
https://htseq.readthedocs.io/en/release_0.11.1/index.html
- I have created a python virtual environment called htseq_new which has Htseq-count installed otherwise follow the directions in the htseq documentation.
```bash
source activate htseq_new

htseq-count -f bam --nonunique=all --stranded=reverse <file/path/to/sample.split.q20.bam <file/path/to/reference_genes_exons.gtf>

source deactivate

source activate htseq_new

for i in *.split.q20.bam; do htseq-count -f bam --nonunique=all --stranded=reverse $i ~/RNASeq/trsht/IM62_SF5_Mguttatus_v2.0_256_gene_exons.gtf > ${i:0:-40}.counts.tsv;done

Source deactivate
```


Useful Commands:
Count read totals: 
```bash 
cat SF5_carpel_1.P.fastq | echo $((`wc -l`/4))
```
Count primary mapped reads: 
```bash
samtools view -F 4 -F 256 -F 1024 SF5_carpel_1.filtered.bam | wc -l
```

