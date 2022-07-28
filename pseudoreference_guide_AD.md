RNASeq Guide Document

Paper: Rampant misexpression in a Mimulus (monkeyflower) introgression line caused by hybrid sterility, not regulatory divergence

RED: These lines represent generalized command line arguments that can be modified to accommodate differences in directory structures and input files.
BLUE: These lines represent the exact commands I used based on my carnation directory structure. 

Part 1: Pseudo-reference genome construction and competitive transcriptome alignment
a.	Need to have a reference to align to that doesn’t bias one species over another 
b.	Aligning SF5, FER and STE reads against the M. guttatus v2.0 reference will likely introduce mapping bias
c.	The goal is to construct a  diploid M. guttatus IM62-M/nasutus SF5 pseudoreference genome and competitively map all RNAseq reads against it

M. nasutus alignment to publicly available SF5 whole genome (gDNA) sequence data:
Goal: To construct a high quality bam file containing reads from SF5 whole genome sequencing data aligned to the M. guttatus v2.0 reference.
2.	Fastq-dump (NCBI toolkit)
 https://ncbi.github.io/sra-tools/fastq-dump.html
a.	Use fastq-dump to download publicly available sequencing data
b.	SF5 fastqs - SRR400478 (75 bp reads)
c.	Fastq-dump -O <path/to/output_dir/ -I --split-files SRA_accesion_number
d.	fastq-dump --split-files SRR400478
i.	Splits dumped fastqs into files containing .1 and .2 paired end reads

3.	Fastqc Check for Adapters, read quality, and quality scores before trimming.
a.	Examine the the html output files to look for read quality issues
b.	for i in *.fastq; do fastqc $i; done
c.	fastqc SRR*
4.	Trimmomatic
 http://usadellab.org/cms/?page=trimmomatic
a.	Trim adapters, low quality bases, and filter reads less than 50 bp
b.	Align 50 -75 bp trimmed reads to M. guttatus v2.0 reference genome using trimmomatic in paired-end mode 
c.	java -jar /file/path/to/trimmomatic-0.35.jar PE \ -threads 8 -phred33(quality encoding) <file/path/to/forward.1.fastq> <file/path/to/reverse.2.fastq> \ <file/path/to/output_dir/trimmed_forward.1.paired.fastq> <file/path/to/output_dir/trimmed_forward.1.unpaired.fastq> \ <file/path/to/output_dir/trimmed_reverse.2.paired.fastq> <file/path/to/output_dir/trimmed_reverse.2.unpaired.fastq> \ ILLUMINACLIP:/file/path/to/adapter/adapter.fa:2:20:10:4 \ LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50 
d.	java -jar /home/thom_nelson/opt/Trimmomatic-0.35/trimmomatic-0.35.jar PE -threads 8 -phred33 SRR400478_1.fastq SRR400478_2.fastq SF5_trimmed.1.P.fastq SF5_trimmed.1.U.fastq SF5_trimmed.2.P.fastq SF5_trimmed.2.U.fastq ILLUMINACLIP:/home/thom_nelson/resources/Illumina/Many.TruSeq.PE.fa:2:20:10:4 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

5.	Fastqc Check for Adapters, read quality, and quality scores for paired reads after trimming.
a.	Examine html output for read quality issues after trimming
b.	for i in *paired.fastq; do fastqc $i; done
c.	fastqc *.P.fastq 

6.	BWA-MEM Align trimmed reads to the M. guttatus v2.0 reference genome.
https://github.com/lh3/bwa
a.	Requires all index files to be present in the same directory as the reference
b.	If reference index files are not present run the following command:
i.	bwa index reference.fa
c.	bwa mem -t 8 <file/path/to/reference.fa> <file/path/to/trimmed_forward.1.paired.fastq> <file/path/to/trimmed_reverse.2.paired.fastq> > <file/path/to/output/aligned.sam>
d.	bwa mem -t 8 /home/thom_nelson/resources/Mimulus/guttatus/IM62_v2.0/Mguttatus_V2refmtcp.fa SF5_trimmed.1.P.fastq SF5_trimmed.2.P.fastq > SF5.sam

7.	Samtools Sort sort aligned sam file.
http://www.htslib.org/doc/samtools-sort.html
a.	samtools sort -T <name_of_file> -o <file/path/to/output/name_of_file.sort.sam> <file/path/to/aligned.sam>
b.	samtools sort -T SF5 -o SF5.sort.sam SF5.sam

8.	AddOrReplaceReadgroups adds read groups to sorted sam file.
https://broadinstitute.github.io/picard/
a.	java -jar /file/path/to/java_jar/picard.jar AddOrReplaceReadGroups \ INPUT=<file/path/to/output/name_of_file.sort.sam> OUTPUT=<file/path/to/output/name_of_file.RG.sort.sam> \ RGSM=<name_of_file> RGLB=<library_type(e.g.TruSeq)> RGPL=<sequencer_type(e.g.Illumina)> RGPU=UNKNOWN VALIDATION_STRINGENCY=LENIENT
b.	java -jar /home/thom_nelson/opt/picard.jar AddOrReplaceReadGroups INPUT=SF5.sort.sam OUTPUT=SF5.RG.sort.sam RGSM=SF5 RGLB=TruSeq RGPL=Illumina RGPU=UNKNOWN VALIDATION_STRINGENCY=LENIENT

9.	Picard Mark Duplicates Removes optical and  PCR duplicates from the SF5 alignment. 
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
a.	Requires coordinate sorted bam or sam as input 
b.	Don’t always specify removal of sequencing duplicates, however in this instance the sequencing_duplicates tag is used to remove optical duplicates 
c.	The remove_duplicates tag uses the Library (LB) read group (@RG) field to remove duplicates from the same library sequenced on different lanes
d.	java -jar /file/path/to/java_jar/picard.jar MarkDuplicates \ INPUT=<file/path/name_of_file.RG.sort.sam> OUTPUT=<file/path/name_of_file.rmdup.RG.sort.sam> \ METRICS_FILE=<name_of_file.rmdup_metrics_fix \ VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE
e.	java -jar /home/thom_nelson/opt/picard.jar MarkDuplicates INPUT=SF5.RG.sort.sam OUTPUT=SF5.rmdup.RG.sort.sam METRICS_FILE=SF5.rmdup_metrics_fix VALIDATION_ST
RINGENCY=LENIENT REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE

10.	Samtools View/Index Filter out reads with a quality score below 20 from the SF5 alignment, convert it to a bam file, and create an index file.
http://www.htslib.org/doc/samtools-view.htm
a.	The bam file format is required for haplotype calling 
b.	samtools view -bS  <name_of_file.sam> -q 20 > <name_of_file.bam>
c.	Samtools index <name_of_file.bam>
d.	samtools view -bS SF5.rmdup.RG.sort.sam -q 20 > SF5.bam
e.	samtools index SF5.bam

Generating SNPs
Goal: generate a set of high-quality SNPs for M. nasutus SF5 pseudoreference genome construction.
11.	GATK HaplotypeCaller (in GVCF mode) Identify phased SNP and insertion/deletion variants from the SF5 alignment.
https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller
a.	Single sample GVCF calling without allele specific annotation was used because biallelic sites will be extracted and used, so multiallelic site calling doesn’t really matter
b.	gatk --java-options "-Xmx4g" HaplotypeCaller  \ -R <file/path/to/reference.fasta> \ -I  <file/path/to/name_of_file.bam> \ -O <file/path/to/output.g.vcf.gz> \ -ERC GVCF
c.	gatk --java-options "-Xmx4g" HaplotypeCaller -R ~/RNASeq/pseudoreference/ref/Mguttatus_V2refmtcp.fa -I SF5.bam -O SF5.g.vcf.gz -ERC GVCF
12.	GATK GenotypeGVCFs 
https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs
a.	gatk --java-options "-Xmx4g" GenotypeGVCFs \ -R <file/path/to/reference.fasta> \  -V <file/path/to/output.g.vcf.gz> \  -O <file/path/to/output.vcf.gz>
b.	gatk --java-options "-Xmx4g" GenotypeGVCFs -R ref/Mguttatus_V2refmtcp.fa -V SF5.g.vcf.gz -O SF5.vcf.gz

13.	GATK SelectVariants Extract only the biallelic SNPs from the vcf file containing the genotype calls for the SF5 alignment.
https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants
a.	gatk SelectVariants \ -R <file/path/to/reference.fasta> \ -V <file/path/to/output.vcf.gz> \ --select-type-to-include SNP \ --restrict-alleles-to BIALLELIC  -O <file/path/to/biallelic_snps.output.vcf>
b.	gatk SelectVariants -R ref/Mguttatus_V2refmtcp.fa -V SF5.vcf.gz --select-type-to-include SNP --restrict-alleles-to BIALLELIC  -O SF5.biallelic_snps.output.vcf

14.	GATK VariantFiltration Mark variants with “QUAL” in the filter column if the quality score (phred-scaled confidence that a variant is present) is less than 40 and “DEPTH” if the read depth in support of a variant site is less than 2.
https://gatk.broadinstitute.org/hc/en-us/articles/360037434691-VariantFiltration
https://commons.apache.org/proper/commons-jexl/reference/syntax.html
a.	gatk VariantFiltration \ -R <file/path/to/reference.fasta> \ -V <file/path/to/biallelic_snps.output.vcf> \ -O <file/path/to/biallelic_snps.filtered.output.vcf> \ --filterExpression "QUAL &lt; 40.0” --filterName “QUAL” \ --filterExpression “DP &lt; 2” --filterName “DEPTH” 
b.	gatk VariantFiltration -R ref/Mguttatus_V2refmtcp.fa -V SF5.biallelic_snps.output.vcf -O SF5.biallelic_snps.filtered.output.vcf --filter-name "QUAL" --filter-expression "MQ < 40.0" --filter-name "DEPTH" --filter-expression "DP < 2"

15.	GATK SelectVariants Remove all variants with anything other than ‘.’ or ‘PASS’ in the Filter field, in effect filtering out sites with mapping quality below 40 or read depth below two.
https://gatk.broadinstitute.org/hc/en-us/articles/360037055952-SelectVariants
a.	gatk SelectVariants \ -R <file/path/to/reference.fasta> \ -V <file/path/to/biallelic_snps.filtered.output.vcf> \ --select-type-to-include SNP \ --exclude-filtered=TRUE  -O <file/path/to/biallelic_snps.filtered.dropped.output.vcf>
b.	gatk SelectVariants -R ref/Mguttatus_V2refmtcp.fa -V SF5.biallelic_snps.filtered.output.vcf --select-type-to-include SNP --exclude-filtered=TRUE  -O SF5.biallelic_snps.filtered.dropped.output.vcf

16.	GATK FastaAlternateReferenceMaker Use filtered SNPs to generate a M. Nasutus SF5 pseudoreference genome.
https://gatk.broadinstitute.org/hc/en-us/articles/360037594571-FastaAlternateReferenceMaker#--snp-mask
a.	Don’t want to use masking file because it replaces SNPs with Ns
b.	Don’t need to provide intervals because we want SNPs across the whole genome
c.	gatk FastaAlternateReferenceMaker \ -R <file/path/to/reference.fasta> \ -O <reference_pseudo.fasta> \ -V <file/path/to/biallelic_snps.filtered.dropped.output.vcf> 
d.	~/gatk-4.1.7.0/gatk FastaAlternateReferenceMaker -R ref/Mguttatus_V2refmtcp.fa -O MguttatusV2refmtcp_SF5.fa -V SF5.biallelic_snps.filtered.dropped.output.vcf

17.	 Appending Allelic Identifiers: Add allelic identifiers (IM62 and SF5) to the chromosome names in the M. guttatus v2.0 reference and M. nasutus SF5 pseudoreference fasta files.
a.	constructed a diploid M. guttatus IM62-M. nasutus SF5 pseudoreference genome by appending allelic (i.e. IM62 or SF5) identifiers onto the chromosome names in the M. guttatus v2.0 reference and M. nasutus SF5 pseudoreference fasta files
b.	M. guttatus allelic reference:
i.	sed '/>/ s/$/_Allelic_ID/' <file/path/to/reference.fasta> > <file/path/to/reference_Allelic_ID.fasta>
ii.	sed '/>/ s/$/_IM62/' ~/mito_project/bams/ref/Mguttatus_V2refmtcp.fasta > Mguttatus_V2refmtcp_IM62.fa
c.	M. nasutus SF5 allelic pseudo reference:
i.	sed '/>/ s/$/_SF5/' <file/path/to/reference_pseudo.fasta> > <file/path/to/reference_pseudo_Allelic_ID.fasta>
ii.	sed 's/:.*/_SF5/' MguttatusV2refmtcp_SF5.fa | sed 's/.*\s/>/' >  MguttatusV2refmtcp_SF5_pseudo.fa

18.	Creating Diploid genome File: Manually merge the M.guttatus allelic reference and the M. nasutus SF5 allelic pseudo reference.
a.	cat <file/path/to/reference_Allelic_ID.fasta> <file/path/to/reference_pseudo_Allelic_ID.fasta> > <file/path/to/reference_Allelic_ID_reference_pseudo_Allelic_ID.fasta
b.	cat ref/Mguttatus_V2refmtcp_IM62.fa MguttatusV2refmtcp_SF5_pseudo.fa > MguttatusV2refmtcp_IM62_SF5_pseudo.fa

19.	Append Allelic Identifiers in GFF3 Files Create a diploid annotation file by appending allelic identifiers onto chromosome, gene, and transcript names, then manually combine the two files and convert the resulting GFF3 file to a GFT file. 
a.	The script to append allelic identifiers is called append_identifiers_GFF3.py and it’s file path on carnation is home/andrew_demaree/scripts/append_identifiers_GFF3.py. A copy has also been added to the RNA2020 Dropbox.The output of this scriptis two GFF3 files (SF5_Mguttatus_v2.0_256_gene_exons.gff3 and IM62_Mguttatus_v2.0_256_gene_exons.gff3)
b.	To change the input file, output file, allelic identifiers, or matched attributes, open the script using VIM (or any text editor that won’t alter the encoding) and change the appropriate string.
c.	After running the script, append the two output files together in the same order as the reference files.
i.	cat <file/path/to/Allelic_ID_1.gff3> <file/path/to/Allelic_ID_2.gff3> > <file/path/to/output/Allelic_ID_1_2.gff3>
ii.	cat IM62_Mguttatus_v2.0_256_gene_exons.gff3 SF5_Mguttatus_v2.0_256_gene_exons.gff3 > IM62_SF5_Mguttatus_v2.0_256_gene_exons.gff3
20.	Gffread Convert the merged gff3 file to a gft file using gffread.
https://github.com/gpertea/gffread
i.	 gffread my.gff3 -T -o my.gtf
ii.	gffread -T -o IM62_SF5_Mguttatus_v2.0_256_gene_exons.gtf

Part 2: RNAseq Alignment Map IM62 and SF5 genotypes to
the diploid pseudoreference genome.
1.	Download Single-end RNAseq Files: Use fastq-dump to download the three carpel and stamen replicates for each parental genotype (IM62 and SF5) listed in Table S1 from the supplemental data. 
a.	Fastq-dump (NCBI toolkit)
https://ncbi.github.io/sra-tools/fastq-dump.html
i.	Fastq-dump -O <path/to/output_dir/>  SRA_accesion_number_1...SRA_accesion_number_N

2.	Fastqc Check for Adapters, read quality, and quality scores before trimming.
a.	Examine the the html output files to look for read quality issues
b.	for i in *.fastq; do fastqc $i; done
c.	fastqc SRR*

3.	Trimmomatic Trim adapter sequences and low-quality bases from the raw RNAseq reads, then filter out reads shorter than 36 bp.
http://usadellab.org/cms/?page=trimmomatic
a.	java -jar /file/path/to/trimmomatic-0.35.jar SE \ -threads 8 -phred33 <file/path/to/forward.1.fastq> <file/path/to/reverse.2.fastq> \ <file/path/to/output_dir/trimmed_forward.1.paired.fastq> <file/path/to/output_dir/trimmed_forward.1.unpaired.fastq> \ <file/path/to/output_dir/trimmed_reverse.2.paired.fastq> <file/path/to/output_dir/trimmed_reverse.2.unpaired.fastq> \ ILLUMINACLIP:file/path/to/adapter/(illumina/Many.TruSeq.PE.fa):2:20:10:4 \ LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
b.	java -jar /home/thom_nelson/opt/Trimmomatic-0.35/trimmomatic-0.35.jar SE -threads 8 -trimlog SF5_stamen_1 SF5_stamen_1.fastq SF5_stamen_1.P.fastq ILLUMINACLIP:/home/thom_nelson/resources/Illumina/Many.TruSeq.PE.fa:2:20:10:4 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

4.	Fastqc Check for Adapters, read quality, and quality scores after trimming.
a.	for i in *.fastq; do fastqc $i; done
b.	fastqc *.P.fastq

5.	STAR Align the 36-50-bp Single-end RNAseq reads to the diploid pseudoreference generated in part 1.
https://github.com/alexdobin/STAR
a.	Generating Genome index files:
i.	--sjdbOverhang: get the read length from the trimmed read fastqc and subtract 1 from it.
ii.	genomeDir is the directory where the genome indices are kept. Need to make this directory before running and should be empty.
iii.	--sjdbGTFfile is the file path to the diploid GTF3 file generated for the pseudoreference.
iv.	~/STAR/STAR/bin/Linux_x86_64_static/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 … --sjdbGTFfile /path/to/annotations.gtf --sjdbOverhang ReadLength-1
v.	~/STAR/STAR/bin/Linux_x86_64_static/STAR --runThreadN 8 --runMode genomeGenerate -genomeDIr GenomeDir/ --genomeFastaFiles MguttatusV2refmtcp_IM62_SF5_pseudo.fa --sjdbGTFfile IM62_SF5_Mguttatus_v2.0_256_gene_exons.gtf --sjdbOverhang 50

b.	Mapping Reads:
i.	-- readFilesIn can take multiple samples in the form sample1,sample2,sample3,...,sampleN
ii.	Reads that overlap a SNP are expected to map uniquely to an SF5 or IM62 allele in the pseudoreference genome, while reads that do not overlap a SNP will map equally well to both alleles
iii.	limit each read to a maximum of two alignments by specifying --outFilterMultimapNmax 2
iv.	Randomly designate one of the two allowed alignments as primary by specifying --outMultimapperOrder Random
v.	~/STAR/STAR/bin/Linux_x86_64_static/STAR --runThreadN NumberOfThreads --genomeDir /path/to/genomeDir --readFilesIn /path/to/read1 /path/to/read2 …  /path/to/readN --genomeDir <file/path/to/genomeDir> --outFilterMultimapNmax 2 --outMultimapperOrder Random
vi.	for i in *.P.fastq; do ~/STAR/STAR/bin/Linux_x86_64_static/STAR --runThreadN 8 --genomeDir GenomeDir/ --readFilesIn $i --outFilterMultimapNmax 2 --outMultimapperOrder Random --outFileNamePrefix ${i:0:-8}; done

6.	SAMtools view Remove secondary alignments and unmapped reads, leaving just unique (i.e. allele-specific) and primary alignments.
a.	-F tells samtools to omit sequences with the int value that follows the tag (e.g. 4 = reads that are unmapped, 100 = reads that are not primary)
b.	samtools view -b -F 4 -F 100 <file/path/to/aligned.bam> > <file/path/to/aligned.filtered.bam>
c.	for i in *.sam; do samtools view -b -F 4 -F 100 $i > ${i:0:-15}.filtered.bam; done

7.	SAMtools sort Convert the filtered bams to coordinates sorted bams.
a.	samtools sort <file/path/to/aligned.filtered.bam> >  <file/path/to/aligned.filtered.sort.bam>
b.	for i in *filtered.bam; do samtools sort $i > ${i:0:-4}.sort.bam; done

8.	Picard Mark Duplicates Remove optical and  PCR duplicates from the RNAseq alignments. 
https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-
a.	Requires coordinate sorted bam or sam as input 
b.	java -jar <file/path/to/java/jar/picard.jar MarkDuplicates \ INPUT=<file/path/name_of_coordinate_sorted.bam> OUTPUT=<file/path/name_of_coordinate_sorted.rmdup.bam> \ METRICS_FILE=<name_of_file.rmdup_metrics_fix \ VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE
c.	for i in *.sort.bam; do java -jar /home/thom_nelson/opt/picard.jar MarkDuplicates INPUT=$i OUTPUT=${i:0:-4}.rmdup.filtered.bam METRICS_FILE=SF5.rmdup_metrics_fix VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE; done

9.	Picard Create Sequence Dictionary Creates a .dict file required by future tools.
https://gatk.broadinstitute.org/hc/en-us/articles/360036729911-CreateSequenceDictionary-Picard-
a.	java -jar picard.jar CreateSequenceDictionary \ R=<file/path/name_of_pseudoreference.fa
b.	java -jar/home/thom_nelson/opt/picard.jar CreateSequenceDictionary R=MguttatusV2refmtcp_IM62_SF5_pseudo.fa

10.	Samtools faidx Creates a .fai index file required by future tools. 
http://www.htslib.org/doc/samtools-faidx.html
a.	samtools faidx <file/path/to/name_of_pseudoreference.fa>
b.	samtools faidx MguttatusV2refmtcp_IM62_SF5_pseudo.fa

11.	GATK SplitNCigarReads Parse intron-spanning reads into exon segments and trim off bases that extend into intronic regions.
a.	Splits reads that contain Ns in their Cigar String into K+1 sequences, where k is the number of Ns.
b.	gatk SplitNCigarReads \ -R <file/path/to/Mguttatus_V2refmtcp_Mnasutus_pseudo.fasta \ -I <file/path/name_of_coordinate_sorted.rmdup.bam> \ -O <file/path/name_of_coordinate_sorted.rmdup.split.bam>
c.	for i in *rmdup.filtered.bam; do ~/gatk-4.1.7.0/gatk SplitNCigarReads -R pseudoreference/MguttatusV2refmtcp_IM62_SF5_pseudo.fa -I $i -O ${i:0:-4}.split.bam; done

12.	SAMtools view filter alignments based on mapping quality.
a.	Obtain allele-specific and primary mapping (i.e. total) reads using a Q20 threshold.
b.	samtools view -b <file/path/name_of_coordinate_sorted.rmdup.split.bam> -q 20 > <file/path/name_of_coordinate_sorted.rmdup.split.Q20.bam>
c.	for i in *.split.bam; do samtools view -b -q 20 -o ${i:0:-4}.q20.bam $i;done 

13.	Htseq-count 
https://htseq.readthedocs.io/en/release_0.11.1/index.html
a.	I have created a python virtual environment called htseq_new which has Htseq-count installed otherwise follow the directions in the htseq documentation.
b.	source activate htseq_new
c.	htseq-count -f bam --nonunique=all --stranded=reverse <file/path/to/sample.split.q20.bam <file/path/to/reference_genes_exons.gtf>
d.	source deactivate
e.	source activate htseq_new
f.	for i in *.split.q20.bam; do htseq-count -f bam --nonunique=all --stranded=reverse $i ~/RNASeq/trsht/IM62_SF5_Mguttatus_v2.0_256_gene_exons.gtf > ${i:0:-40}.counts.tsv;done
g.	Source deactivate


Useful Commands:
Count read totals: 
1.	cat SF5_carpel_1.P.fastq | echo $((`wc -l`/4))
Count primary mapped reads: 
2.	samtools view -F 4 -F 256 -F 1024 SF5_carpel_1.filtered.bam | wc -l

