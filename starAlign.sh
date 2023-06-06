names="AHQF1A_YB_19
AHQF1A_YB_20"

dataDir="/home/colette_berg/YNP/RNA_2021/IM/fastqs"
genomeDir="/home/colette_berg/YNP/NT_fastqs/genome"
outDir="/home/colette_berg/YNP/NT_fastqs/bams"

for name in ${names}

do 

gunzip ~/YNP/RNA_2021/IM/fastqs/${name}.1.paired.fastq.gz
gunzip ~/YNP/RNA_2021/IM/fastqs/${name}.2.paired.fastq.gz

~/resources/packages/STAR-2.7.9a/bin/Linux_x86_64_static/STAR \
    --runThreadN 8 \
    --genomeDir ${genomeDir}  \
    --readFilesIn ${dataDir}/${name}.1.paired.fastq ${dataDir}/${name}.2.paired.fastq \
    --sjdbGTFfeatureExon exon \
    --sjdbGTFtagExonParentTranscript ID \
    --sjdbGTFtagExonParentGene Parent \
    --sjdbGTFfile ${genomeDir}/Tv1_N_pseudo.gff3 \
    --outFilterMultimapNmax 2 \
    --outMultimapperOrder Random \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${outDir}/${name}_Tv1

gzip ~/YNP/RNA_2021/IM/fastqs/${name}.1.paired.fastq
gzip ~/YNP/RNA_2021/IM/fastqs/${name}.2.paired.fastq

done
