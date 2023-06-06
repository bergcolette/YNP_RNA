
#for i in *.bam; do samtools view -b -F 4 -F 100 $i > ${i}.filtered.bam; done

#for i in *filtered.bam; do samtools sort $i > ${i}.sort.bam; done

for i in *.sort.bam; do java -jar /home/colette_berg/resources/packages/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar MarkDuplicates INPUT=$i OUTPUT=${i:0:-4}.rmdup.filtered.bam METRICS_FILE=SF5.rmdup_metrics_fix VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE REMOVE_SEQUENCING_DUPLICATES=TRUE; done

java -jar /home/colette_berg/resources/packages/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar CreateSequenceDictionary -R ~/YNP/NT_fastqs/genome/Tv1_N_pseudo.fasta

for i in *rmdup.filtered.bam; do /home/colette_berg/resources/packages/gatk-4.2.6.1/gatk-package-4.2.6.1-local.jar SplitNCigarReads -R ~/YNP/NT_fastqs/genome/Tv1_N_pseudo.fasta -I $i -O ${i:0:-4}.split.bam; done

for i in *.split.bam; do samtools view -b -q 20 -o ${i:0:-4}.q20.bam $i;done
