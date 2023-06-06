genomeDir="/home/colette_berg/YNP/NT_fastqs/genome"

~/resources/packages/STAR-2.7.9a/bin/Linux_x86_64_static/STAR --genomeSAindexNbases 13 --runThreadN 8 --runMode genomeGenerate --genomeDir ${genomeDir} --genomeFastaFiles ${genomeDir}/Tv1_N_pseudo.fasta --sjdbGTFfile ${genomeDir}/Tv1_N_pseudo.gff3