#!/bin/bash


#Nombre de gènes de résistance à certains antibiotiques dans le génome pour l'échantillon G est de 31 environ



#2 arguments positionnels: wgs.sh dossier_reads_bruts dossier_sortie
dossier_reads_bruts=$1
dossier_sortie=$2

mkdir -p $dossier_sortie/bwt2

./soft/bowtie2-build databases/all_genome.fasta all_genome
./soft/bowtie2 -x all_genome -1 $dossier_reads_bruts/EchG_R1.fastq.gz -2 $dossier_reads_bruts/EchG_R2.fastq.gz --end-to-end --fast -S $dossier_sortie/bwt2/echG.sam

./soft/samtools-1.6/samtools view -S -b $dossier_sortie/bwt2/echG.sam > $dossier_sortie/bwt2/echG.bam
./soft/samtools-1.6/samtools sort $dossier_sortie/bwt2/echG.bam > $dossier_sortie/bwt2/echG.sorted.bam
./soft/samtools-1.6/samtools index $dossier_sortie/bwt2/echG.sorted.bam
./soft/samtools-1.6/samtools idxstats $dossier_sortie/bwt2/echG.sorted.bam
grep ">" databases/all_genome.fasta | cut -f 2 -d ">" > association.tsv

./soft/megahit -1 $dossier_reads_bruts/EchG_R1.fastq.gz -2 $dossier_reads_bruts/EchG_R2.fastq.gz --k-list 21 --mem-flag 0 -o ./megahit_out

./soft/prodigal -i megahit_out/final.contigs.fa -d ./prodigal.out.fna

sed "s:>:*\n>:g" prodigal.out.fna | sed -n "/partial=00/,/*/p" | grep -v "*" > genes_full.fna

./soft/blastn -query genes_full.fna -db databases/resfinder.fna -outfmt '6 qseqid sseqid pident qcovs evalue' -out blast.out -evalue 0.001 -qcov_hsp_perc 80 -perc_identity 80 -best_hit_score_edge 0.001
