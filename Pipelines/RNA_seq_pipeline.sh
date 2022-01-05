#!/bin/bash
#author : @Julien Nguinkal
#set -euxo  
 ## Wrapper script for maping paired end data with STAR aligner ##
 ## This script accepts the absolute folder containing the paired reads, ###
 ## the genome file, a prefixfolder for the index  and Outfile Prefix  ####
 ##########################################################################
 
 
# if [ $# -lt 4 ]; then
 # echo "Usage: $0 PathInputReads genonmeFile(in fasta) PrefixGenomeIndex OutfilePrefix "
 # exit 1
 #fi
 
#date 
printf "######## STAR RNA-Seq Mapping pipeline ###### \n"
#VERS=`STAR --version`

 
 
 ### setting important variables ##
 
  #ReadsDir="$1"  # Path to reads files
  #GENOME="$2" # The genome file in working directory
  #Ge#nIndexDir="$3" # Index folder
  #OutfilePrefix="$4"  #Prefix for output index
  
  STAR="/disk2/nguinkal/ZanderProject/pipelines/STAR-MAPPING/STAR-2.7.5a/bin/Linux_x86_64_static/STAR"  # STAR executable
  HTSEQCOUNT="/home/nguinkal/anaconda3/envs/ragooEnv/bin/htseq-count" ## htseq-count executable
  ANNOT="GCF_008315115.2_SLUC_FBN_1.2_genomic.gtf"
  CPUS=128 # number of threads 
  READLENGTH=149 ## read leng -1 
  GENOME="SLUC_refGenome.fa"
  ### Create the genome index ###
  
  mkdir -p GenomeDir # create dir if not exists
  cp $GENOME GenomeDir #copy the genome inside that folder
 
  
  printf "###### GENERATING GENOME in GenIndexDir... ###### \n"
# $STAR --runThreadN 32 \
 #  --runMode genomeGenerate \
  #--genomeDir GenomeDir \
  #--genomeFastaFiles $GENOME \
  #--sjdbOverhang $READLENGTH  \
  #--sjdbGTFfile $ANNOT \
  #--genomeSAindexNbases 13;
   
  echo  "#### GENOME INDEXING WAS SUCESSFUL! ##### "
  
  echo  "#### STARTING MAPPING READS .. ##### "
  
  #mkdir -p aligned ## save aligned reads in bam files here
  $STAR --genomeLoad LoadAndExit --genomeDir GenIndexDir
  



for a in $(seq 1 2 $(wc -l <(ls SL*.fq.?) | cut -f1 -d " "));

   do

      i=$(cut -f $a <(sort <(ls  SL*.fq.?) | tr -s \\n \\t ));
      j=$(cut -f $(( $a + 1 )) <(sort <(ls SL*.fq.?) | tr -s \\n \\t ));


  printf "###### MAPPING SAMPLE ${i%%.blacklist_paired_unaligned.fq.1}  ###### \n"
 
 # mkdir -p ${i%%.blacklist_paired_unaligned.fq.1} ## alignment folder
 $STAR --genomeDir GenomeDir \
--readFilesIn $i,$j \
--runThreadN 32 \
--outFileNamePrefix ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1} \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFilterType BySJout \
--sjdbGTFfile $ANNOT \
--twopassMode Basic ; 

### sort with  samtools ## 
 samtools sort -@80 ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}*.out.bam ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}.sorted


  printf "###### READS QUANTIFICATION  WITH HTSEQ-COUNT FOR SAMPLE ${i%%.blacklist_paired_unaligned.fq.1} ... ###### \n"


### counting with htseq-count 
$HTSEQCOUNT -f bam -s yes -m union  -i gene_id \
           ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}*.bam $ANNOT \
			&> ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}.htseq.log
 grep gene ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}.htseq.log | sed 's/gene://g' > ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}.counts.csv 

### quantification with feature count
#featureCounts -T 32  -G  -p  -g gene_id -o ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}.featuresCount.txt ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}.sorted.bam -a $ANNOT
#cut -f1,7,8,9,10,11,12 ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}.featuresCount.txt > ${i%%.blacklist_paired_unaligned.fq.1}/${i%%.blacklist_paired_unaligned.fq.1}.featureCount.Rmatrix.txt
	  




echo "HTSEQ COUNTING FOR ${i%%.blacklist_paired_unaligned.fq.1}  SUCESSFUL! "
done 

	 printf "######  MAPPING AND SUMMARIZATION  FOR ALL SAMPLES FINISHED... GENERATING TABLE...###### \n"

mkdir -p COUNTS
cp SL*/SL*.counts.csv COUNTS && cd COUNTS

for i in *csv; do sed -i "1igene\t${i%%.counts.csv}" ${i%%.counts.csv} ; done # add column names
N=$(($(ls -l *.csv | wc -l)*2)) # count number of files
paste *csv | cut -f 1,$(seq -s, 2 2 $N) > all_HTSeq.csv # merge and keep only one column with gene names
cd ..

printf "######  ALL DONE! BYE. ###### \n"

