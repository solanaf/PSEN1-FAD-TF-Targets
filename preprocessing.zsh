
################################# PREPROCESSING PIPELINE! #######################################
# Project: PSEN1 FAD TFs
# Creator: Solana Fernandez
##################################################################################################
##################################################################################################
################    NOTE: export $GET_DATA and $DIRECTORYbeforehand 	  ########################
################	             				          ########################
################                                                          ########################
################  EXAMPLE USAGE (FEEL FREE TO UNCOMMENT AND RUN with):    ########################
################       (no spaces between '=' when defining)    	  ########################
# export DIRECTORY="/Users/solanafernandez/Desktop/UCSD_local/BENG204/project" # desired location of file dump
# export GET_DATA="test1.sh" # 	location of wget bash script

##################################################################################################
######## --------------> $DIRECTORY should be the path that contains zsh data retrieval folder
cd $DIRECTORY

mkdir raw_reads

##########################  Retrieve Data: Test (make sure executable) ############################
chmod +x zsh/$GET_DATA
cd raw_reads
zsh zsh/$GET_DATA
cd ..

##########################  FastQC (input: fastq; output: fastq) ##################################  
mkdir fastqc
find /raw_reads -name "*.fastq.gz" | xargs fastqc -t 8 -o fastqc

################################## TrimGalore #####################################################
######## --------------> prep fastq files to iterate through IN ORDER

find . -type f -name "*_1.fastq.gz" > 1_end.txt
find . -type f -name "*_2.fastq.gz" > 2_end.txt
echo "1 end to trim: " && cat 1_end.txt
sort 1_end.txt -o 1_end.txt
echo "2 end to trim: " && cat 2_end.txt
sort 2_end.txt -o 2_end.txt

mkdir trimmed_reads
paste -d ' ' 1_end.txt 2_end.txt | xargs -n 2 trim_galore --paired --fastqc -o trimmed_reads

######## --------------> Delete raw reads for space
rm -rf raw_reads
###################################  ALIGNMENT #####################################################
######## --------------> Retrieve hg38
if [ ! -f "hg38.fa" ]; then
    echo "Indexed hg38 file does not exist. Retrieving..."
    mkdir hg38_index
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz 
    gunzip hg38.fa.gz 
    rm -rf hg38.fa.gz 
    bowtie2-build hg38.fa hg38_index/hg38
    rm -rf hg38.fa
fi

######## --------------> Retrieve pair-end reads
find . -type f -name "*_1.fq.gz" > 1_end.txt
find . -type f -name "*_2.fq.gz" > 2_end.txt
echo "1 end to align: " && cat 1_end.txt
sort 1_end.txt -o 1_end.txt
echo "2 end to align: " && cat 2_end.txt
sort 2_end.txt -o 2_end.txt

######## --------------> For multiple reads, will output one SAM/BAM file
bowtie2 -x hg38_index/hg38 \
    -1 $(tr '\n' ',' < 1_end.txt | sed 's/,$//') \
    -2 $(tr '\n' ',' < 2_end.txt | sed 's/,$//') \
    -S aligned.sam --threads 8 | samtools view -@8 -bS | samtools sort -@8 -o aligned_reads/aligned.sorted.bam


######## --------------> Delete trimmed reads & SAM files for space
rm -rf aligned_reads/*.sam trimmed_reads 1_end.txt 2_end.txt

######## --------------> Get GTF/GFF
if [ ! -f "hg38.refGene.gtf" ]; then
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz
	gunzip hg38.refGene.gtf.gz
	rm -rf hg38.refGene.gtf.gz
fi

######## --------------> Count reads
featureCounts -T 8 -p -t exon -g gene_id -a hg38.refGene.gtf -o gene_counts.txt aligned_reads/*.bam 
 
###################################  Move to New Folder #############################################
GET_DATA=${GET_DATA//./_}
mkdir $GET_DATA
sudo mv aligned_reads gene_counts.txt gene_counts.txt.summary fastqc $GET_DATA/

########################### COMPLETION NOTIFICATION #################################################
Echo "Process Complete."
afplay /System/Library/Sounds/Purr.aiff
