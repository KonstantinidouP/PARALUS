#!/bin/bash
# PARAlus is a bedtools and awk based program to analyze PAR-CLIP or CLIP data
# created by Dimitrios Anastasakis
# PARalus 1.5.1 this version was with addition of --outSAMmultNmax 1 --outMultimapperOrder Random 
# PARalus 1.5.3 this version was with addition of rRNA mapper toget rRNA maping stats and also with a cluster.csv cleaner with priority in annotations as 3ss > 5ss > start_codon > stop_codon> 3utr> 5utr > CDS > exon > intron 

unset ALL_PROXY

prog=$0
scriptdir=$(dirname "${prog}")

if [ $# -eq 0 ]; then
  echo 'PARalus v1.4.3
Please give me something to do 
No input files provided
Please type "--help" for instructions'
  exit 0
fi

if [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
  cat << 'EOF'
  
          ||
       ##################
       ##################   
       ###Paralus1.5.1###
       ##################
       ##################
      ##################
           ||         
   &  __@_/  ||       &
   \ \  ||    ||       //
    \__\_/\___@___@___@___@____//   
      \   /  /  /  /  /
  ~~~~~~~~\____/___/___/___/____))))~~~~~~~~~~
   ~ ~~   /  /  /  /   ~~   ~
   ~~    /  /  /  /  ~    ~~~
 
 Usage: ./Paralus -x genome_index_base -i input_file -G genome.gtf 
 
 Required fields:
  
  -i  Your input file (fastqz/fastq/fasta/bam/SRA) 
  -x  Genome index folder
    		May build one like this:
  			make dir /path/Genome_index
  			STAR --runMode genomeGenerate --genomeDir /path/Genome_index --genomeFastaFiles rRNA_genome.fasta
  -A  Annotation. A 6 or 7 fieled bed file containing 
 			Can make this with the masterchef script using a gtf (preferable gencode to cotnain gene_type infomration
 			You may add a bed file for your favourite organism from tRNA Data Base
 			bash /path/maserchef.sh -g genes.gtf -t tRNA_DB.bed
 			

  *Make sure you have loaded samtools bedtools and cutadapt acording to your enviroment
  
  
  
  
 Optional parameters:
  
  Preprocess:
  
  -n  Minimum length allowed after trimming. Default: 19
  -m  Minimum number of aligned reads to consider. Default: 3
  -a  Your 3' adapter. Default: No adapter
  -d  Your 5' adapter. Default: GTTCAGAGTTCTACAGTCCGACGATC
  -u  Number of wildcards anchored to the 5' end of the sequence. Default: 0
     (do not add this if 5' adapter is barcoded and ends with UMIs)
  -p  Threads. Default: 8 dont use more than 12!! else star wont run
  -e  rRNA genome to pre map sequnces and provide % of rRNA reads 
  		May build one like this:
  		make dir /path/rRNA_genome
  		STAR --runMode genomeGenerate --genomeDir /path/rRNA_genome --genomeFastaFiles rRNA_genome.fasta
  
  PAR-CLIP binding site identification - parameters:
  
  -c  Conversion location count. How many TtoC location to be considered cluster.(1-2) Default: 1
  -r  A group will be considered a cluster if TtoC containing reads are equal or more than this fraction: 0.3 (30%) 
  -f  If genome.fa (<Genome index basename>.fa) does not exist in the subfolder of genome index,
     please provide a the genome.fa using the -f option or copy the genome.fa into the genome index folder
  
#######################
EOF
  exit 0
fi

minlen=19
minRead=3
five_adapt=GTTCAGAGTTCTACAGTCCGACGATC
umis=0
f=fq
threads=8
convloc=1
fraction=0.3

while getopts hi:x:A:p:r:v:n:m:a:d:u:s:e:c:f:L: OPT; do
  case "$OPT" in
    i)
      input=$OPTARG
      ;;
    x)
		index=$OPTARG
      ;;
    A)
      annotation=$OPTARG
      ;;
    p)
      threads=$OPTARG
      ;;
    r)
      fraction=$OPTARG
      ;;
    n)
      minlen=$OPTARG
      ;;
    m)
      minRead=$OPTARG
      ;;
    a)
      three_adapt=$OPTARG
      ;;
    d)
      five_adapt=$OPTARG
      ;;
    u)
      umis=$OPTARG
      ;;

    c)
    	convloc=$OPTARG
      ;;
    f)
    	genomefasta=$OPTARG
      ;;
     e)
    	rRNA_genome=$OPTARG
      ;;

    h) usage
    exit 0;; 
  esac
  
done
input=$(basename "${input}")
name=${input%.fast*}

file_dir=$(dirname "${input}")


ext="${input##*.}"


genomedir=$(dirname "${index}")
rRNA_genomedir=$(dirname "${rRNA_genome}")


cd $file_dir

if [[ ! -f "${index}/SAindex" && ! -f "${genomefasta}" ]]; then
echo "error:
please provide a genome.fa using the -f option or copy the genome.fa into the genome index folder"
 exit 0

elif [ ! -f "${annotation}" ]; then
echo "error:
please provide a annotation file"
 exit 0
elif [ ! -f "${input}" ]; then
echo "error:
please provide an input file (fastqz/fastq/fasta/"
 exit 0
 elif [[ ! -f "${rRNA_genome}/SAindex" ]]; then
echo "Warning:
rRNA_genome was not provided"
else
echo "Input files provided"
fi
if [ -f "${name}"_stats.txt ] ; then rm "${name}"_stats.txt ; fi 
if [ -f "${name}".log ] ; then rm "${name}".log ; fi 
echo "Your parameters:" >> "${name}"_stats.txt
echo "Input file : "${input}"" >> "${name}"_stats.txt
echo "genome file : "${index}"" >> "${name}"_stats.txt
echo "gtf file : "${gtf}"" >> "${name}"_stats.txt
echo "threads : "${threads}"" >> "${name}"_stats.txt
echo "Minimun TtoC fraction : "${fraction}"" >> "${name}"_stats.txt
echo "Minimun conversion locations : "${convloc}"" >> "${name}"_stats.txt
echo "Minimun reads for group definition : "${minRead}"" >> "${name}"_stats.txt
echo "Minumun length allowed by cutadapt : "${minlen}"" >> "${name}"_stats.txt
echo "Multimaping allowed (m option in bowtie) : "${tmapa}"" >> "${name}"_stats.txt
echo "Your 3' adapter : "${three_adapt}"" >> "${name}"_stats.txt
echo "Your 5' adapter : "${five_adapt}"" >> "${name}"_stats.txt
echo "bases trimed from the 5' end of the read ( wildcards after 5' adapter) : "${umis}"" >> "${name}"_stats.txt
echo "genome for Homer : "${genome}"" >> "${name}"_stats.txt
num_of_lines=$(wc -l < "${input}")
echo "number of lines in fastq: ${num_of_lines}" >> "${name}"_stats.txt
num_of_seq=$(echo $(( num_of_lines / 4 )))
echo "number of sequences : ${num_of_seq}" >> "${name}"_stats.txt

##########################################################################################################################
## preprocessing #in theory adapters should be removed prior 
## collapsing but error rates of illumina are < 0.1% making it insignificant #
##########################################################################################################################


if [[ ! -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -gt 0 && $umis -gt 0 ]]
	then
		echo "You have a fastq.gz file and have UMIs and adapters and did not provide a rRNA genome"
		zcat "${name}".fastq.gz | awk 'NR%4==2' | sort --parallel="$threads" \
		| uniq -c | sed 's/\s\s//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -u "$umis" -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		 -o "${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt	

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn "${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 

elif [[ ! -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -gt 0 && $umis -eq 0 ]]
	then
		echo "You have a fastq.gz file and dont provide UMIs and rRNA genome  and provide adapters a"
		zcat "${name}".fastq.gz | awk 'NR%4==2' | sort --parallel="$threads" \
		| uniq -c | sed 's/\s\s//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		 -o "${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt	

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn "${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 


elif [[ ! -f ${rRNA_genome}/SAindex && $ext == "fastq" && ${#three_adapt} -gt 0 && $umis -gt 0 ]]
	then
		echo "You have a fastq file and have UMIs and adapters and did not provide a rRNA genome"
		awk 'NR%4==2' "${name}".fastq | sort --parallel="$threads" | uniq -c | sed 's/\s\s//g' \
		| awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -u "$umis" -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		-o "${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt	

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn "${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 

elif [[ ! -f ${rRNA_genome}/SAindex && $ext == "fasta" && ${#three_adapt} -eq 0 && $umis -eq 0 ]]
	then
		echo "You have a fasta file and did not provide UMIs and adapters and a rRNA genome"
  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn "${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 

elif [[ -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -gt 0 && $umis -gt 0 ]]
	then
				echo "You have a fastq.gz file and have UMIs and adapters and provided a rRNA genome"
		zcat "${name}".fastq.gz | awk 'NR%4==2' | sort --parallel="$threads" | uniq -c \
		| sed 's/\s\s//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -u "$umis" -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		-o pre_rRNA_dep_"${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt

	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" --readFilesIn pre_rRNA_dep_"${name}"_trimmed.fasta \
	--outFileNamePrefix rRNA_aligned_"${name}" --outSAMtype BAM SortedByCoordinate --outSAMattributes NH MD HI AS \
	--runThreadN "$threads" --outReadsUnmapped Fastx
mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta 

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log
elif [[ -f ${rRNA_genome}/SAindex && $ext == "fastq" && ${#three_adapt} -gt 0 && $umis -gt 0 ]]
	then
			echo "You have a fastq file and have UMIs and adapters and provided a rRNA genome"
			
		awk 'NR%4==2' "${name}".fastq | sort --parallel="$threads" | uniq -c | sed 's/\s\s//g' \
		| awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -u "$umis" -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		-o pre_rRNA_dep_"${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt	

	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" \
	--readFilesIn pre_rRNA_dep_"${name}"_trimmed.fasta --outFileNamePrefix rRNA_aligned_"${name}" \
	--outSAMtype BAM SortedByCoordinate --outSAMattributes NH MD HI AS --runThreadN "$threads" --outReadsUnmapped Fastx

mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." \
	--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log
elif [[ -f ${rRNA_genome}/SAindex && $ext == "fasta" && ${#three_adapt} -eq 0 && $umis -eq 0 ]]
	then
		echo "You have a fasta file and no UMIs and adapters and provided a rRNA genome"
		
	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" --readFilesIn "${input}" \
	--outFileNamePrefix rRNA_aligned_"${name}" --outSAMtype BAM SortedByCoordinate \
	--outSAMattributes NH MD HI AS --runThreadN "$threads" --outReadsUnmapped Fastx
mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta 


	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." \
	--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log
elif [[ ! -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -eq 0 && $umis -eq 0 ]]
	then
		echo "You have a fastq.gz file and not provie UMIs, adapters and a rRNA genome"
	
		zcat "${name}".fastq.gz | awk 'NR%4==2' | sort --parallel="$threads" \
		| uniq -c | sed 's/\s\s//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		echo "aligning sequences to genome" >> "${name}".log
  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn ${name}_collapsed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log

elif [[ -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -eq 0 && $umis -eq 0 ]]
	then
	echo "You have a fastq.gz file and no UMIs and adapters and provided a rRNA genome"
		zcat "${name}".fastq.gz | awk 'NR%4==2' | sort --parallel="$threads" \
		| uniq -c | sed 's/\s\s//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		echo "aligning sequences to rRNA_genome" >> "${name}".log
	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" \
	--readFilesIn ${name}_collapsed.fasta --outFileNamePrefix rRNA_aligned_"${name}" \
	--outSAMtype BAM SortedByCoordinate --outSAMattributes NH MD HI AS --runThreadN "$threads" --outReadsUnmapped Fastx
mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta
		echo "aligning sequences to genome" >> "${name}".log
  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." \
	--outSAMtype BAM SortedByCoordinate --outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 

elif [[ -f ${rRNA_genome}/SAindex && $ext == "gz" && ${#three_adapt} -gt 0 && $umis -eq 0 ]]
	then
				echo "You have a fastq.gz file and have UMIs and adapters and provided a rRNA genome"
		zcat "${name}".fastq.gz | awk 'NR%4==2' | sort --parallel="$threads" | uniq -c \
		| sed 's/\s\s//g' | awk '{OFS= "\n"; print ">"NR"-"$1,$2}' > ${name}_collapsed.fasta
		cutadapt -m "$minlen" -a "$three_adapt" -g "$five_adapt;min_overlap=10" \
		-o pre_rRNA_dep_"${name}"_trimmed.fasta ${name}_collapsed.fasta >> "${name}"_stats.txt

	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --genomeDir "${rRNA_genome}" --readFilesIn pre_rRNA_dep_"${name}"_trimmed.fasta \
	--outFileNamePrefix rRNA_aligned_"${name}" --outSAMtype BAM SortedByCoordinate --outSAMattributes NH MD HI AS \
	--runThreadN "$threads" --outReadsUnmapped Fastx
mv rRNA_aligned_"${name}"Unmapped.out.mate1 rRNA_depleted_"${name}"_trimmed.fasta 

  	STAR --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterMultimapScoreRange 0 \
	--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignEndsType EndToEnd \
	--alignSJDBoverhangMin 999 --alignIntronMax 1 --runThreadN "$threads" --genomeDir "${index}" \
	--readFilesIn rRNA_depleted_"${name}"_trimmed.fasta --outFileNamePrefix "${name}." --outSAMtype BAM SortedByCoordinate \
	--outSAMmultNmax 1 --outMultimapperOrder Random --outSAMattributes NH MD HI AS 
		echo "done!" >> "${name}".log

echo "Could not find any combinations of options for the input"


fi

##############
# index bam ##
##############

samtools index "${name}".Aligned.sortedByCoord.out.bam


###############################
## converting bamtobed    ##
## good for bedtools coverage #
###############################


bedtools bamtobed -i "${name}".Aligned.sortedByCoord.out.bam > "${name}".aligned_sorted.bed

###############################
## get TtoC bam and index it ##
###############################

samtools calmd -e "${name}".Aligned.sortedByCoord.out.bam "${genomefasta}" | awk -F "\t" '(($1 ~"^@") || ($2 == "0" || $2 == "256") && ($10 ~ "C" && $13 ~ "T")) || (($1 ~"^@") || ($2 == "16" || $2 == "272") && ($10 ~ "G" && $13 ~ "A"))' | samtools view -Sbh --threads "$threads" - | samtools sort --threads "$threads" - > "${name}".aligned_TtoC.bam


samtools index "${name}".aligned_TtoC.bam

########################################################################################
### make TtoC bed to run coverage faster (bedtools coverage with bams would collapse ###
########################################################################################


bedtools bamtobed -i "${name}".aligned_TtoC.bam > "${name}".aligned_TtoC.bed


##############################
#### extract TtoC locationg ##
##############################


samtools view --threads "$threads" "${name}".aligned_TtoC.bam | awk '{OFS="\t" } { print $1,$2,$3,$4,$10,$13}' | sed -e 's/MD:Z://g' -e 's/\^//g' -e 's/=*//g' | awk 'BEGIN{FS=OFS="\t"} {gsub("A|G|C|T", "\t&\t", $6)} {gsub("A|G|C|T", "&\t", $5)}1' | awk 'NF >5' | awk '{OFS="\t"} \
  { if ((NF == 8) && ($2 == 0 || $2 == 256) && ($5~"C" && $7~"T")) \
    { print $3,$4+$6-1,$4+$6,$1,".","+"} \
  else if ((NF == 8) && ($2 == 16 || $2 == 272) && ($5~"G" && $7~"A")) \
    { print $3,$4+$6-1,$4+$6,$1,".","-"}
  else if ((NF == 11) && ($2 == 0 || $2 == 256) && ($5~"C" && $8~"T") && ($6~"C" && $10~"T")) \
    { print $3,$4+$7-1,$4+$7,$1,".","+","\n"$3,$4+$7+$9,$4+$7+$9+1,$1,".","+"}
  else if ((NF == 11) && ($2 == 0 || $2 == 256) && ($5~"C" && $8~"T") && ($6!~"C" || $10!~"T")) \
    { print $3,$4+$7-1,$4+$7,$1,".","+"}
  else if ((NF == 11) && ($2 == 0 || $2 == 256) && ($5!~"C" || $8!~"T") && ($6~"C" || $10~"T")) \
    { print $3,$4+$7+$9,$4+$7+$9+1,$1,".","+"}
  else if ((NF == 11) && ($2 == 16 || $2 == 272) && ($5~"G" && $8~"A") && ($6~"G" && $10~"A")) \
    { print $3,$4+$7-1,$4+$7,$1,".","-","\n"$3,$4+$7+$9,$4+$7+$9+1,$1,".","-"}
  else if ((NF == 11) && ($2 == 0 || $2 == 256) && ($5~"G" && $8~"A") && ($6!~"G" || $10!~"A")) \
    { print $3,$4+$7-1,$4+$7,$1,".","-"}
  else if ((NF == 11) && ($2 == 0 || $2 == 256) && ($5!~"G" || $8!~"A") && ($6~"G" || $10~"A")) \
    { print $3,$4+$7+$9,$4+$7+$9+1,$1,".","-"}}'  | bedtools sort | bedtools merge -d -1 -s -i - -c 6 -o distinct \
| awk '{print $1,$2,$3,"TtoC",".",$4}' | sed 's/\s/\t/g' | bedtools sort > "${name}".aligned_TtoC_loc.bed


######################
## merging bedfile ##
#####################


## the ($2 != $4 || $3 != $5 ) means that group of reads that have exact same start and end are discarded 
## these are artifacts
#########################################################################################################

bedtools merge -s -i "${name}".aligned_sorted.bed -c 2,3,4,5,6,2 -o mean,mean,count,distinct,distinct,collapse \
	| awk -v var=${minRead} '($6 >= var ) && ($3-$2 > 1) && ($2 != $4 || $3 != $5 )' \
	| sed 's/-[0-9]/1/g' > "${name}".aligned_merged.bed


####################################
##### first round of statistics ####
####################################


cat "${name}".aligned_merged.bed | awk '{OFS="\t"} { gsub(",","\t", $9)}1' | awk '{OFS="\t"} {sum=0; var=0 ; for (i=9; i<=NF; i++) {sum=sum+$i} mean=sum/(NF-8); for (i=9; i<=NF; i++) {var+=($i-mean)^2/(NF-9)} ; sd=sqrt(var) ; print $1,int($4-sd),int($5+sd),$6,"Cln"NR,$8,sd,$5-$4,$2,$3 }' | sed 's/-[0-9]/1/g' | bedtools sort | bedtools merge -s -c 4,5,6,7,8,9,10 -o sum,distinct,distinct,mean,mean,min,max -delim "" > "${name}".aligned_merged_stats.bed 
  
######################################
## coverage on mean or max locations #
######################################

awk '{OFS="\t"} {print $1,$9,$10,$4,$5,$6,$7,$8}' "${name}".aligned_merged_stats.bed \
	| bedtools coverage -counts -s -a - -b "${name}".aligned_TtoC_loc.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_TtoC.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_sorted.bed > "${name}".MinMax_Cov.bed 


awk '{OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8}' "${name}".aligned_merged_stats.bed \
	| bedtools coverage -counts -s -a - -b "${name}".aligned_TtoC_loc.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_TtoC.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_sorted.bed > "${name}".MeanStartMeanEnd_Cov.bed



##############################################################################
# Filter by distribution SDV of reads set > 68% for normal at sd-1 - sd +1g 
# Peaks that either will not follow normal distribution or the SDV is as
##############################################################################


paste "${name}".MinMax_Cov.bed "${name}".MeanStartMeanEnd_Cov.bed | awk -F "\t" '$11 > 0' > "${name}"_pasted.txt


awk -F "\t" '(($22/$11 >= 0.68) && ($7 < $8))' "${name}"_pasted.txt | awk -F "\t" '{OFS="\t"} \
  { if (($2 < $13) && ($3 > $14)) \
    { print $1,$13,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 > $14)) \
    { print $1,$2,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 < $14)) \
    { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 < $13) && ($3 < $14)) \
    { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11}
  else if (($2 == $13) && ($3 > $14)) \
    { print $1,$13,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 == $14)) \
    { print $1,$2,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 == $13) && ($3 < $14)) \
    { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 < $13) && ($3 == $14)) \
    { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
   else { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11} }' > "${name}".groups_normal.bed
    
awk '{OFS="\t"} {if (($22/$11 < 0.68) || ($7 > $8)) print $1,$2,$3,sqrt($4)/$4,$5,$6,$7,$8,$9,$10,$11}' ${name}_pasted.txt > "${name}".groups_NonNormal.bed 


rm ${name}_pasted.txt


#############################################
# subsample bam file for NonNormal clusters #
#############################################

number_of_NonNormal=`wc -l < "${name}".groups_NonNormal.bed`


samtools view -H "${name}".Aligned.sortedByCoord.out.bam > "${name}".aligned_sqrt.sam


bedtools intersect -abam "${name}".Aligned.sortedByCoord.out.bam -b "${name}".groups_NonNormal.bed > "${name}".aligned_NonNormal.bam

samtools index "${name}".aligned_NonNormal.bam


while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6; do samtools view -s $f4 "${name}".aligned_NonNormal.bam $f1":"$f2"-"$f3 >> "${name}".aligned_sqrt.sam ; done < "${name}".groups_NonNormal.bed 


samtools view -Shb "${name}".aligned_sqrt.sam --threads ${threads} \
	| samtools sort - --threads ${threads} | bedtools bamtobed -i - \
	| bedtools merge -s -i - -c 2,3,4,5,6,2 -o mean,mean,count,distinct,distinct,collapse \
	| awk -v var=${minRead} '($6 >= var ) && ($3-$2 > 1)' | sed 's/-[0-9]/1/g' > "${name}".aligned_merged_sqrt.bed



######################################################
# second round on subsamples sam round of statistics #
######################################################

cat "${name}".aligned_merged_sqrt.bed | awk '{OFS="\t"} { gsub(",","\t", $9)}1' \
	| awk '{OFS="\t"} {sum=0; var=0 ; for (i=9; i<=NF; i++) {sum=sum+$i} mean=sum/(NF-8); for (i=9; i<=NF; i++) {var+=($i-mean)^2/(NF-9)} ; sd=sqrt(var) ; print $1,int($4-sd),int($5+sd),$6,"Cla"NR,$8,sd,$5-$4,$2,$3 }' \
	| sed 's/-[0-9]/1/g' | bedtools sort | bedtools merge -s -c 4,5,6,7,8,9,10 -o sum,distinct,distinct,mean,mean,min,max -delim "" > "${name}".aligned_merged_stats_sqrt.bed 
  
awk '{OFS="\t"} {print $1,$9,$10,$4,$5,$6,$7,$8}' "${name}".aligned_merged_stats_sqrt.bed \
	| bedtools coverage -counts -s -a - -b "${name}".aligned_TtoC_loc.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_TtoC.bed | bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_sorted.bed > "${name}".MinMax_Cov_sqrt.bed 

awk '{OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8}' "${name}".aligned_merged_stats_sqrt.bed \
	| bedtools coverage -counts -s -a - -b "${name}".aligned_TtoC_loc.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_TtoC.bed \
	| bedtools coverage -counts -s -F 0.66 -a - -b "${name}".aligned_sorted.bed > "${name}".MeanStartMeanEnd_Cov_sqrt.bed

paste "${name}".MinMax_Cov_sqrt.bed "${name}".MeanStartMeanEnd_Cov_sqrt.bed | awk -F "\t" '$11 > 0' > "${name}"_pasted_sqrt.txt


awk -F "\t" '(($22/$11 >= 0.68) && ($7 < $8))' "${name}"_pasted_sqrt.txt |awk -F "\t" '{OFS="\t"} \
   { if (($2 < $13) && ($3 > $14)) \
    { print $1,$13,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 > $14)) \
    { print $1,$2,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 < $14)) \
    { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 < $13) && ($3 < $14)) \
    { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11}
  else if (($2 == $13) && ($3 > $14)) \
    { print $1,$13,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 > $13) && ($3 == $14)) \
    { print $1,$2,$14,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 == $13) && ($3 < $14)) \
    { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
  else if (($2 < $13) && ($3 == $14)) \
    { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11} \
   else { print $1,$13,$3,$4,$5,$6,$7,$8,$9,$10,$11} }' > "${name}".groups_normal_sqrt.bed
   
#################################################################################################################
# Merge normal stage 1 and normal stage 2 and filter satellites by tukeys factor is set by nt from peak to peak #
#################################################################################################################

factor=50

cat "${name}".groups_normal.bed "${name}".groups_normal_sqrt.bed | bedtools sort | awk -F "\t" -v var=$factor '{if (NR == 1) {prev1=$1;prev3=$3;prev6=$6;prev11=$11; {OFS="\t";print $0}} \
  else if ( $1 != prev1 ||$6 != prev6 || sqrt(($2-prev3)^2) > var*log((prev11/$11))/log(10) ) {OFS="\t";print $0} \
  ;prev1=$1;prev3=$3;prev6=$6;prev11=$11 }' | tac | awk -F "\t" -v var=$factor '{if (NR == 1) {prev1=$1;prev3=$3;prev6=$6;prev11=$11; {OFS="\t";print $0}} \
  else if ( $1 != prev1 ||$6 != prev6 || sqrt(($2-prev3)^2) > var*log((prev11/$11))/log(10) ) {OFS="\t";print $0} \
  ;prev1=$1;prev3=$3;prev6=$6;prev11=$11 }' | tac | awk -F "\t" -v var=$factor '{if (NR == 1) {prev1=$1;prev3=$3;prev6=$6;prev11=$11; {OFS="\t";print $0}} \
  else if ( $1 != prev1 ||$6 != prev6 || sqrt(($2-prev3)^2) > var*log((prev11/$11))/log(10) ) {OFS="\t";print $0} \
  ;prev1=$1;prev3=$3;prev6=$6;prev11=$11 }' | tac | awk -F "\t" -v var=$factor '{if (NR == 1) {prev1=$1;prev3=$3;prev6=$6;prev11=$11; {OFS="\t";print $0}} \
  else if ( $1 != prev1 ||$6 != prev6 || sqrt(($2-prev3)^2) > var*log((prev11/$11))/log(10) ) {OFS="\t";print $0} \
  ;prev1=$1;prev3=$3;prev6=$6;prev11=$11 }' | tac | awk -F "\t" -v var=$factor '{if (NR == 1) {prev1=$1;prev3=$3;prev6=$6;prev11=$11; {OFS="\t";print $0}} \
  else if ( $1 != prev1 ||$6 != prev6 || sqrt(($2-prev3)^2) > var*log((prev11/$11))/log(10) ) {OFS="\t";print $0} \
  ;prev1=$1;prev3=$3;prev6=$6;prev11=$11 }' | tac | awk -F "\t" -v var=$factor '{if (NR == 1) {prev1=$1;prev3=$3;prev6=$6;prev11=$11; {OFS="\t";print $0}} \
  else if ( $1 != prev1 ||$6 != prev6 || sqrt(($2-prev3)^2) > var*log((prev11/$11))/log(10) ) {OFS="\t";print $0} \
  ;prev1=$1;prev3=$3;prev6=$6;prev11=$11 }' | tac | awk -F "\t" -v var=$factor '{OFS="\t"; print $1,$2,$3,$5,".",$6,$7,$8,$9,$10,$11}' > "${name}"_filt.groups
 
 
 
###########################################
# annotate peaks using bedtools intersect #
###########################################

bedtools intersect -s -wa -wb -loj -a "${name}"_filt.groups -b "${annotation}" \
	| awk -F "\t" '{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,($10/$11),$15,$16,$18}' \
	| bedtools sort \
	| bedtools merge -s -c 4,5,6,7,8,9,10,11,12,13,14,15 \
	-o distinct,distinct,distinct,max,max,max,max,max,max,distinct,distinct,distinct -delim "|" \
	| sed -e 's/CDS|exon/CDS/g' -e 's/5utr|exon/5utr/g' -e 's/3utr|exon/3utr/g' > "${name}".groups.tsv


cat "${name}".groups.tsv | awk -F "\t" '{OFS="\t"; print $1,$2,$3,$4,$5,$6 }' \
	| bedtools getfasta -s -fi ${genomefasta} -bed - -tab > "${name}".groups.fa


echo "chromosome,start,end,cluster_id,dot,strand,gene_name,annotation,gene_type,sequence,crosslinked_reads,total_reads,fraction,TtoC_locations,avr_readlength,read_distance_SD,chr:start-end(strand)" > "${name}".groups.csv


paste "${name}".groups.tsv "${name}".groups.fa | awk '{OFS=","; print $1,$2,$3,$4,$5,$6,$13,$14,$15,$17,$10,$11,$12,$9,$8,$7,$16 }' \
	| sed -e 's/CDS|exon/CDS/g' -e 's/5utr|exon/5utr/g' -e 's/3utr|exon/3utr/g' -e 's/,\.,\.,\.,/,NA,No_Annotation,NA,/g' -e 's/protein_coding|protein_coding/protein_coding/g' >> "${name}".groups.csv


cat "${name}".groups.csv | awk -F "," -v convloc="${convloc}" -v frac="${fraction}" '{ if ($14 >= convloc && $13 >=frac) print $0 } ' \
	| awk -F, 'BEGIN {
  IGNORECASE=1
# Priority for annotation (column 8) — exon moved to highest
  priority_annot[1]="exon"; priority_annot[2]="3ss"; priority_annot[3]="5ss"; 
  priority_annot[4]="start_codon"; priority_annot[5]="stop_codon"; 
  priority_annot[6]="3utr"; priority_annot[7]="5utr"; 
  priority_annot[8]="CDS"; priority_annot[9]="intron";

  # Priority for biotype (column 9)
  priority_biotype[1]="snoRNA"; priority_biotype[2]="snRNA"; priority_biotype[3]="rRNA";
  priority_biotype[4]="miRNA"; priority_biotype[5]="tRNA"; priority_biotype[6]="Mt_tRNA";
  priority_biotype[7]="scaRNA"; priority_biotype[8]="protein_coding"; priority_biotype[9]="lncRNA";

  # Prefix priority for gene names
  gene_prefix[1]="SNOR"; gene_prefix[2]="RN"; gene_prefix[3]="RNA";
  gene_prefix[4]="MIR"; gene_prefix[5]="tRNA"; gene_prefix[6]="Mt_tRNA"; gene_prefix[7]="SCARN";
}
{
  if(NR==1) {print; next} # Print header line unchanged

  # --- Column 7: gene name ---
  if (index($7, "|") > 0) {
    split($7, genes, "|");
    best_gene = genes[1];
    found = 0;

    for (p in gene_prefix) gene_prefix_lc[p] = tolower(gene_prefix[p]);

    for (p in gene_prefix) {
      for (g in genes) {
        if (index(tolower(genes[g]), gene_prefix_lc[p]) == 1) {
          best_gene = genes[g];
          found = 1;
          break;
        }
      }
      if (found) break;
    }

    if (!found) {
      for (g in genes) {
        if (index(genes[g], "-") == 0) {
          best_gene = genes[g];
          break;
        }
      }
    }

    $7 = best_gene;
  }

  # --- Column 8: annotation ---
  if ($8 != "No_Annotation") {
    for (i=1; i<=9; i++) {
      if (index($8, priority_annot[i]) > 0) {
        $8 = priority_annot[i];
        break;
      }
    }
  }

  # --- Column 9: biotype ---
  if (index($9, "|") > 0) {
    split($9, types, "|");
    found_bio = 0;
    for (i=1; i<=9 && !found_bio; i++) {
      for (t in types) {
        if (types[t] == priority_biotype[i]) {
          $9 = priority_biotype[i];
          found_bio = 1;
          break;
        }
      }
    }
    if (!found_bio) $9 = types[1];
  }

  # Print modified line
  for(i=1; i<NF; i++) printf "%s,", $i;
  printf "%s\n", $NF;
}' | sed 's/3utr,tRNA/exon,tRNA/g' | sed 's/3utr,rRNA/exon,rRNA/g' | sed 's/5utr,tRNA/exon,tRNA/g' | sed 's/5utr,rRNA/exon,rRNA/g' | sed 's/3utr,snoRNA/exon,snoRNA/g' | sed 's/3utr,snRNA/exon,snRNA/g' > "${name}".clusters.csv


grep -v chromosome "${name}".clusters.csv \
	| awk -F "," '{OFS="\t"; print $1,$2,$3,$5,$4,$6 }' \
	| bedtools sort > "${name}".clusters.bed


echo "Binding statistics" >> "${name}"_stats.txt

groups=$(wc -l < "${name}".groups.csv)
echo "groups: ${groups}" >> "${name}"_stats.txt
clusters=$(wc -l < "${name}".clusters.csv)
echo "clusters: ${clusters}" >> "${name}"_stats.txt

Protein_coding_5UTR_XL_reads=$(grep 5utr < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "Protein_coding-5UTR_XL_reads: ${Protein_coding_5UTR_XL_reads}" >> "${name}"_stats.txt
Protein_coding_CDS_XL_reads=$(grep CDS < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "Protein_coding-CDS_XL_reads: ${Protein_coding_CDS_XL_reads}" >> "${name}"_stats.txt
Protein_coding_intron_XL_reads=$(grep intron < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "Protein_coding-introns_XL_reads: ${Protein_coding_intron_XL_reads}" >> "${name}"_stats.txt
Protein_coding_3UTR_XL_reads=$(grep 3utr < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "Protein_coding-3UTR_XL_reads: ${Protein_coding_3UTR_XL_reads}" >> "${name}"_stats.txt
lincRNA_XL_reads=$(grep lincRNA < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "lincRNA_XL_reads: ${lincRNA_XL_reads}" >> "${name}"_stats.txt
snoRNA_XL_reads=$(grep snoRNA < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "snoRNA_XL_reads: ${snoRNA_XL_reads}" >> "${name}"_stats.txt
snRNA_XL_reads=$(grep snRNA < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "snRNA_XL_reads: ${snRNA_XL_reads}" >> "${name}"_stats.txt
rRNA_XL_reads=$(grep rRNA < "${name}".clusters.csv | awk -F',' '{sum+=$11;}END{print sum;}')
echo "rRNA_XL_reads: ${rRNA_XL_reads}" >> "${name}"_stats.txt

Protein_coding_5UTR_clusters=$(grep 5utr < "${name}".clusters.csv | wc -l)
echo "Protein_coding-5UTR_clusters: ${Protein_coding_5UTR_clusters}" >> "${name}"_stats.txt
Protein_coding_CDS_clusters=$(grep CDS < "${name}".clusters.csv | wc -l)
echo "Protein_coding-CDS_clusters: ${Protein_coding_CDS_clusters}" >> "${name}"_stats.txt
Protein_coding_intron_clusters=$(grep intron < "${name}".clusters.csv | wc -l)
echo "Protein_coding-introns_clusters: ${Protein_coding_intron_clusters}" >> "${name}"_stats.txt
Protein_coding_3UTR_clusters=$(grep 3utr < "${name}".clusters.csv | wc -l)
echo "Protein_coding-3UTR_clusters: ${Protein_coding_3UTR_clusters}" >> "${name}"_stats.txt
lncRNA_clusters=$(grep lncRNA < "${name}".clusters.csv | wc -l)
echo "lncRNA_clusters: ${lncRNA_clusters}" >> "${name}"_stats.txt
snoRNA_clusters=$(grep snoRNA < "${name}".clusters.csv | wc -l)
echo "snRNA_clusters: ${snRNA_clusters}" >> "${name}"_stats.txt
snRNA_clusters=$(grep snRNA < "${name}".clusters.csv | wc -l)
echo "snoRNA_clusters: ${snoRNA_clusters}" >> "${name}"_stats.txt
rRNA_clusters=$(grep rRNA < "${name}".clusters.csv | wc -l)
echo "rRNA_clusters: ${rRNA_clusters}" >> "${name}"_stats.txt



###################
## cleaning files #
###################

rm "${name}"_collapsed.fasta
#rm "${name}"_trimmed.fasta
rm "${name}".SJ.out.tab
rm "${name}".Log.progress.out
rm "${name}".Log.final.out
rm "${name}".Log.out
rm "${name}".aligned_sorted.bed
rm "${name}".aligned_TtoC.bed
rm "${name}".aligned_merged.bed
rm "${name}".aligned_merged_stats.bed
rm "${name}".MinMax_Cov.bed
rm "${name}".MeanStartMeanEnd_Cov.bed
rm "${name}".groups_normal.bed
rm "${name}".aligned_sqrt.sam
rm "${name}".aligned_merged_sqrt.bed
rm "${name}".aligned_merged_stats_sqrt.bed
rm "${name}".MinMax_Cov_sqrt.bed
rm "${name}".MeanStartMeanEnd_Cov_sqrt.bed
rm "${name}"_pasted_sqrt.txt
rm "${name}".groups_normal_sqrt.bed
rm "${name}"_filt.groups
rm "${name}".log
rm "${name}".groups.tsv
rm "${name}".groups.fa
rm "${name}".aligned_NonNormal.bam
rm "${name}".aligned_NonNormal.bam.bai
rm "${name}".groups_NonNormal.bed 


###########################
# get the 5mer enrichment #
###########################

awk -F "\t" '{OFS="\t"; print $1,$2-200,$3-200,$4,$5,$6}' "${name}".clusters.bed \
	| awk '$2 > 0' > background_"${name}".clusters.bed

bedtools getfasta -s -fi "${genomefasta}" -bed "${name}".clusters.bed -fo bound_"${name}".fasta
bedtools getfasta -s -fi "${genomefasta}" -bed background_"${name}".clusters.bed -fo background_"${name}".fasta

##############################################
#create all combinations of 5mers and 3x_3mers
##############################################

perl -le '@c = ("A","T","C","G"); for $a (@c){for $b(@c){for $c(@c){for $d(@c){for $e(@c){print "$a$b$c$d$e"}}}}}' > "${name}"_5mers.txt
perl -le '@c = ("A","T","C","G"); for $a (@c){for $b(@c){for $c(@c){for $d(@c){for $e (@c){for $f(@c){for $g(@c){for $h(@c){for $i(@c){print "$a$b$c\[\A\G\T\C\]\\{0,15\\}$d$e$f\[\A\G\T\C\]\\{0,15\\}$g$h$i"}}}}}}}}}' > ${input}_tri_3mers.txt


#count each mer in the fasta files

grep -v ">" background_"${name}".fasta | grep -v "@" > filt_background_"${name}".fasta

grep -v ">" bound_"${name}".fasta | grep -v "@" > filt_bound_"${name}".fasta

rm background_"${name}".fasta
rm bound_"${name}".fasta


#parallel 'echo {1} ; grep -i -o {1} {2} | wc -l' :::: "${name}"_5mers.txt ::: filt_background_"${name}".fasta \
#	| paste - - | sort > background_"${name}"_5mers_out1.txt
#wait
#parallel 'echo {1} ; grep -i -o {1} {2} | wc -l' :::: "${name}"_5mers.txt ::: filt_bound_"${name}".fasta \
#	| paste - - | sort > bound_"${name}"_5mers_out1.txt
#wait




while read mer; do
    echo "$mer"
    grep -i -o "$mer" filt_background_"${name}".fasta | wc -l
done < "${name}"_5mers.txt | paste - - | sort > background_"${name}"_5mers_out1.txt



while read mer; do
    echo "$mer"
    grep -i -o "$mer" filt_bound_"${name}".fasta | wc -l
done < "${name}"_5mers.txt | paste - - | sort > bound_"${name}"_5mers_out1.txt




while read mer; do
    echo "$mer"
    grep -i -o "$mer" filt_background_"${name}".fasta | wc -l
done < ${input}_tri_3mers.txt | paste - - | sort > background_"${name}"_tri_3mers_out1.txt



while read mer; do
    echo "$mer"
    grep -i -o "$mer" filt_bound_"${name}".fasta | wc -l
done < ${input}_tri_3mers.txt | paste - - | sort > bound_"${name}"_tri_3mers_out1.txt




rm filt_background_"${name}".fasta
rm filt_bound_"${name}".fasta


paste background_"${name}"_5mers_out1.txt bound_"${name}"_5mers_out1.txt | awk '{OFS="\t"; print $1,$2,$4}' > "${name}"_5mers.txt
paste background_"${name}"_tri_3mers_out1.txt bound_"${name}"_tri_3mers_out1.txt | awk '{OFS="\t"; print $1,$2,$4}' > "${name}"_tri_3mers.txt

rm background_"${name}"_5mers_out1.txt
rm bound_"${name}"_5mers_out1.txt

rm background_"${name}"_tri_3mers_out1.txt
rm bound_"${name}"_tri_3mers_out1.txt


rm "${name}"_5mers.txt
rm "${name}"_tri_3mers.txt


#############
## THE END ##
#############
