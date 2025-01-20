
cat ../../samples.txt | while read line
do
bwa mem -t 8 /path/to/reference/merged.fna ../../03decom/$line'_paired_1'.fastq ../../03decom/$line'_paired_2'.fastq > 01bwa/$line.sam
samtools view --threads 8 -b 01bwa/$line.sam > 02bam/$line.raw.bam

time picard AddOrReplaceReadGroups \
    -I 02bam/$line.raw.bam \
    -O 02bam/$line.AddOrReplaceReadGroups.bam \
    -SO "coordinate" \
    --RGID foo \
    --RGLB bar \
    --RGPL illumina \
    --RGSM $line \
    --RGPU lib1 

time picard MarkDuplicates \
    -I 02bam/$line.AddOrReplaceReadGroups.bam \
    -O 02bam/$line.MarkDuplicates.bam \
    -M 02bam/$line.MarkDuplicates.matrics \
    --VALIDATION_STRINGENCY LENIENT 

time picard SortSam \
    -I 02bam/$line.MarkDuplicates.bam \
    -O 02bam/$line.SortSam.bam \
    -SO "coordinate" \
    --CREATE_INDEX True 
    
time gatk --java-options "-Xmx32g" HaplotypeCaller  \
   -R /path/to/reference/merged.fna \
   -I 02bam/$line.SortSam.bam \
   --indel-size-to-eliminate-in-ref-model 50 \
   --native-pair-hmm-threads 8 \
   --sample-ploidy 1 \
   -O 03vcf/$line.raw.vcf 

time gatk --java-options "-Xmx32g" SelectVariants \
	-R /path/to/reference/merged.fna \
	-V 03vcf/$line.raw.vcf \
    --max-indel-size 50 \
	-select-type INDEL \
	-O 04indel/$line.raw.indel.vcf

done
