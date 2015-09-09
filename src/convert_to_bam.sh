#! /bin/sh

SAM_PATH=$1

echo "SAMtools view..."
echo "samtools view -Sb $SAM_PATH.sam > $SAM_PATH.bam"
samtools view -Sb $SAM_PATH.sam > $SAM_PATH.bam
echo " "

echo "SAMtools sort..."
echo "samtools sort $SAM_PATH.bam $SAM_PATH-sorted"
samtools sort $SAM_PATH.bam $SAM_PATH-sorted
echo " "

echo "SAMtools index..."
echo "samtools index $SAM_PATH-sorted.bam $SAM_PATH-sorted.bam.bai"
samtools index $SAM_PATH-sorted.bam $SAM_PATH-sorted.bam.bai
echo " "

echo "SAMtools flagstat..."
echo "samtools flagstat $SORTED_BAM_FILE"
samtools flagstat $SAM_PATH-sorted.bam
echo " "
