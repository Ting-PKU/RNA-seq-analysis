## 1. 去接头
```
cat config |while read line;do arr=($line); fq1=${arr[0]}; fq2=${arr[1]};trim_galore -q 25 --phred33  -e 0.1 --stringency 2 --paired $fq1 $fq2 & done
```
## 2. 准备 star genome index
```
STAR --runThreadN 38 \
--runMode genomeGenerate --genomeDir star_genome \
--genomeFastaFiles /data/ting/genome/hg38/hg38.fa \
--sjdbGTFfile /data/ting/genome/hg38/hg38.refGene.gtf \
--sjdbOverhang 149 &
```
## 3. 比对
```
cat config | while read line; do array=($line);sample=${array[0]};fq1=${array[1]};fq2=${array[2]};\
nohup STAR --runThreadN 10 \
--genomeDir ~/genome/star_genome_v2.7.11a/ --readFilesCommand zcat \
--readFilesIn $fq1 $fq2 --outSAMtype BAM Unsorted \
--outFileNamePrefix ./align_results/$sample & done
```
## 4. htseq-count
```
ls *.bam | while read id; do htseq-count  -q -f bam -r name -s no $id \
/data/ting/Annotation/Genes/genes.gtf > ${id%%.*}.count  & done
```
