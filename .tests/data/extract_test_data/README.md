# Make big test dataset

Extract region and reads for big test dataset.

```bash
export PATH="~/programs/bedtools-2.29.2/bin:$PATH"
```

Get SNP density across 10K windows
```bash
D="/scratch/timothy/projects/0022_Coral_Genotype_Analysis/03_Analysis/2022-02-05/Pocillopora_acuta"

zcat $D/03_polulation_structure/vcftools_relatedness2/GVCFall.filtered.recode.vcf.gz | grep -v '^#' | awk -F'\t' '{print $1"\t"$2-1"\t"$2}' | sort -k1,1 -k2,2n -k3,3n > SNPs.bed
awk -F'\t' '{print $1"\t0\t"$2}' $D/00_databases/Pocillopora_acuta_KBHIv2.assembly.fasta.fai > genome.txt

bedtools makewindows -b genome.txt -w 10000 -s 1000 > windows.bed

bedtools coverage -a windows.bed -b SNPs.bed > overlap.bed

sort -k4,4nr overlap.bed | less
```

Make output data directory.
```bash
mkdir -p ../data
```

Region with highest SNP density: `Pocillopora_acuta_KBHIv2___xfSc0000005:600000-800000`
```bash
bedtools getfasta -fi $D/00_databases/Pocillopora_acuta_KBHIv2.assembly.fasta -bed <(echo -e 'Pocillopora_acuta_KBHIv2___xfSc0000005\t600000\t800000') | sed -e 's/>.*/>Contig1/' > ../data/ref.fa
~/programs/bbmap/bbmap.sh ref=../data/ref.fa
```

Get annotations
```bash
cut -f5-11 $D/samples_Pacuta.annotations.txt | sed -e 's/plugid/id/' > ../data/samples.tsv
```

Extract reads from raw data on our server
```bash
while read ID;
do
  R1=$(ls /scratch/timothy/TranscriptomeAnalysis_Mcapitata_Pacuta_TimeSeriesDE/RNA_Seq_data/raw_reads_renamed/*${ID}_R1.fastq.gz)
  R2=$(ls /scratch/timothy/TranscriptomeAnalysis_Mcapitata_Pacuta_TimeSeriesDE/RNA_Seq_data/raw_reads_renamed/*${ID}_R2.fastq.gz)
  echo -e "${ID}\t${R1}\t${R2}"
done < <(awk 'NR>1{print $1}' samples.tsv) > reads.tsv
```




```bash
bedtools getfasta -fi $D/00_databases/Pocillopora_acuta_KBHIv2.assembly.fasta -bed <(echo -e 'Pocillopora_acuta_KBHIv2___xfSc0000005\t100000\t300000') | gzip -c > ../ref_bad.fa.gz
```



