######################################################

## Tools and scripts
```
liftOver=/srv/gevaertlab/tools/UCSCtools/liftOver
hg19ToHg38overchain=/srv/gevaertlab/tools/UCSCtools/hg19ToHg38.over.chain  
```

- mergeBed and intersectBed are commands from bedtools. 
- Extract.pl: in-house script. Refer to zhengh42/LittlestScripts

######################################################

## Download 
	
### GENCODE
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
gzip -d gencode.v27.annotation.gtf.gz
```

### NONCODE
```
wget http://www.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz
gzip -d NONCODEv5_human_hg38_lncRNA.gtf.gz
```

### MiTranscriptome
```
wget http://mitranscriptome.org/download/mitranscriptome.gtf.tar.gz
tar -zxf mitranscriptome.gtf.tar.gz
gzip -d mitranscriptome.gtf/mitranscriptome.v2.gtf.gz
```

######################################################

## Stats

```
less gencode.v27.annotation.gtf | grep -v ^# | awk '$3~/transcript/' | cut -f9 | awk 'OFS="\t"{print $2,$4,$6,$12}' | sed 's/[";]//g' | awk 'OFS="\t"{$5="other";$6="other";if($4~/non_coding|3prime_overlapping_ncRNA|antisense|bidirectional_promoter_lncRNA|lincRNA|macro_lncRNA|sense_intronic|sense_overlapping/)$6="lncRNA";if($4~/protein_coding/)$6="proteincoding";if($3~/protein_coding/)$5="proteincoding";if($3~/non_coding|3prime_overlapping_ncRNA|antisense|bidirectional_promoter_lncRNA|lincRNA|macro_lncRNA|sense_intronic|sense_overlapping/)$5="lncRNA";  print $0}'| sed '1i gene_id\ttranscript\tgenetype1\ttranscripttype1\tgenetype\ttranscripttype' > gencode.v27.gene_transcript_type
less gencode.v27.gene_transcript_type | cut -f1,5 | sed 1d | sort | uniq | cut -f2 | sort  | uniq -c # gene type in GENCODE
less NONCODEv5_human_hg38_lncRNA.gtf | awk '{print $10}' | sort | uniq | wc -l # gene type in NONCODE
les mitranscriptome.gtf/mitranscriptome.v2.gtf | awk '$3~/transcript/' | awk '{print $10,$12}' | sort | uniq | awk '{print $1}' | sort | uniq -c # gene type in MiTranscriptome
```

**Composition of each transcriptome source**

|      | GENCODE  |      NONCDOE      |  MiTranscriptome |
|:----------:|:----------:|:-------------:|:------:|
| Protein-coding |  19,836 | - | 19,922 |
| lncRNA |  14,168 | 96,308 | 63,615 |
| Other |  24,284 | - | 15,445 |


######################################################

## Get bed file for lncRNAs transcripts

### GENCODE
```
less gencode.v27.annotation.gtf | grep -v ^# | awk '$3!="gene"' | awk '{print $12,$0}' | sed 's/^"//;s/";\s\+/\t/' > gencode.v27.annotation.gtf.tmp
less gencode.v27.gene_transcript_type | awk '$5=="lncRNA"' | sort | uniq | Extract.pl - 2 gencode.v27.annotation.gtf.tmp 1 | awk '$4=="transcript"' | cut -f2,5,6 | sort -k1,1 -k2,2n  | mergeBed -i - > gencode.v27.annotation.lncRNA.bed
```

### NONCODE
```
less NONCODEv5_human_hg38_lncRNA.gtf  | grep -v ^# | awk '{print $12,$0}' | sed 's/^"//;s/";\s\+/\t/' > NONCODEv5_human_hg38_lncRNA.gtf.tmp
less NONCODEv5_human_hg38_lncRNA.gtf.tmp | awk '$4=="transcript"' | cut -f2,5,6 | sort -k1,1 -k2,2n  | mergeBed -i - > NONCODEv5_human_hg38.lncRNA.bed
```
### MiTranscriptome
```
$liftOver -gff mitranscriptome.gtf/mitranscriptome.v2.gtf $hg19ToHg38overchain mitranscriptome.gtf/mitranscriptome.v2.hg38.gtf mitranscriptome.gtf/mitranscriptome.v2.hg19tohg38.unmapped # MiTranscriptome is based on hg19, so before comparison it is lifted over to hg38
less mitranscriptome.gtf/mitranscriptome.v2.hg38.gtf | awk '$3=="transcript" && $10~"lncrna"' | cut -f1,4,5 | sort -k1,1 -k2,2n |   mergeBed -i - > mitranscriptome.v2.hg38.lncRNA.bed
```

######################################################

## Comparison

### Total length of lncRNA transcripts in base
```
awk '{print $3-$2}' gencode.v27.annotation.lncRNA.bed | awk 'BEGIN{sum=0}{sum+=$1}END{print sum/1000000}' 
awk '{print $3-$2}' NONCODEv5_human_hg38.lncRNA.bed | awk 'BEGIN{sum=0}{sum+=$1}END{print sum/1000000}'
awk '{print $3-$2}' mitranscriptome.v2.hg38.lncRNA.bed | awk 'BEGIN{sum=0}{sum+=$1}END{print sum/1000000}'
```

### Total length of overlapping lncRNA transcripts 
```
intersectBed -a gencode.v27.annotation.lncRNA.bed -b NONCODEv5_human_hg38.lncRNA.bed | awk '{print $3-$2}' | awk 'BEGIN{sum=0}{sum+=$1}END{print sum/1000000}'
intersectBed -a gencode.v27.annotation.lncRNA.bed -b mitranscriptome.v2.hg38.lncRNA.bed | awk '{print $3-$2}' | awk 'BEGIN{sum=0}{sum+=$1}END{print sum/1000000}'
intersectBed -a NONCODEv5_human_hg38.lncRNA.bed -b mitranscriptome.v2.hg38.lncRNA.bed | awk '{print $3-$2}' | awk 'BEGIN{sum=0}{sum+=$1}END{print sum/1000000}'
```

######################################################

## Summarize

|      | GENCODE  |      NONCDOE      |  MiTranscriptome |
|:----------:|:----------:|:-------------:|:------:|
| GENCODE |  331 | 325 | 268 |
| NONCDOE |       |   1028 |  641|
| MiTranscriptome |  | | 1201|

*The number indicates the total length (Mb) of overlapping lncRNA transcripts bewteen two source.*
