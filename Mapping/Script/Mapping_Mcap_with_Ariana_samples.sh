

cd /lustre1/home/mass/amalia/montipora_mapping/

mkdir Ariana_samples

cd /lustre1/home/mass/amalia/montipora_mapping/Ariana_samples

#Making an executable script in the interactive session
nano download_mcap.sh

#!/bin/bash

# List of SRR accessions
SRR_LIST=(
SRR14864072
SRR14864071
SRR14864070
SRR14864069
SRR14864068
SRR14864067
SRR14864066
SRR14864065
SRR14864064
)

# Output directory for FASTQ files
OUTDIR="/lustre1/home/mass/amalia/montipora_mapping/fastqreads"
mkdir -p "$OUTDIR"

# Download + Convert loop
for SRR in "${SRR_LIST[@]}"; do
    echo " Downloading $SRR..."
    prefetch "$SRR"

    echo "Converting $SRR to FASTQ..."
    fasterq-dump "$SRR" --threads 4 --progress --outdir "$OUTDIR"
done

echo "All downloads and conversions completed."

chmod +x download_mcap.sh
conda activate Mcap #if not already activated
./download_mcap.sh

#When done
cd fastqreads
md5sum *.fastq > raw_checksum.md5

md5sum -c raw_checksum.md5 #check it's OK for all files
#All OK
SRR14864064_1.fastq: OK
SRR14864064_2.fastq: OK
SRR14864065_1.fastq: OK
SRR14864065_2.fastq: OK
SRR14864066_1.fastq: OK
SRR14864066_2.fastq: OK
SRR14864067_1.fastq: OK
SRR14864067_2.fastq: OK
SRR14864068_1.fastq: OK
SRR14864068_2.fastq: OK
SRR14864069_1.fastq: OK
SRR14864069_2.fastq: OK
SRR14864070_1.fastq: OK
SRR14864070_2.fastq: OK
SRR14864071_1.fastq: OK
SRR14864071_2.fastq: OK
SRR14864072_1.fastq: OK
SRR14864072_2.fastq: OK
SRR22293445.fastq: OK
SRR22293446.fastq: OK
SRR22293447.fastq: OK
SRR22293448.fastq: OK
SRR22293449.fastq: OK
SRR22293450.fastq: OK
SRR22293451.fastq: OK
SRR22293452.fastq: OK
SRR22293453.fastq: OK
SRR22293454.fastq: OK
SRR22293455.fastq: OK
SRR22293456.fastq: OK
SRR22293457.fastq: OK
SRR22293458.fastq: OK
SRR22293459.fastq: OK
SRR22293460.fastq: OK
SRR22293461.fastq: OK
SRR22293462.fastq: OK
SRR22293463.fastq: OK
SRR22293464.fastq: OK
SRR22293465.fastq: OK
SRR22293466.fastq: OK
SRR22293467.fastq: OK
SRR22293468.fastq: OK
SRR22293469.fastq: OK
SRR22293470.fastq: OK
SRR22293471.fastq: OK
SRR22293472.fastq: OK
SRR22293473.fastq: OK
SRR22293474.fastq: OK
SRR22293475.fastq: OK
SRR22293476.fastq: OK
SRR22293477.fastq: OK
SRR22293478.fastq: OK
SRR22293479.fastq: OK
SRR22293480.fastq: OK
SRR22293481.fastq: OK
SRR22293482.fastq: OK
SRR22293483.fastq: OK
zgrep -c "@SRR" *fastq #Count number of reads per file


# Navigate to the Ariana_samples folder
cd /lustre1/home/mass/amalia/montipora_mapping/Ariana_samples


#Symbolic links
cd /lustre1/home/mass/amalia/montipora_mapping/Ariana_samples
rm -rf raw && mkdir raw
for f in fastqreads/*.fastq; do
  base=$(basename "$f")
  command ln -s "$(pwd)/$f" "raw/$base"
done

cd /lustre1/home/mass/amalia/montipora_mapping/Ariana_samples

# sanity: how many FASTQs?
ls raw/*.fastq | wc -l    # expect 57

# run FastQC on raw reads
mkdir -p fastqc_raw
fastqc raw/*.fastq --outdir fastqc_raw --threads 8

# summarize with MultiQC (pre-trim only)
mkdir -p outputs/multiqc_fastqc_raw
multiqc fastqc_raw -o outputs/multiqc_fastqc_raw -n multiqc_fastqc_raw.html

cd /lustre1/home/mass/amalia/montipora_mapping/Ariana_samples
multiqc fastqc_raw -o multiqc_fastqc_raw



#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

cd /lustre1/home/mass/amalia/montipora_mapping/Ariana_samples

# make sure fastp is available (uncomment if needed)
# conda activate Mcap

THREADS=8
mkdir -p trimmed_fastq fastp_reports

# ---- paired-end ----
for r1 in raw/*_1.fastq; do
  sample=${r1##*/}; sample=${sample%_1.fastq}
  r2="raw/${sample}_2.fastq"
  # skip if already trimmed
  if [[ -f "trimmed_fastq/${sample}_1.trim.fastq.gz" && -f "trimmed_fastq/${sample}_2.trim.fastq.gz" ]]; then
    echo "skip (PE exists): $sample"; continue
  fi
  if [[ -f "$r2" ]]; then
    echo "fastp (PE): $sample"
    fastp \
      -i "$r1" -I "$r2" \
      -o "trimmed_fastq/${sample}_1.trim.fastq.gz" \
      -O "trimmed_fastq/${sample}_2.trim.fastq.gz" \
      --thread "$THREADS" \
      --detect_adapter_for_pe \
      --html "fastp_reports/${sample}_fastp.html" \
      --json "fastp_reports/${sample}_fastp.json" \
      --report_title "$sample fastp report"
  fi
done

# ---- single-end ----
for se in raw/*.fastq; do
  b=${se##*/}
  [[ "$b" == *_1.fastq || "$b" == *_2.fastq ]] && continue
  sample=${b%.fastq}
  # skip if already trimmed
  if [[ -f "trimmed_fastq/${sample}.trim.fastq.gz" ]]; then
    echo "skip (SE exists): $sample"; continue
  fi
  echo "fastp (SE): $sample"
  fastp \
    -i "$se" \
    -o "trimmed_fastq/${sample}.trim.fastq.gz" \
    --thread "$THREADS" \
    --html "fastp_reports/${sample}_fastp.html" \
    --json "fastp_reports/${sample}_fastp.json" \
    --report_title "$sample fastp report"
done

# optional summary (if installed)
if command -v multiqc >/dev/null 2>&1; then
  multiqc fastp_reports -o fastp_reports
fi

echo "fastp finished. Outputs: trimmed_fastq/ ; Reports: fastp_reports/"

#Downloading reports: 
scp amalia@hive02.haifa.ac.il:/lustre1/home/mass/amalia/montipora_mapping/Ariana_samples/multiqc_fastqc_raw/multiqc_report.html "C:\Users\amurg\Downloads\multiqc_fastqc_raw.html"
scp amalia@hive02.haifa.ac.il:/lustre1/home/mass/amalia/montipora_mapping/Ariana_samples/fastp_reports/multiqc_report.html "C:\Users\amurg\Downloads\multiqc_fastp.html"

#Checking if reads are stranded
[I have no name!@bee76 Ariana_samples]$ IDX=ref/hisat2/Mcap_ref
[I have no name!@bee76 Ariana_samples]$ hisat2 -p 8 --dta --rna-strandness F -x "$IDX" -U test_SE.fastq \
  -S /dev/null 2> se_F.log
[I have no name!@bee76 Ariana_samples]$ hisat2 -p 8 --dta --rna-strandness R -x "$IDX" -U test_SE.fastq \
  -S /dev/null 2> se_R.log
[I have no name!@bee76 Ariana_samples]$ grep "overall alignment rate" se_F.log
grep "overall alignment rate" se_R.log
55.60% overall alignment rate
55.60% overall alignment rate
[I have no name!@bee76 Ariana_samples]$ hisat2 -p 8 --dta -x "$IDX" -U test_SE.fastq -S /dev/null 2> se_unstr.log
[I have no name!@bee76 Ariana_samples]$ grep "overall alignment rate" se_unstr.log
55.60% overall alignment rate
[I have no name!@bee76 Ariana_samples]$ hisat2 -p 8 --dta -x "$IDX" -U test_SE.fastq -S /dev/null 2> se_unstr.log
[I have no name!@bee76 Ariana_samples]$ grep "overall alignment rate" se_unstr.log
55.60% overall alignment rate

mkdir hisat2
cd hisat2

#Checking the design 
conda activate rseqc_env
infer_experiment.py -i bam/SRR14864064.bam -r "$BED"
#Reading reference gene model /lustre1/home/mass/amalia/montipora_mapping/ref/Mcap_transcripts_clean.bed ... Done
#Loading SAM/BAM file ...  Total 200000 usable reads were sampled


#This is PairEnd Data
#Fraction of reads failed to determine: 0.0006
#Fraction of reads explained by "1++,1--,2+-,2-+": 0.4953
#Fraction of reads explained by "1+-,1-+,2++,2--": 0.5042
#(rseqc_env) [I have no name!@bee76 Ariana_samples]$ infer_experiment.py -i bam/SRR22293445.bam -r "$BED"
#Reading reference gene model /lustre1/home/mass/amalia/montipora_mapping/ref/Mcap_transcripts_clean.bed ... Done
#Loading SAM/BAM file ...  Total 200000 usable reads were sampled


#This is SingleEnd Data
#Fraction of reads failed to determine: 0.0032
#Fraction of reads explained by "++,--": 0.9577
#Fraction of reads explained by "+-,-+": 0.0391
#(rseqc_env) [I have no name!@bee76 Ariana_samples]$ BED=/lustre1/home/mass/amalia/montipora_mapping/ref/Mcap_transcripts_clean.bed

# from /lustre1/home/mass/amalia/montipora_mapping/Ariana_samples
mkdir -p counts logs/featurecounts

featureCounts -T 8 \
  -a /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.gtf \
  -t exon -g gene_id \
  -p -B -C \
  -s 0 \
  -o counts/counts_pe.txt \
  bam/SRR14864*.bam \
  2> logs/featurecounts/featureCounts_pe.log

featureCounts -T 8 \
  -a /lustre1/home/mass/amalia/montipora_mapping/ref/Montipora_capitata_HIv3.gtf \
  -t exon -g gene_id \
  -s 0 \
  -o counts/counts_se.txt \
  bam/SRR22293*.bam \
  2> logs/featurecounts/featureCounts_se.log

ls -lh counts/counts_pe.txt counts/counts_se.txt
tail -n 30 logs/featurecounts/featureCounts_pe.log
tail -n 30 logs/featurecounts/featureCounts_se.log

#Cleaning the output
# 0) start fresh
rm -f pe.* se.* counts_all.*

# 1) Clean + label PE (skip comment lines, strip path/.bam, prefix PE_)
awk -F'\t' 'BEGIN{OFS="\t"}
  $1 ~ /^#/ {next}
  NR==1 {
    printf "Geneid";
    for (i=7; i<=NF; i++) { n=$i; gsub(/"/,"",n); sub(/^.*\//,"",n); sub(/\.bam$/,"",n); printf "\tPE_%s", n }
    print ""; next
  }
  { printf "%s", $1; for (i=7; i<=NF; i++) printf OFS "%s", $i; print "" }
' counts_pe.txt > pe.named

# 2) Clean + label SE (prefix SE_)
awk -F'\t' 'BEGIN{OFS="\t"}
  $1 ~ /^#/ {next}
  NR==1 {
    printf "Geneid";
    for (i=7; i<=NF; i++) { n=$i; gsub(/"/,"",n); sub(/^.*\//,"",n); sub(/\.bam$/,"",n); printf "\tSE_%s", n }
    print ""; next
  }
  { printf "%s", $1; for (i=7; i<=NF; i++) printf OFS "%s", $i; print "" }
' counts_se.txt > se.named

# 3) Sanity: make sure no commas snuck in
grep -n ',' pe.named se.named || true

# 4) Sort for a stable join
(head -n1 pe.named; tail -n +2 pe.named | sort -k1,1) > pe.sorted
(head -n1 se.named; tail -n +2 se.named | sort -k1,1) > se.sorted

# 5) Merge headers & bodies (fill missing with 0)
{ printf "%s", "$(head -n1 pe.sorted)"; printf "\t%s\n", "$(head -n1 se.sorted | cut -f2-)"; } > counts_all.header
join -t $'\t' -a1 -a2 -e 0 -o auto -1 1 -2 1 \
  <(tail -n +2 pe.sorted) <(tail -n +2 se.sorted) > counts_all.body

cat counts_all.header counts_all.body > counts_all.tsv

# 6) Verify the first row now
head -n 2 counts_all.tsv | cut -f-12

# make a clean version: strip commas, remove paths/.bam, fix leading comma on data lines
awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{
  gsub(/,/,"");                        # kill commas in header
  printf "Geneid";
  for(i=2;i<=NF;i++){ n=$i; sub(/^.*\//,"",n); sub(/\.bam$/,"",n); printf "\t%s", n }
  print ""; next
}
{
  gsub(/^,/,""); gsub(/,/,"");         # kill any stray commas in data lines
  print
}' counts_all.tsv > counts_all.clean.tsv

# quick peek
head -n 2 counts_all.clean.tsv | cut -f-12
