# Trim the adapters from small RNA datasets.

# Convert BAM to fastq files.
for BAM in /mnt/data/ccbc_environment/project/hCRISPR/externalData/Celllines_Trimmed/*bam;
  do echo "bedtools bamtofastq -i ${BAM} -fq ${BAM/bam/fq}";
done

# Trim using fastp (v0.20.0)
for fq in *fq;
  do echo "/home/jvanriet/.local/bin/cutadapt -a file:/mnt/data/ccbc_environment/project/hCRISPR/externalData/GSE80400/adapters.fa --minimum-length 15 -j 20 -o ${fq/fq/trimmed.fq} ${fq}";
done  
# 

# ENCODE 

# PCa
