#!/bin/bash
# =============================================================================
# Pipeline scRNA-seq — PBMC Humain (Homo sapiens, chr21)
# =============================================================================
# Prérequis : conda activate scrnaseq
# Auteur    : noebr
# Technologie : 10x Genomics Chromium v3
#               R1 = 28bp (16bp barcode + 12bp UMI)
#               R2 = 91bp (cDNA)
# =============================================================================
set -euo pipefail

# Création des dossiers
mkdir -p ~/scrnaseq_exercise/{data/{raw,trimmed,genome},results/{qc,alignments,matrix},scripts,logs}
#Recup des données 10x genomics : pbmc est un jeu classique
cd ~/scrnaseq_exercise/data/raw/
wget -q --show-progress \
    "https://zenodo.org/records/3457880/files/subset_pbmc_1k_v3_S1_L001_R1_001.fastq.gz" \
    -O pbmc_R1.fastq.gz
wget -q --show-progress \
    "https://zenodo.org/records/3457880/files/subset_pbmc_1k_v3_S1_L001_R2_001.fastq.gz" \
    -O pbmc_R2.fastq.gz

# Contrôle qualité (FastQC + MultiQC)
#décompresse les fastq en temporaire, fais le qc puis les effaces pour en garder que les ficheirs comrpéssés
zcat ~/scrnaseq_exercise/data/raw/pbmc_R1.fastq.gz > ~/scrnaseq_exercise/data/raw/pbmc_R1_tmp.fastq
zcat ~/scrnaseq_exercise/data/raw/pbmc_R2.fastq.gz > ~/scrnaseq_exercise/data/raw/pbmc_R2_tmp.fastq
fastqc \
    ~/scrnaseq_exercise/data/raw/pbmc_R1_tmp.fastq \
    ~/scrnaseq_exercise/data/raw/pbmc_R2_tmp.fastq \
    --outdir ~/scrnaseq_exercise/results/qc/ \
    --threads 4
rm ~/scrnaseq_exercise/data/raw/pbmc_R1_tmp.fastq \
   ~/scrnaseq_exercise/data/raw/pbmc_R2_tmp.fastq

# MultiQC (environnement séparé)
conda run -n multiqc multiqc \
    ~/scrnaseq_exercise/results/qc/ \
    --outdir ~/scrnaseq_exercise/results/qc/ \
    --filename multiqc_report

# explorer.exe . ~/scrnaseq_exercise/results/qc/   

# Trimming avec fastp
# NOTE : R1 (barcode+UMI) ne doit PAS être modifié
#        Seul R2 (cDNA) est trimmé des adaptateurs

fastp \
    --in1  ~/scrnaseq_exercise/data/raw/pbmc_R1.fastq.gz \
    --in2  ~/scrnaseq_exercise/data/raw/pbmc_R2.fastq.gz \
    --out1 ~/scrnaseq_exercise/data/trimmed/pbmc_R1_trimmed.fastq.gz \
    --out2 ~/scrnaseq_exercise/data/trimmed/pbmc_R2_trimmed.fastq.gz \
    --json ~/scrnaseq_exercise/results/qc/pbmc_fastp.json \
    --html ~/scrnaseq_exercise/results/qc/pbmc_fastp.html \
    --thread 4 \
    --detect_adapter_for_pe \
    --disable_quality_filtering \
    --length_required 20 \                                        
    #length_required 15 pour du singel nuclei
    2> ~/scrnaseq_exercise/logs/pbmc_fastp.log

# Téléchargement du génome et de l'annotation
#  On utilise uniquement chr21 pour rester léger sur un portable

# Génome hg38 chr21 (UCSC)
wget -q --show-progress \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz" \
    -O ~/scrnaseq_exercise/data/genome/hg38_chr21.fa.gz

gunzip ~/scrnaseq_exercise/data/genome/hg38_chr21.fa.gz

# UCSC utilise "chr21", Ensembl utilise "21" 
wget -q --show-progress \
    "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz" \
    -O ~/scrnaseq_exercise/data/genome/hg38.gtf.gz

zcat ~/scrnaseq_exercise/data/genome/hg38.gtf.gz | \
    awk '$1=="21"' | \
    sed 's/^21\t/chr21\t/' \
    > ~/scrnaseq_exercise/data/genome/hg38_chr21.gtf

rm ~/scrnaseq_exercise/data/genome/hg38.gtf.gz



# =============================================================================
#  Génération de l'index STAR
# NOTE : genomeSAindexNbases 11 est adapté à un seul chromosome
#        Pour un génome humain complet → valeur par défaut (14)



mkdir -p ~/scrnaseq_exercise/data/genome/star_index

STAR --runMode genomeGenerate \
     --runThreadN 4 \
     --genomeDir ~/scrnaseq_exercise/data/genome/star_index \
     --genomeFastaFiles ~/scrnaseq_exercise/data/genome/hg38_chr21.fa \
     --sjdbGTFfile ~/scrnaseq_exercise/data/genome/hg38_chr21.gtf \
     --genomeSAindexNbases 11




#Whitelist 10x

wget -q --show-progress \
    "https://zenodo.org/records/3457880/files/3M-february-2018.txt.gz" \
    -O ~/scrnaseq_exercise/data/genome/whitelist_10x_v3.txt.gz

gunzip ~/scrnaseq_exercise/data/genome/whitelist_10x_v3.txt.gz


# Alignement et comptage avec STARsolo

# Différences vs bulk RNAseq :
#   - readFilesIn : R2 (cDNA) en PREMIER, R1 (barcode+UMI) en SECOND
#   - soloType CB_UMI_Simple : mode 10x Genomics
#   - soloCBlen 16 : longueur du barcode cellulaire (10x v3)
#   - soloUMIlen 12 : longueur de l'UMI (10x v3)


STAR \
    --runThreadN 4 \
    --genomeDir ~/scrnaseq_exercise/data/genome/star_index \
    --readFilesIn ~/scrnaseq_exercise/data/trimmed/pbmc_R2_trimmed.fastq.gz \                
    #R2 = sequence à aligner
                  ~/scrnaseq_exercise/data/trimmed/pbmc_R1_trimmed.fastq.gz \               
                #R1 = barcode + UMI
    --readFilesCommand zcat \                                                                
    #décompresse els fichiers
    --soloType CB_UMI_Simple \                                                               
    # car singel cell
    --soloCBwhitelist ~/scrnaseq_exercise/data/genome/whitelist_10x_v3.txt \
    --soloCBstart 1 --soloCBlen 16 \                                                          
    #positon du barcode : commence en 1 et et comprend 16 base
    --soloUMIstart 17 --soloUMIlen 12 \                                                        
    #posiiton UMI : commence en 17, comprends 12 abses 
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM CB UB \
    --outFileNamePrefix ~/scrnaseq_exercise/results/alignments/pbmc_ \
    2> ~/scrnaseq_exercise/logs/pbmc_star.log



# Pour singel nuclei :

#STAR \
 #   --runThreadN 4 \
  #  --genomeDir ~/scrnaseq_exercise/data/genome/star_index \
   # --readFilesIn ~/scrnaseq_exercise/data/trimmed/pbmc_R2_trimmed.fastq.gz \
    #              ~/scrnaseq_exercise/data/trimmed/pbmc_R1_trimmed.fastq.gz \
    #--readFilesCommand zcat \
    #--soloType CB_UMI_Simple \
    #--soloCBwhitelist ~/scrnaseq_exercise/data/genome/whitelist_10x_v3.txt \
    #--soloCBstart 1 --soloCBlen 16 \
    #--soloUMIstart 17 --soloUMIlen 12 \
    #--soloFeatures GeneFull \
    #--outSAMtype BAM SortedByCoordinate \
    #--outSAMattributes NH HI AS NM CB UB \
    #--outFileNamePrefix ~/scrnaseq_exercise/results/alignments/snrnaseq_ \
    #2> ~/scrnaseq_exercise/logs/snrnaseq_star.log

cat ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/Summary.csv


#echo "--- Matrice filtrée (vraies cellules) ---"
#echo "Cellules détectées :"
#wc -l ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/filtered/barcodes.tsv
#echo "Gènes détectés :"
#wc -l ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/filtered/features.tsv
#echo "Dimensions de la matrice :"
#head -n 3 ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/filtered/matrix.mtx

# Copie de la matrice finale dans results/matrix/
cp -r ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/filtered/ \
      ~/scrnaseq_exercise/results/matrix/                                                                            
      #GeneFull pour singlenuclei


gzip ~/scrnaseq_exercise/results/matrix/filtered/barcodes.tsv
gzip ~/scrnaseq_exercise/results/matrix/filtered/features.tsv
gzip ~/scrnaseq_exercise/results/matrix/filtered/matrix.mtx
