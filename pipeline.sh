#!/bin/bash
# =============================================================================
# Pipeline scRNA-seq — PBMC Humain (Homo sapiens, chr21)
# =============================================================================
# Prérequis : conda activate scrnaseq
# Auteur    : noebr
# Dataset   : PBMC 1k 10x Genomics v3 (sous-samplé ~300 cellules)
#             Source : https://zenodo.org/records/3457880
# Technologie : 10x Genomics Chromium v3
#               R1 = 28bp (16bp barcode + 12bp UMI)
#               R2 = 91bp (cDNA)
# =============================================================================

set -euo pipefail   # Arrêt immédiat en cas d'erreur

# =============================================================================
# ÉTAPE 0 — Arborescence du projet
# =============================================================================
mkdir -p ~/scrnaseq_exercise/{data/{raw,trimmed,genome},results/{qc,alignments,matrix},scripts,logs}

# =============================================================================
# ÉTAPE 1 — Téléchargement des données brutes
# =============================================================================
echo "=== ÉTAPE 1 : Téléchargement des données ==="

cd ~/scrnaseq_exercise/data/raw/

# Fichiers fastq.gz PBMC 1k 10x v3 (sous-samplé)
wget -q --show-progress \
    "https://zenodo.org/records/3457880/files/subset_pbmc_1k_v3_S1_L001_R1_001.fastq.gz" \
    -O pbmc_R1.fastq.gz

wget -q --show-progress \
    "https://zenodo.org/records/3457880/files/subset_pbmc_1k_v3_S1_L001_R2_001.fastq.gz" \
    -O pbmc_R2.fastq.gz

# Vérification des longueurs
echo "Longueur R1 (attendu : 28bp) :"
zcat pbmc_R1.fastq.gz | awk 'NR==2{print length($0), "bp"; exit}'
echo "Longueur R2 (attendu : 91bp) :"
zcat pbmc_R2.fastq.gz | awk 'NR==2{print length($0), "bp"; exit}'

echo "Done étape 1"

# =============================================================================
# ÉTAPE 2 — Contrôle qualité (FastQC + MultiQC)
# =============================================================================
echo "=== ÉTAPE 2 : Contrôle qualité ==="

# Décompression temporaire pour FastQC
zcat ~/scrnaseq_exercise/data/raw/pbmc_R1.fastq.gz \
    > ~/scrnaseq_exercise/data/raw/pbmc_R1_tmp.fastq
zcat ~/scrnaseq_exercise/data/raw/pbmc_R2.fastq.gz \
    > ~/scrnaseq_exercise/data/raw/pbmc_R2_tmp.fastq

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

# explorer.exe ~/scrnaseq_exercise/results/qc/   # WSL : ouvrir dans l'explorateur

echo "Done étape 2"

# =============================================================================
# ÉTAPE 3 — Trimming avec fastp
# NOTE : R1 (barcode+UMI) ne doit PAS être modifié
#        Seul R2 (cDNA) est trimmé des adaptateurs
# =============================================================================
echo "=== ÉTAPE 3 : Trimming ==="

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
    2> ~/scrnaseq_exercise/logs/pbmc_fastp.log

echo "Done étape 3"

# =============================================================================
# ÉTAPE 4 — Téléchargement du génome et de l'annotation
# NOTE : On utilise uniquement chr21 pour rester léger sur un portable
#        Pour un usage en production, télécharger le génome complet
# =============================================================================
echo "=== ÉTAPE 4 : Téléchargement du génome ==="

# Génome hg38 chr21 (UCSC)
wget -q --show-progress \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz" \
    -O ~/scrnaseq_exercise/data/genome/hg38_chr21.fa.gz
gunzip ~/scrnaseq_exercise/data/genome/hg38_chr21.fa.gz

# GTF Ensembl release 109 — extraction chr21
# UCSC utilise "chr21", Ensembl utilise "21" → on corrige avec sed
wget -q --show-progress \
    "https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.chr.gtf.gz" \
    -O ~/scrnaseq_exercise/data/genome/hg38.gtf.gz

zcat ~/scrnaseq_exercise/data/genome/hg38.gtf.gz | \
    awk '$1=="21"' | \
    sed 's/^21\t/chr21\t/' \
    > ~/scrnaseq_exercise/data/genome/hg38_chr21.gtf

rm ~/scrnaseq_exercise/data/genome/hg38.gtf.gz

echo "Done étape 4"

# =============================================================================
# ÉTAPE 5 — Génération de l'index STAR
# NOTE : genomeSAindexNbases 11 est adapté à un seul chromosome
#        Pour un génome humain complet → valeur par défaut (14)
# =============================================================================
echo "=== ÉTAPE 5 : Génération de l'index STAR ==="

mkdir -p ~/scrnaseq_exercise/data/genome/star_index

STAR --runMode genomeGenerate \
     --runThreadN 4 \
     --genomeDir ~/scrnaseq_exercise/data/genome/star_index \
     --genomeFastaFiles ~/scrnaseq_exercise/data/genome/hg38_chr21.fa \
     --sjdbGTFfile ~/scrnaseq_exercise/data/genome/hg38_chr21.gtf \
     --genomeSAindexNbases 11

echo "Done étape 5"

# =============================================================================
# ÉTAPE 6 — Téléchargement de la whitelist 10x v3
# =============================================================================
echo "=== ÉTAPE 6 : Whitelist 10x v3 ==="

wget -q --show-progress \
    "https://zenodo.org/records/3457880/files/3M-february-2018.txt.gz" \
    -O ~/scrnaseq_exercise/data/genome/whitelist_10x_v3.txt.gz

gunzip ~/scrnaseq_exercise/data/genome/whitelist_10x_v3.txt.gz

echo "Nombre de barcodes dans la whitelist :"
wc -l ~/scrnaseq_exercise/data/genome/whitelist_10x_v3.txt

echo "Done étape 6"

# =============================================================================
# ÉTAPE 7 — Alignement et comptage avec STARsolo
# =============================================================================
# Différences clés vs bulk RNAseq :
#   - readFilesIn : R2 (cDNA) en PREMIER, R1 (barcode+UMI) en SECOND
#   - soloType CB_UMI_Simple : mode 10x Genomics
#   - soloCBlen 16 : longueur du barcode cellulaire (10x v3)
#   - soloUMIlen 12 : longueur de l'UMI (10x v3)
#   - La matrice de comptage est générée directement (pas besoin de featureCounts)
# =============================================================================
echo "=== ÉTAPE 7 : Alignement STARsolo ==="

STAR \
    --runThreadN 4 \
    --genomeDir ~/scrnaseq_exercise/data/genome/star_index \
    --readFilesIn ~/scrnaseq_exercise/data/trimmed/pbmc_R2_trimmed.fastq.gz \
                  ~/scrnaseq_exercise/data/trimmed/pbmc_R1_trimmed.fastq.gz \
    --readFilesCommand zcat \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist ~/scrnaseq_exercise/data/genome/whitelist_10x_v3.txt \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM CB UB \
    --outFileNamePrefix ~/scrnaseq_exercise/results/alignments/pbmc_ \
    2> ~/scrnaseq_exercise/logs/pbmc_star.log

echo "Done étape 7"

# =============================================================================
# ÉTAPE 8 — Vérification de la matrice de sortie
# =============================================================================
echo "=== ÉTAPE 8 : Vérification de la matrice ==="

echo "--- Résumé STARsolo ---"
cat ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/Summary.csv

echo ""
echo "--- Matrice filtrée (vraies cellules) ---"
echo "Cellules détectées :"
wc -l ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/filtered/barcodes.tsv
echo "Gènes détectés :"
wc -l ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/filtered/features.tsv
echo "Dimensions de la matrice :"
head -n 3 ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/filtered/matrix.mtx

# Copie de la matrice finale dans results/matrix/
cp -r ~/scrnaseq_exercise/results/alignments/pbmc_Solo.out/Gene/filtered/ \
      ~/scrnaseq_exercise/results/matrix/

echo ""
echo "============================================="
echo " Pipeline terminé avec succès !"
echo " Matrice finale dans : results/matrix/filtered/"
echo "   - barcodes.tsv  : une cellule par ligne"
echo "   - features.tsv  : un gène par ligne"
echo "   - matrix.mtx    : counts UMI (format sparse)"
echo " Prochaine étape : Seurat (R) ou Scanpy (Python)"
echo "============================================="
