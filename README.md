# Pipeline scRNA-seq — PBMC Humain (*Homo sapiens*, chr21)

Pipeline scRNA-seq pédagogique complète, des fichiers FASTQ bruts jusqu'à la matrice de comptage cellules × gènes. Conçue pour tourner en local (4 cœurs, 8 Go RAM) avec un jeu de données sous-samplé avant passage sur cluster.

## Contexte biologique

Cellules mononucléées du sang périphérique humain (**PBMC**) séquencées en **10x Genomics Chromium v3**.  
Dataset sous-samplé à ~300 cellules à des fins pédagogiques.

| Fichier | Contenu | Source |
|---------|---------|--------|
| `pbmc_R1.fastq.gz` | Barcode (16bp) + UMI (12bp) = 28bp | Zenodo 3457880 |
| `pbmc_R2.fastq.gz` | cDNA = 91bp | Zenodo 3457880 |

> **Note :** On utilise uniquement le chromosome 21 (hg38) pour réduire le temps de calcul. Pour une analyse en production, utiliser le génome complet.

## Pipeline

```
FASTQ bruts (R1 = barcode+UMI, R2 = cDNA)
    │
    ├─ FastQC + MultiQC     Contrôle qualité
    ├─ fastp                Trimming (R2 uniquement)
    ├─ STAR genomeGenerate  Génération de l'index
    ├─ STARsolo             Alignement + démultiplexage cellules
    │
    └─ filtered/
         ├── barcodes.tsv   1 cellule par ligne
         ├── features.tsv   1 gène par ligne
         └── matrix.mtx     Counts UMI (format sparse)
```

## Différences clés avec le bulk RNAseq

| | Bulk RNAseq | scRNA-seq (ce pipeline) |
|---|---|---|
| **Unité d'analyse** | Échantillon entier | Cellule individuelle |
| **Structure R1/R2** | cDNA + cDNA | Barcode+UMI + cDNA |
| **Outil de comptage** | featureCounts | STARsolo (intégré) |
| **Format de sortie** | TSV dense (gènes × échantillons) | MTX sparse (gènes × cellules) |
| **Gestion des duplicats** | Non | UMI dédupliqués automatiquement |
| **Filtrage** | N/A | raw → filtered (vraies cellules) |

## Arborescence

```
.
├── pipeline.sh             Script principal
├── environment.yml         Environnement conda
├── data/
│   ├── raw/                FASTQ bruts (.fastq.gz)
│   ├── trimmed/            FASTQ nettoyés
│   └── genome/             Génome hg38 chr21 + GTF + index STAR + whitelist
├── results/
│   ├── qc/                 Rapports FastQC / MultiQC / fastp
│   ├── alignments/         BAM + Solo.out (matrice brute)
│   └── matrix/             Matrice filtrée finale
└── logs/                   Logs de chaque étape
```

## Installation

### Prérequis

- Miniconda
- Linux ou WSL (Ubuntu)
- 4 cœurs, 8 Go RAM minimum

### Créer les environnements

```bash
# Environnement principal
conda env create -f environment.yml
conda activate scrnaseq

# MultiQC (environnement séparé — conflits de dépendances)
conda create -n multiqc -c conda-forge python=3.11 -y
conda activate multiqc
pip install multiqc
```

## Utilisation

```bash
conda activate scrnaseq
bash pipeline.sh
```

## Résultats obtenus

| Étape | Métrique | Valeur |
|-------|----------|--------|
| FastQC | Q30 R2 (cDNA) | 92.4% |
| fastp | Reads conservés | 99.8% |
| STARsolo | Reads avec barcode valide | 99.2% |
| STARsolo | Cellules détectées | **138** |
| STARsolo | Gènes détectés (chr21) | 316 |
| STARsolo | UMI médian par cellule | 280 |

## Prochaine étape

La matrice `results/matrix/filtered/` peut être chargée directement dans :

```r
# Seurat (R)
library(Seurat)
pbmc <- Read10X("results/matrix/filtered/")
seurat_obj <- CreateSeuratObject(counts = pbmc)
```

```python
# Scanpy (Python)
import scanpy as sc
adata = sc.read_10x_mtx("results/matrix/filtered/")
```

## Notes pédagogiques

- **R1 en second dans STARsolo** : contrairement au bulk, STARsolo attend `--readFilesIn R2 R1` (cDNA d'abord, barcode ensuite)
- **Whitelist 10x v3** : liste des 6,7 millions de barcodes valides 10x Chromium v3 — indispensable pour le démultiplexage
- **raw vs filtered** : `raw/` contient tous les barcodes de la whitelist ; `filtered/` ne conserve que les vraies cellules après l'algorithme CellRanger 2.2
- **chr21 uniquement** : le faible taux d'alignement (~9.4%) est normal — les reads provenant des 22 autres chromosomes ne s'alignent pas

## Références

- [STARsolo](https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md)
- [fastp](https://github.com/OpenGene/fastp)
- [10x Genomics PBMC 1k v3](https://www.10xgenomics.com/datasets/1-k-pbm-cs-from-a-healthy-donor-v-3-chemistry-3-standard-3-0-0)
- Dataset : [Zenodo 3457880](https://zenodo.org/records/3457880)
- Génome : [hg38 UCSC](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/)
- Annotation : [Ensembl release 109](https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/)
