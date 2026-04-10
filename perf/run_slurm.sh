#!/bin/bash
#SBATCH --job-name=agat_perf
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=nextflow_%j.out
#SBATCH --error=nextflow_%j.err

# Script pour lancer le pipeline Nextflow sur SLURM
# Ce script soumet le pipeline principal, qui créera ensuite des jobs SLURM individuels

# Charger les modules nécessaires (adapter selon votre cluster)
# module load nextflow
# module load singularity

# Variables
GFF_FILE="${1:-Homo_sapiens.GRCh38.114.chr.4171206.gff3}"
SIZES="${2:-100000,500000,1000000,2000000,4171206}"
CPUS="${3:-0,1,2,4,8}"

echo "==================================="
echo "AGAT Performance Benchmark on SLURM"
echo "==================================="
echo "GFF file: $GFF_FILE"
echo "Sizes: $SIZES"
echo "CPUs: $CPUS"
echo "==================================="

# Lancer le pipeline avec le profil SLURM
nextflow run performance.nf \
    -profile slurm_singularity \
    --gff "$GFF_FILE" \
    --sizes "$SIZES" \
    --cpus "$CPUS" \
    -resume

echo "Pipeline terminé avec le code: $?"
