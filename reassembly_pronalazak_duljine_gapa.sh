#!/bin/bash

# Defaultne vrijednosti
BAM_PATH="/storage4/icacic/RAK_GUSTERACE_SASTAVLJANJE_GENOMA/herro_assembly_2/igv/simplex_ont.to_asm_herro.sorted.bam"
REF_PATH="/storage4/icacic/RAK_GUSTERACE_SASTAVLJANJE_GENOMA/herro_assembly_2/herro_hifiasm/HG008_T.asm.herro.fa"
REGION="h1tg000027l:9600000-10000000"
THREADS=64

usage() {
    echo "Upotreba: $0 <output_dir> [--bam PATH] [--ref PATH] [--region REGION] [--threads N]"
    exit 1
}

if [ "$#" -lt 1 ]; then
    usage
fi

OUT_DIR="$1"
shift
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --bam)
            BAM_PATH="$2"
            shift 2
            ;;
        --ref)
            REF_PATH="$2"
            shift 2
            ;;
        --region)
            REGION="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        *)
            echo "Nepoznata opcija: $1"
            usage
            ;;
    esac
done

OUT_DIR="$(realpath "$OUT_DIR")"
mkdir -p "$OUT_DIR"
LOGFILE="$OUT_DIR/gap_finding.log"
touch "$LOGFILE"

REGION_FILE="$OUT_DIR/region.txt"

echo "$REGION" > "$REGION_FILE"

cd "$OUT_DIR" || exit 1

run_cmd() {
    echo "[INFO] Izvršavam: $1"
    eval "$1" >> "$LOGFILE" 2>&1
}

cd "$OUT_DIR" || exit 1

# 1. Izdvajanje subregiona
run_cmd "samtools view -hb -F256 \"$BAM_PATH\" \"$REGION\" > subregions_herro.bam"

# 2. Pretvorba u FASTQ
run_cmd "samtools fastq subregions_herro.bam > subregions_herro.fastq"

# 3. De novo sastavljanje s Flye
run_cmd "flye --nano-corr subregions_herro.fastq --out-dir reassembly --threads $THREADS"

# 4. Poravnanje nove sekvence na referencu
run_cmd "minimap2 -ax asm5 \"$REF_PATH\" reassembly/assembly.fasta | samtools sort -o herro_aligned.bam"

# 5. Indexiranje rezultirajućeg BAM-a
run_cmd "samtools index herro_aligned.bam"

echo "[INFO] Gotovo! Detalji su zapisani u $LOGFILE"

