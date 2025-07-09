#!/bin/bash

REF_FA="/storage4/icacic/RAK_GUSTERACE_SASTAVLJANJE_GENOMA/herro_assembly_2/herro_hifiasm/HG008_T.asm.herro.fa"
BAM_FILE="/storage4/icacic/RAK_GUSTERACE_SASTAVLJANJE_GENOMA/herro_assembly_2/igv/simplex_ont.to_asm_herro.sorted.bam"
OUTDIR="./out"

usage() {
    echo "Usage: $0 -c CONTIG -r REGIJE [-f REF_FA] [-b BAM_FILE] [-o OUTDIR]"
    echo "  -c CONTIG         Ime contiga"
    echo "  -r REGIJE         Lista regija za uklanjanje (oblik: \"start1-end1 start2-end2 ...\")"
    echo "  -f REF_FA         Putanja do referentnog FASTA (default: $REF_FA)"
    echo "  -b BAM_FILE       Putanja do BAM datoteke (default: $BAM_FILE)"
    echo "  -o OUTDIR         Izlazni direktorij (default: ./out)"
    exit 1
}

while getopts "c:r:f:b:o:" opt; do
  case $opt in
    c) CONTIG=$OPTARG ;;
    r) REGIONS=$OPTARG ;;
    f) REF_FA=$OPTARG ;;
    b) BAM_FILE=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    *) usage ;;
  esac
done

if [ -z "$CONTIG" ] || [ -z "$REGIONS" ]; then
    usage
fi

mkdir -p "$OUTDIR"
LOG_FILE="$OUTDIR/gap_closing.log"

echo "[INFO] Pokrećem gap closing za contig $CONTIG, regije: $REGIONS" > "$LOG_FILE"
echo "[INFO] Izlazni direktorij: $OUTDIR" >> "$LOG_FILE"

# 1. uklanjanje regije iz contig reference
echo "[STEP] Uklanjam regije iz contiga..." >> "$LOG_FILE"
python3 ./remove_gap.py "$REF_FA" "$CONTIG" $REGIONS >> "$LOG_FILE" 2>&1

MODIFIED_FA="${CONTIG}.modified.fasta"
if [ ! -f "$MODIFIED_FA" ]; then
    echo "[ERROR] Modificirani fasta ($MODIFIED_FA) nije pronađen!" >> "$LOG_FILE"
    exit 1
fi

mv "$MODIFIED_FA" "$OUTDIR/" >> "$LOG_FILE" 2>&1

# 2. Indeksiranje faste
echo "[STEP] Indeksiram modificirani fasta..." >> "$LOG_FILE"
samtools faidx "$OUTDIR/$MODIFIED_FA" >> "$LOG_FILE" 2>&1

# 3. Izdvajanje readova od contiga
echo "[STEP] Izdvajam readove iz BAM-a..." >> "$LOG_FILE"
samtools view -@64 -b "$BAM_FILE" "$CONTIG" > "$OUTDIR/contig_reads.bam" 2>> "$LOG_FILE"
samtools index "$OUTDIR/contig_reads.bam" >> "$LOG_FILE" 2>&1

# 4. Pretvorba BAM -> FASTQ
echo "[STEP] Pretvaram BAM u FASTQ..." >> "$LOG_FILE"
samtools fastq -@64 "$OUTDIR/contig_reads.bam" > "$OUTDIR/contig_reads.fastq" 2>> "$LOG_FILE"

# 5. Mapiranje readova na novu referencu
echo "[STEP] Mapiram readove na novi contig..." >> "$LOG_FILE"
minimap2 -ax map-ont -t 64 "$OUTDIR/$MODIFIED_FA" "$OUTDIR/contig_reads.fastq" \
| samtools view -hb -@64 -F256 - \
| samtools sort -@64 -o "$OUTDIR/realigned.sorted.bam" 2>> "$LOG_FILE"

# 6. Indeksiranje BAM-a
samtools index -@64 "$OUTDIR/realigned.sorted.bam" >> "$LOG_FILE" 2>&1

echo "[INFO] Gotovo! Rezultati su u: $OUTDIR" >> "$LOG_FILE"
