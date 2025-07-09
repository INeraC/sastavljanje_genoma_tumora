from Bio import SeqIO
import sys

def find_gaps(fasta_file):
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        scaffold_name = record.id
        gaps = []
        
        in_gap = False
        gap_start = None
        
        for i, base in enumerate(seq):
            if base == "N":
                if not in_gap:
                    gap_start = i  # FASTA indeksi su 1-bazirani
                    in_gap = True
            else:
                if in_gap:
                    gap_length = i - gap_start
                    gaps.append((gap_start, i, gap_length))  # Dodaj trojku (Pocetak, Kraj, Duljina)
                    in_gap = False

        # Ako sekvenca zavrsava gapom
        if in_gap:
            gap_length = len(seq) - gap_start + 1
            gaps.append((gap_start, len(seq), gap_length))

        # Ispisuje samo ako ima gapova
        if gaps:
            print(f"{scaffold_name}")
            print(len(gaps))
            print(gaps)
            print("\n\n")  # Tri prazna retka za razdvajanje scaffolda

# Pokreni funkciju nad FASTA datotekom
#fasta_file = "supercontigs.fa"
#find_gaps(fasta_file)
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python gaps_positions.py <fasta_file>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    find_gaps(fasta_file)
