from Bio import SeqIO
import sys
import os

def parse_regions(region_list):
    """Parsiraj regije oblika 'start1-end1 start2-end2 ...' u listu tuple-ova (start, end)."""
    regions = []
    for region in region_list:
        try:
            start, end = map(int, region.strip().split("-"))
            if start >= end:
                raise ValueError
            regions.append((start, end))
        except ValueError:
            print(f"Neispravan format regije: {region}. Očekuje se oblik 'start-end'.")
            sys.exit(1)
    return regions

def remove_multiple_regions(fasta_path, contig_name, regions):
    output_filename = f"{contig_name}.modified.fasta"

    with open(fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == contig_name:
                sequence = record.seq

                # Sortiramo regije silazno kako bismo izbjegli pomake zbog brisanja
                regions_sorted = sorted(regions, reverse=True)

                for start, end in regions_sorted:
                    sequence = sequence[:start] + sequence[end:]

                with open(output_filename, "w") as out:
                    out.write(f">{record.id}\n{str(sequence)}\n")

                print(f"Spremljeno u: {output_filename}")
                return

    print(f"Contig '{contig_name}' nije pronađen u datoteci.")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Upotreba: python remove_regions_from_fasta.py <putanja.fasta> <contig> <start1-end1> [<start2-end2> ...]")
        sys.exit(1)

    fasta_file = sys.argv[1]
    contig = sys.argv[2]
    region_args = sys.argv[3:]

    regions = parse_regions(region_args)
    remove_multiple_regions(fasta_file, contig, regions)
