import re
import matplotlib.pyplot as plt
import argparse

def parse_fasta(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as file:
        sequence_id = None
        sequence = ''
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id is not None:
                    sequences[sequence_id] = sequence
                sequence_id = line[1:]  # Remove the '>' character
                sequence = ''
            else:
                sequence += line
        if sequence_id is not None:
            sequences[sequence_id] = sequence
    return sequences

def find_n_stretches(sequence):
    return [(m.start(), m.end()) for m in re.finditer(r'N+', sequence)]

def generate_histogram(n_stretches, output_file):
    lengths = [end - start for start, end in n_stretches]
    plt.hist(lengths, bins=range(1, max(lengths) + 2), edgecolor='black')
    plt.title('Histogram of N Stretches')
    plt.xlabel('Length of N Stretches')
    plt.ylabel('Frequency')
    plt.savefig(output_file)
    #plt.show()

def write_bed_file(n_stretches, sequence_id, bed_file):
    with open(bed_file, 'w') as file:
        for start, end in n_stretches:
            file.write(f'{sequence_id}\t{start}\t{end}\n')

def process_fasta(fasta_file, histogram_output, bed_output):
    sequences = parse_fasta(fasta_file)
    for sequence_id, sequence in sequences.items():
        n_stretches = find_n_stretches(sequence)
        generate_histogram(n_stretches, histogram_output)
        write_bed_file(n_stretches, sequence_id, bed_output)

def main():
    parser = argparse.ArgumentParser(description="Process a FASTA file to find 'N' stretches, generate a histogram, and output a BED file.")
    parser.add_argument('fasta_file', type=str, help="Input FASTA file")
    parser.add_argument('histogram_output', type=str, help="Output filename for the histogram (PNG format)")
    parser.add_argument('bed_output', type=str, help="Output filename for the BED file")

    args = parser.parse_args()

    process_fasta(args.fasta_file, args.histogram_output, args.bed_output)

if __name__ == "__main__":
    main()
