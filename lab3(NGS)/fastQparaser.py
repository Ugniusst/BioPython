from io import RawIOBase
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import matplotlib.pyplot as plt

RANGES = {
    'Sanger': (33, 73),
    'Illumina-1.8': (33, 74),
    'Solexa': (59, 104),
    'Illumina-1.3': (64, 104),
    'Illumina-1.5': (66, 105)
}


def get_quality_range(qual_str):
    min_base_qual = min(qual_str)
    max_base_qual = max(qual_str)

    return (min_base_qual, max_base_qual)


def get_encodings_in_range(rmin, rmax, ranges=RANGES):
    valid_encodings = []
    for encoding, (emin, emax) in ranges.items():
        if ord(rmin) >= emin and ord(rmax) <= emax:
            valid_encodings.append(encoding)
    return valid_encodings

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def calculate_gc_frequency(sequence):
    gc_counter = sequence.count('G') + sequence.count('C')
    return gc_counter/len(sequence)

valid_encoding = RANGES
count = 0
total_len = 0
gc_frequencies = []
reads = []
big_gc_frequencies = []
big_reads = []

reads_freqs = [0 for i in range(100)]

with open("reads_for_analysis.fastq") as in_handle:
    for title, seq, qual in FastqGeneralIterator(in_handle):
        count += 1
        reads.append(count)
        gc_freq = calculate_gc_frequency(seq)
        gc_frequencies.append(round(gc_freq*100))
        reads_freqs[round(gc_freq*100)] = reads_freqs[round(gc_freq*100)] + 1
        if(round(gc_freq*100) == 70):
            big_gc_frequencies.append(seq)
            big_reads.append(count)
            print(str(count))
            print(title)
            print(seq+"\n\n")
        total_len += len(seq)
        qual_range = get_quality_range(qual)
        #print(qual_range)
        valid_encoding = intersection(valid_encoding,get_encodings_in_range(qual_range[0], qual_range[1]))
            #print(get_encodings_in_range(qual_range[0], qual_range[1]))

print(valid_encoding)

percentages = [i for i in range(100)]
plt.plot(percentages, reads_freqs, label=percentages)
for i,j in zip(percentages,reads_freqs):
    plt.annotate(str(percentages[i]),xy=(i,j))
plt.show()

print("%i records with total sequence length %i" % (count, total_len))

