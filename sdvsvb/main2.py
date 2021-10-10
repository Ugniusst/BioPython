from Bio import SeqIO
from Bio.Seq import Seq
from tabulate import tabulate
import os 


dir_path = os.path.dirname(os.path.realpath(__file__)) 
my_seq = any
my_seq2 = any
for seq_record in SeqIO.parse(dir_path+"\\data\\mamalian1.fasta", "fasta"):
    my_seq = seq_record.seq
for seq_record in SeqIO.parse(dir_path+"\\data\\bacterial4.fasta", "fasta"):
    my_seq2 = seq_record.seq

start_Codon = Seq("ATG")
stop_Codons = [Seq("TAA"),Seq("TAG"),Seq("TGA")]
orfs = []
heading = ["X"]

def search_for_ORFS(frame, seq_input):
    seq_start = 0
    seq_stop = 0
    seq_started = False

    for x in range(frame,len(seq_input),3):
        codon = seq_input[x:x+3]
        if seq_started == False and codon == start_Codon:
            seq_start = x
            seq_started = True
        if seq_started == True and codon in stop_Codons:
            seq_stop = x+3
            if len(seq_input[seq_start:seq_stop]) >= 100:
                orfs.append(seq_input[seq_start:seq_stop])
            seq_started = False 

def concatenate_dna_sequences(dna_sequencies):
    concatenated = Seq("")
    for dna_sequence in dna_sequencies:
        concatenated += dna_sequence
    return concatenated

def count_codon_frequency(searchable_codon, dna_sequence):
    counter = 0
    for x in range(0, len(dna_sequence),3):
        codon = dna_sequence[x:x+3]
        if codon == searchable_codon:
            counter += 1
    frequency = (counter/(len(dna_sequence)/3))*1000
    return frequency

def count_all_codon_frequencies(dna_sequence):
    frequencies = [[[0 for j in range(4)] for k in range(4)] for l in range(4)]
    symbols = ["T", "C", "A", "G"]
    for s1 in range(4):
        for s2 in range(4):
            for s3 in range(4):
                codon = symbols[s1]+symbols[s2]+symbols[s3]
                frequencies[s1][s2][s3] = count_codon_frequency(Seq(codon), dna_sequence)
    return frequencies

def count_dicodon_frequency2(searchable_dicodon, dna_sequence):
    counter = 0
    for x in range(0, len(dna_sequence)-3,3):
        dicodon = dna_sequence[x:x+6]
        dicodon = get_protein(dicodon)
        if dicodon == searchable_dicodon:
            counter += 1
    frequency = (counter/((2 * len(dna_sequence))/6 -1))*1000
    return frequency

def count_dicodon_frequency(searchable_dicodon, dna_sequence):
    counter = 0
    dna_sequence = dna_sequence.translate()
    for x in range(0, len(dna_sequence)-1):
        dicodon = dna_sequence[x:x+2]
        if dicodon == searchable_dicodon:
            counter += 1
    frequency = (counter/((2 * len(dna_sequence))/2 -1))*1000
    return frequency

def get_protein(codon):
    codon = codon.translate()
    return codon

def get_all_proteins_list():
    proteins_list = []
    symbols = ["T", "C", "A", "G"]
    for s1 in range(4):
        for s2 in range(4):
            for s3 in range(4):
                codon = symbols[s1]+symbols[s2]+symbols[s3]
                if codon not in stop_Codons:
                    protein = get_protein(Seq(codon))
                    if protein not in proteins_list:
                        proteins_list.append(protein)
    return proteins_list

def count_all_dicodon_frequencies(dna_sequence, proteins_list):
    dicodon_frequencies = [[0 for j in range(len(proteins_list)+1)] for k in range(len(proteins_list))]
    for p1 in range(len(proteins_list)):
        heading.append(str(proteins_list[p1]))
        dicodon_frequencies[p1][0] = proteins_list[p1]
        for p2 in range(len(proteins_list)):
            dicodon = proteins_list[p1] + proteins_list[p2]
            dicodon_frequencies[p1][p2+1] = count_dicodon_frequency(dicodon, dna_sequence)
    return dicodon_frequencies

def count_all_dicodon_distances(matrix1, matrix2, proteins_list):
    dicodon_distances = [[0 for j in range(len(proteins_list)+1)] for k in range(len(proteins_list))]
    for p1 in range(len(proteins_list)):
        heading.append(str(proteins_list[p1]))
        dicodon_distances[p1][0] = proteins_list[p1]
        for p2 in range(len(proteins_list)):
            fraction1 = matrix1[p1][p2+1]/(matrix1[p1][p2+1]+matrix2[p1][p2+1])
            fraction2 = matrix2[p1][p2+1]/(matrix1[p1][p2+1]+matrix2[p1][p2+1])
            dicodon_distances[p1][p2+1] = round(abs(fraction1-fraction2),1)
    return dicodon_distances

prot = get_all_proteins_list()
dna = []

for frame in range(0,3): 
    search_for_ORFS(frame, my_seq)
    search_for_ORFS(frame, my_seq.reverse_complement())
    #dna = searchFurthestStart(orfs)

#print(dna)
print("Aaaa")
#cccc = concatenate_dna_sequences(orfs)
print(len(orfs))
my_seq = concatenate_dna_sequences(orfs)
#print(my_seq)
#frequency = count_codon_frequency(Seq("ATG"), my_seq)
matrix1 = count_all_dicodon_frequencies(my_seq, prot)
#print(tabulate(matrix1))

orfs = []
dna = []
print("1 failas")
for frame in range(0,3):
    print(frame)
    search_for_ORFS(frame, my_seq2)
    search_for_ORFS(frame, my_seq2.reverse_complement())
    #dna = searchFurthestStart(orfs)
print(len(orfs))
print(orfs)
my_seq2 = concatenate_dna_sequences(orfs)
#frequency = count_codon_frequency(Seq("ATG"), my_seq)
matrix2 = count_all_dicodon_frequencies(my_seq2, prot)
print("2 failas")

print(tabulate(count_all_dicodon_distances(matrix1, matrix2, prot),heading))
  