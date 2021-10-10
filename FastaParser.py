from Bio import SeqIO
import os 

class FastaParser:

    def parse_fasta_file(self, file_name):
        dir_path = os.path.dirname(os.path.realpath(__file__))  + "\\data\\" + file_name
        my_seq = any
        for seq_record in SeqIO.parse(dir_path, "fasta"):
            my_seq = seq_record.seq
        return my_seq

