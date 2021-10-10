from Bio.Seq import Seq

class ORFsFinder:
    start_Codon = Seq("ATG")
    stop_Codons = [Seq("TAA"),Seq("TAG"),Seq("TGA")]
    orfs = []

    def search_for_ORFS(self, frame, seq_input):
        seq_start = 0
        seq_stop = 0
        seq_started = False

        for x in range(frame,len(seq_input),3):
            codon = seq_input[x:x+3]
            if seq_started == False and codon == self.start_Codon:
                seq_start = x
                seq_started = True
            if seq_started == True and codon in self.stop_Codons:
                seq_stop = x+3
                if len(seq_input[seq_start:seq_stop]) >= 100:
                    self.orfs.append(seq_input[seq_start:seq_stop])
                seq_started = False
            

    def search_all_frames(self, sequence):
        for frame in range(0,3): 
            self.search_for_ORFS(frame, sequence)
            self.search_for_ORFS(frame, sequence.reverse_complement())


    def get_all_ORFS(self, sequence):
        self.orfs = []
        self.search_all_frames(sequence)
        return self.orfs