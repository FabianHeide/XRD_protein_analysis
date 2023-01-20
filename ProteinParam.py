#Various analysis classes and functions for the generation of additional properties from a given protein amino acid sequence

#imports (biopython)
import Bio
from Bio.SeqUtils import ProtParam

#Amino acid groupings
amino_acids = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','O','S','T','W','Y','V']
hpho_res = ['A','G','F','I','L','M','P','V','W']
hphi_res = ['R','K','H','D','E','S','T','N','Q','C','Y']
pos_res = ['R','K','H']
neg_res = ['D','E']
aro_res = ['F','Y','W']

class SequenceAnalysis():
    """This class analyzes basic content of a given protein sequence"""

    def __init__(self):
        pass

    def seq_length(self, sequence):
        sequence = sequence.upper()
        return len(sequence)

    def hphobic_res(self, sequence):
        pho_value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in hpho_res and i in amino_acids:
                pho_value += 1    
        return pho_value

    def hphilic_res(self, sequence):
        value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in hphi_res and i in amino_acids:
                value += 1     
        return value

    def positive_res(self, sequence):
        value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in pos_res and i in amino_acids:
                value += 1     
        return value

    def negative_res(self, sequence):
        value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in neg_res and i in amino_acids:
                value += 1     
        return value

    def aromatic_res(self, sequence):
        value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in aro_res and i in amino_acids:
                value += 1     
        return value

class ProteinProperties(SequenceAnalysis):
    """This class calculates various ratios and theoretical properties of the given protein sequence"""

    def __init__(self):
        super().__init__()

    #calculates the hydrophobic ratio of the entire sequence
    def hydrophobic_ratio(self, sequence):
        hphobic_res = self.hphobic_res(sequence)
        ratio = float(hphobic_res/len(sequence))
        return ratio

    #calculates the hydrophilic ratio of the entire sequence
    def hydrophilic_ratio(self, sequence):
        hphilic_res = self.hphilic_res(sequence)
        ratio = float(hphilic_res/len(sequence))
        return ratio

    #calculates the hydrophobic/hydrophilic ratio of the sequence
    def hydro_ratio(self, sequence):
        hphobic_res = self.hphobic_res(sequence)
        hphilic_res = self.hphilic_res(sequence)
        ratio = float(hphobic_res/hphilic_res)
        return ratio

    #calculates the aromatic ratio of the entire sequence
    def aro_ratio(self, sequence):
        aro_res = self.aromatic_res(sequence)
        ratio = float(aro_res/len(sequence))
        return ratio
    
    #calculates the isoelectric point of give sequence
    def isoelectric_point(self, sequence):
        x = ProtParam.ProteinAnalysis(sequence)
        pI = x.isoelectric_point()
        return pI

    #calculates the alipathic index of given sequence, higher index indicates a more stable protein
    def alipathic_index(self, sequence):
        sequence = sequence.upper()
        ali_index = float(100*sequence.count('A')/len(sequence) + 2.9*(100*sequence.count('V')/len(sequence)) + 3.9*(100*sequence.count('I')/len(sequence) + 100*sequence.count('L')/len(sequence)))
        return ali_index
