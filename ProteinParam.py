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

    def hydrophobic_res(self, sequence):
        pho_value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in hpho_res:
                pho_value += 1    
        return pho_value

    def hydrophilic_res(self, sequence):
        value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in hphi_res:
                value += 1     
        return value

    def positive_res(self, sequence):
        value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in pos_res:
                value += 1     
        return value

    def negative_res(self, sequence):
        value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in neg_res:
                value += 1     
        return value

    def aromatic_res(self, sequence):
        value = 0
        sequence = sequence.upper()
        for i in sequence:
            if i in aro_res:
                value += 1     
        return value

class ProteinProperties(SequenceAnalysis):
    """This class calculates various ratios and theoretical properties of the given protein sequence"""

    def __init__(self):
        super().__init__()

    #calculates the hydrophobic ratio of the entire sequence
    def hydrophobic_ratio(self, sequence):
        hphobic_res = self.hydrophobic_res(sequence)
        ratio = float(hphobic_res/len(sequence))
        return ratio

    #calculates the hydrophilic ratio of the entire sequence
    def hydrophilic_ratio(self, sequence):
        hphilic_res = self.hydrophilic_res(sequence)
        ratio = float(hphilic_res/len(sequence))
        return ratio

    #calculates the hydrophobic/hydrophilic ratio of the sequence
    def hydro_ratio(self, sequence):
        hphobic_res = self.hydrophobic_res(sequence)
        hphilic_res = self.hydrophilic_res(sequence)
        ratio = float(hphobic_res/hphilic_res)
        return ratio

    #calculates the aromatic ratio of the entire sequence
    def aromatic_ratio(self, sequence):
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

    #estimates overall charge at specified pH (theoretical and highly dependent on physicochemical environment)
    def charge_ph(self, sequence, ph=7.0):
        carboxy_values = {'D':3.65,'E':4.25}
        amino_values = {'K':10.53,'R':12.48,'H':6.0}
        carboxy_list = [3.3]
        amino_list = [7.7]
        for i in sequence.upper():
            if i in carboxy_values:
                carboxy_list.append(carboxy_values[i])
            elif i in amino_values:
                amino_list.append(amino_values[i])
        below = len([1 for i in carboxy_list if i >= ph])
        above = len([1 for i in amino_list if i < ph])
        overall_charge = above-below
        return overall_charge

    #calculates the MW in Daltons for sequence
    def molecular_weight(self, sequence):
        amino_acid_weights = {'A':89.09, 'C':121.16, 'D':133.1, 'E':147.13, 'F':165.19, 'G':75.07, 'H':155.16, 
                                'I':131.18, 'K':146.19, 'L':131.18, 'M':149.21, 'N':132.12, 'P':115.13, 'Q':146.15,
                                'R':174.2, 'S':105.09, 'T':119.12, 'V':117.15, 'W':204.23, 'Y':181.19}
        mw = 16.02 - (18.02*(len(sequence)-1))
        for i in sequence.upper():
            if i in amino_acid_weights:
                mw += amino_acid_weights[i]
        return mw

class ProteinMotifs(SequenceAnalysis):
    """
    This class contains methods for finding certain structure, glycosylation or ligand binding motifs in the protein 
    sequence. Predictions are purely theoretical and it should be noted that sequence motifs can occur randomly. Also, 
    this list of motifs is by no means exhaustive.
    """

    def __init__(self):
        super().__init__()

    #checks sequence for absence (0) or presence (1) of Cardin-Weintraub motif (heparin binding motif)
    def cardin_weintraub(self, sequence):
        basic_res = ['R','K','H']
        seq_mod = []
        for i in sequence.upper():
            if i in basic_res:
                seq_mod.append('B')
            else:
                seq_mod.append('X')
        if 'XBBXBX' in ''.join(str(i) for i in seq_mod) or 'XBXBBX' in ''.join(str(i) for i in seq_mod):
            motif = 1
        elif 'XBBBXXBX' in ''.join(str(i) for i in seq_mod) or 'XBXXBBBX' in ''.join(str(i) for i in seq_mod):
            motif = 1
        else:
            motif = 0
        return motif

    #checks sequence for the any N-glycosylation patterns
    def n_glycosylation(self, sequence):
        x_residues = ['A','R','D','C','Q','E','G','H','I','L','K','M','F','O','W','Y','V']
        motif_seq = ['NXS','NXT','NSS','NTS','NST','NTT']
        seq_mod = []
        for i in sequence.upper():
            if i in x_residues:
                seq_mod.append('X')
            else:
                seq_mod.append(i)
        if i in ''.join(str(i) for i in motif_seq):
            motif = 1
        else:
            motif = 0
        return motif
        