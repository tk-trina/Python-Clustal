"""
input: list of unligned sequences of DNA/proteins,
    or list of aligned sequences of DNA/proteins => dealign (user must confirm it)
output:
    list of unaligned sequences
"""

from Bio import SeqIO

def read_seqs(f, alignment_mode):
    if alignment_mode == 'unaligned':
        fasta_sequences = SeqIO.parse(open(f), 'fasta')
        seqs = []
        names = []
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            seqs.append(sequence)
            names.append(name)
        
    else:
        clustal_sequences = SeqIO.parse(open(f), 'clustal')
        seqs = []
        names = []
        for seq in clustal_sequences:
            name, sequence = seq.id, str(seq.seq)
            seqs.append((name, sequence))
            names.append(name)
        
    return seqs, names
