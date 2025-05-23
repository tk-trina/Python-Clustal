
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from io import StringIO

def read_seqs(f, alignment_mode):
    '''
    f (str): an input filename
    alignment_mode {unaligned, aligned}: whether sequences need to be 
                                        preprocessed before alignment

    Returns sequences and their names from the original file
    '''
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
            name, sequence = seq.id, str(seq.seq).replace('-','')
            seqs.append(sequence)
            names.append(name)
        
    return seqs, names


def fasta_to_clustal(ids, names, sequences):
    '''
    ids (list): list of ids derived from the guide tree
    names (list): names of sequences from the original file
    sequences (list): list of aligned sequences

    Prints in stdout the alignment in clustal format 
    '''
    records = [SeqRecord(Seq(s), id=names[int(i)], description="") 
            for i, s in zip(ids, sequences)]
    alignment = MultipleSeqAlignment(records)
    output=StringIO()
    AlignIO.write(alignment, output, "clustal")
    clustal_str = output.getvalue()
    clustal_str = clustal_str[:7] + clustal_str[16:]
    print(clustal_str)