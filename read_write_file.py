
from Bio import SeqIO

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

def fasta_to_clustal(ids, names, sequences, output, line_length=60):
    '''
    ids (list): list of ids derived from the guide tree
    names (list): names of sequences from the original file
    sequences (list): list of aligned sequences
    output(str):  file name for writing the output

    Prints in stdout the alignment in clustal format 
    '''

    
    consensus = []
    for col in zip(*sequences):
        char = {c for c in col}
        if len(char) == 1:
            consensus.append('*')
        else:
            consensus.append(' ')
    
    max_name = max(names, key=lambda x: len(x))
    if output is None:
        print("CLUSTAL multiple sequence alignment\n")
        for i in range(0, len(sequences[0]), line_length):
            block_end = min(i+line_length, len(sequences[0]))
            
            for id, seq in zip(ids, sequences):
                print(f"{names[int(id)].ljust(len(max_name) + 5)}{seq[i:block_end]}")
            
            print(f"{' '.ljust(len(max_name) + 5)}{''.join(consensus[i:block_end])}\n")
    else:
        with open(output, "w", encoding="utf-8") as f:
            f.write("CLUSTAL multiple sequence alignment\n\n\n")
            for id, seq in zip(ids, sequences):
                f.write(f"{names[int(id)].ljust(len(max_name) + 5)}{seq}\n")
            f.write(f"{' '.ljust(len(max_name) + 5)}{''.join(consensus)}")