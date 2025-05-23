import time
import blosum as bl
from itertools import product

from upgma import (
    create_distance_matrix,
    upgma,
)
from parser import parse_args
from progressive_alignment import progressive_alignment
from read_write_file import read_seqs, fasta_to_clustal


def get_weight_matrix():
    matrix = bl.BLOSUM(62, default=0)
    acids = [acid for acid in matrix.keys() if acid != "*"]
    
    return {
        **{
            (acid1, acid2): float(matrix[acid1][acid2]) + 5
            for acid1 in acids
            for acid2 in acids
        },
        **{(acid, "-"): 0.0 for acid in acids},
        **{("-", acid): 0.0 for acid in acids},
        ("-", "-"): 0.0,
    }

def get_dna_matrix(match, mismatch, gap):
    matrix ={}
    for i, j in product(['A', 'T', 'G', 'C','-'], repeat = 2):
        if i == '-' and j == '-':
             matrix[(i,j)] = 0
             continue
        elif i == '-' or j == '-':
            matrix[(i,j)] = -gap
            matrix[(j,i)] = -gap
        elif i == j:
            matrix[(i,j)] = match
            matrix[(j, i)] = match
        else:
            matrix[(i,j)] = -mismatch
            matrix[(j,i)] = -mismatch
    return matrix

def get_ids_from_guide_tree(ids):
    stack = list(ids)
    res = []
    while stack:
        item = stack.pop()
        if isinstance(item, tuple):
            stack.extend(item) 
        elif item=='(' or item==')' or item==',':
            continue
        else:
            res.append(item)
    return list(reversed(res))


def main_():
    args = parse_args()
    molecule = args.molecule
    if molecule == 'DNA':
        weight_matrix = get_dna_matrix(args.match, args.mismatch, args.gap_open)
    else:
        weight_matrix = get_weight_matrix()


    sequences, names = read_seqs(args.filename, args.alignment_mode)
    distances, _ = create_distance_matrix(sequences=sequences)

    node = upgma(dist_matrix=distances)
    aligned_sequences = progressive_alignment(
        sequences=sequences, 
        guide_tree_root=node,
        weight_matrix=weight_matrix,
        gap_open = args.gap_open,
        gap_extend = args.gap_extension
    )

    ids = get_ids_from_guide_tree(node.id)

    fasta_to_clustal(ids,names,aligned_sequences)

if __name__ == "__main__":
    main_()
