from itertools import combinations
import numpy as np

from pairwise_alignment import (
    needleman_wunsch,
    needleman_wunsch_affine,
)


def create_distance_matrix(sequences, match=1, mismatch=1, gap=1, gap_extend=None):
    """
    sequences (list): list of sequences
    match (int): score for matching characters
    mismatch (int): penalty for mismatching characters
    gap (int): penalty for gaps

    Returns pairwise distance matrix

    """
    n = len(sequences)
    dist_matrix = np.zeros((n, n))

    align_seqs = {}

    for i, j in combinations(range(n), 2):
        if gap_extend is None:
            align_seq1, align_seq2, _, score = needleman_wunsch(sequences[i][1], sequences[j][1], match, mismatch, gap)
        else:
            align_seq1, align_seq2, _, score = needleman_wunsch_affine(sequences[i][1], sequences[j][1], match, mismatch, gap, gap_extend)
        
        align_seqs[(i, j)] = (align_seq1, align_seq2)
        align_seqs[(j, i)] = (align_seq2, align_seq1)

        dist = 1 - (score / max(len(sequences[i]), len(sequences[j])))
        dist_matrix[i][j], dist_matrix[j][i] = dist, dist

    return dist_matrix, align_seqs


class UPGMA_Node:
    def __init__(self, id, children=None, height=0):
        self.id = id
        self.children = children if children else []
        self.height = height
        self.size = 1 if not children else sum(child.size for child in children)


def upgma(dist_matrix: np.ndarray):
    n = len(dist_matrix)
    clusters = [UPGMA_Node(i) for i in range(n)]
    distances = dist_matrix.copy()
    np.fill_diagonal(distances, np.inf)

    while len(clusters) > 1:
        i, j = np.unravel_index(np.argmin(distances), distances.shape)
        min_dist = distances[i][j] / 2

        new_cluster = UPGMA_Node(
            id=f"({clusters[i].id},{clusters[j].id})",
            children=[(clusters[i]), (clusters[j])],
            height=min_dist
        )
        new_distances = np.zeros((distances.shape[0] + 1, distances.shape[1] + 1))
        new_distances[:-1, :-1] = distances
        new_distances[-1, -1] = np.inf

        for k in range(len(distances)):
            if k != i and k != j:
                d = (distances[i][k] * clusters[i].size +
                     distances[j][k] * clusters[j].size) / (clusters[i].size + clusters[j].size)
                new_distances[k][-1] = d
                new_distances[-1][k] = d

        mask = np.ones(len(new_distances), dtype=bool)
        mask[[i, j]] = False
        distances = new_distances[mask][:, mask]

        i, j = sorted([i, j])
        del clusters[j]
        del clusters[i]
        clusters.append(new_cluster)

    return clusters[0]
