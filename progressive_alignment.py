import numpy as np

from upgma import UPGMA_Node
from pairwise_alignment import base_needleman_wunsch_affine


class Cluster:
    def __init__(self, seqs: list[str]):
        self.seqs = seqs


def cluster_alignment(
    first: Cluster,
    second: Cluster,
    weight_matrix: dict,
    gap_open: float,
    gap_extend: float,
) -> Cluster:
    def weight_gen(i: int, j: int) -> float:
        return np.array([
            weight_matrix[(seq1[i], seq2[j])]
            for seq1 in first.seqs
            for seq2 in second.seqs
        ]).mean()


    align1, align2, _, _ = base_needleman_wunsch_affine(
        seq1_len=len(first.seqs[0]),
        seq2_len=len(second.seqs[0]),
        weight_gen=weight_gen,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )

    return Cluster(
        seqs=[
            "".join([
                seq1[idx] if isinstance(idx, int) else '-'
                for idx in align1
            ])
            for seq1 in first.seqs
        ] + [
            "".join([
                seq2[idx] if isinstance(idx, int) else '-'
                for idx in align2
            ])
            for seq2 in second.seqs
        ]
    )




def progressive_alignment(
        sequences: list[str],
        guide_tree_root: UPGMA_Node,
        weight_matrix=None,
        gap_open: float = 1.0,
        gap_extend: float = 0.5,
):
    """
    Args:
        sequences (list): list of sequences
        guide_tree_root (UPGMA_Node): root node of the guide tree
        weight_matrix (dict): if weight matrix is None, aligns DNA sequences
    Returns:
        aligned sequences (list[str])
    """

    def inner(node: UPGMA_Node):
        if not node.children:
            return Cluster(seqs=[sequences[node.id]])

        if weight_matrix is None:
            return profile_alignment(
                Cluster1=inner(node=node.children[0]),
                Cluster2=inner(node=node.children[1]), gap_open=gap_open, gap_extend=gap_extend)

        return cluster_alignment(
            first=inner(node=node.children[0]),
            second=inner(node=node.children[1]),
            weight_matrix=weight_matrix,
            gap_open=gap_open,
            gap_extend=gap_extend,
        )

    return inner(node=guide_tree_root).seqs