import time
import blosum as bl
from lolviz import *

from upgma import (
    create_distance_matrix,
    upgma,
)
from progressive_alignment import progressive_alignment


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


def main():
    weight_matrix = get_weight_matrix()

    with open("example.txt", "r") as fp:
        sequences = [s.strip() for s in fp.readlines()]

    distances, _ = create_distance_matrix(sequences=sequences)

    print(distances)

    node = upgma(dist_matrix=distances)

    aligned_sequences = progressive_alignment(
        sequences=sequences,
        guide_tree_root=node,
        weight_matrix=weight_matrix,
    )

    print(*aligned_sequences, sep="\n")

    g = objviz(node)
    g.view()

    time.sleep(60)


if __name__ == "__main__":
    main()
