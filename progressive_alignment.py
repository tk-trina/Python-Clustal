def progressive_alignment(sequences, guide_tree_root):
    """
    sequences (list): list of sequences
    guide_tree_root(UPGMA_Node): root node of the guide tree

    Returns pairwise distance matrix
    """

    aligned_seq = {i: sequences[i] for i in range(len(sequences))}

    def align_nodes(node):
        if not node.children:
            return aligned_seq[node.id]

        left = align_nodes(node.children[0])
        right = align_nodes(node.children[1])

        aligned_profile1, aligned_profile2  = profile_alignment(left, right)
        return aligned_profile1, aligned_profile2

    return align_nodes(guide_tree_root)
