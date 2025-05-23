
import numpy as np
from progressive_alignment import Cluster
from itertools import product

def add_value_to_profile(profile, value_column):
    """
       Args:
           profile (list(list)): new sequence profile
           value_column (np.ndarray): column of the previous sequence profile

      Adds column from the  previous sequence profile to a new one
    """
    for p, value in zip(profile, map(str,value_column)):
        p.append(value)

def add_gap_to_profile(profile):
    """
          Args:
              profile (list(list)): new sequence profile

         Adds gap column  to a new profile
       """
    for p in profile:
        p.append('-')

def count_columns_score(column1, column2, match, mismatch, gap):
    """
        Args:
            column1 (np.ndarray): column of the first sequence profile
            column2 (np.ndarray): column of the second sequence profile
            match (int): score for matching characters
            mismatch (int): penalty for mismatching characters
            gap (int): penalty for gaps
        Returns:
            score(float):  score for the column from the first profile and the column from the second profile
    """
    count = 0
    score = 0
    for i, j in itertools.product(column1, column2):
        if i == '-' and j =='-':
            continue
        elif i == j:
            score += match
        elif i != j and (i == '-' or j =='-' ):
            score -= gap
        else:
            score -= mismatch
        count += 1
    if score == 0:
        return 0
    return score /count


def profile_alignment_simple(Cluster1, Cluster2, match=1, mismatch=1, gap=2):
    """
        Args:
            Cluster1 (Cluster): the first sequence profile
            Cluster2 (Cluster): the second sequence profile
            match (int): score for matching characters
            mismatch (int): penalty for mismatching characters
            gap (int): penalty for gaps

        Returns a new aligned profile (Cluster)

    """
    seq_list1 = Cluster1.seqs
    seq_list2 = Cluster2.seqs
    seq_list1 = np.array([list(seq) for seq in seq_list1])
    seq_list2 = np.array([list(seq) for seq in seq_list2])
    n1, l1 = seq_list1.shape
    n2, l2 = seq_list2.shape
    dp = np.zeros((l1 + 1, l2 + 1))
    dp[:, 0] = np.linspace(0, -l1 * gap, l1 + 1)
    dp[0, :] = np.linspace(0, -l2 * gap, l2 + 1)
    traceback = np.zeros((l1 + 1, l2 + 1))

    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            diag_score = dp[i - 1][j - 1] + count_columns_score(seq_list1[:, i - 1], seq_list2[:, j - 1], match,
                                                                mismatch, gap)
            up_score = dp[i - 1][j] - gap
            left_score = dp[i][j - 1] - gap

            scores = (diag_score, up_score, left_score)
            max_direction = np.argmax(scores)

            dp[i][j] = scores[max_direction]
            traceback[i][j] = max_direction

    align_prof1 = [[] for i in range(n1)]
    align_prof2 = [[] for i in range(n2)]
    i, j = l1, l2

    while i > 0 and j > 0:
        if traceback[i][j] == 0:
            add_value_to_profile(align_prof1, seq_list1[:, i - 1])
            add_value_to_profile(align_prof2, seq_list2[:, j - 1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 1:
            add_value_to_profile(align_prof1, seq_list1[:, i - 1])
            add_gap_to_profile(align_prof2)
            i -= 1
        else:
            add_gap_to_profile(align_prof1)
            add_value_to_profile(align_prof2, seq_list2[:, j - 1])
            j -= 1
    while i > 0:
        add_value_to_profile(align_prof1, seq_list1[:, i - 1])
        add_gap_to_profile(align_prof2)
        i -= 1
    while j > 0:
        add_gap_to_profile(align_prof1)
        add_value_to_profile(align_prof2, seq_list2[:, j - 1])
        j -= 1
    profile = list(map(lambda x: ''.join(reversed(x)), align_prof1)) + list(
        map(lambda x: ''.join(reversed(x)), align_prof2))
    return Cluster(profile)

def count_columns_score_affine(column1, column2, match, mismatch, gap_open, gap_extend):
    """
        Args:
            column1 (np.ndarray): column of the first sequence profile
            column2 (np.ndarray): column of the second sequence profile
            match (int): score for matching characters
            mismatch (int): penalty for mismatching characters
            gap_open (int): penalty for opening gaps
            gap_extend (float): penalty for extending gaps
        Returns:
            score(float):  score for the column from the first profile and the column from the second profile
    """
    count = 0
    score = 0
    prev_gap_i = False
    prev_gap_j = False
    for i, j in itertools.product(column1, column2):
        if i == '-' and j =='-':
           prev_gap_i = True
           prev_gap_j = True
           continue
        elif i == j:
            score += match
            prev_gap_i = False
            prev_gap_j = False
        elif i != j and (i == '-' or j =='-' ):
            if (i == '-' and prev_gap_j) or (j == '-' and prev_gap_i):
                score -= gap_extend
            else:
                score -= gap_open
            prev_gap_i = (i == '-')
            prev_gap_j = (j == '-')
        else:
            score -= mismatch
            prev_gap_i = False
            prev_gap_j = False

        count += 1
    if score == 0:
        return 0
    return score /count


def profile_alignment_affine(Cluster1, Cluster2, match=1, mismatch=1, gap_open=1, gap_extend=0.5):
    """
        Args:
            Cluster1 (Cluster): the first sequence profile
            Cluster2 (Cluster): the second sequence profile
            match (int): score for matching characters
            mismatch (int): penalty for mismatching characters
            gap_open (int): penalty for opening gaps
            gap_extend (float): penalty for extending gaps

        Returns a new aligned profile (Cluster)

    """
    seq_list1 = Cluster1.seqs
    seq_list2 = Cluster2.seqs
    seq_list1 = np.array([list(seq) for seq in seq_list1])
    seq_list2 = np.array([list(seq) for seq in seq_list2])
    n1, l1 = seq_list1.shape
    n2, l2 = seq_list2.shape

    dp = np.zeros((l1 + 1, l2 + 1))
    dp[1:, 0] = np.linspace(-gap_open, -(l1 - 1) * gap_extend, l1)
    dp[0, 1:] = np.linspace(-gap_open, -(l2 - 1) * gap_extend, l2)

    I1, I2 = np.zeros((l1 + 1, l2 + 1)), np.zeros((l1 + 1, l2 + 1))
    I1[1:, 0] = np.linspace(-gap_open, -(l1 - 1) * gap_extend, l1)
    I2[0, 1:] = np.linspace(-gap_open, -(l2 - 1) * gap_extend, l2)
    traceback = np.zeros((l1 + 1, l2 + 1))

    for i in range(1, l1 + 1):
        for j in range(1, l2 + 1):
            diag_score = dp[i - 1][j - 1] + count_columns_score_affine(seq_list1[:, i - 1], seq_list2[:, j - 1], match,
                                                                       mismatch, gap_open, gap_extend)
            I1[i][j] = max(I1[i - 1][j] - gap_extend, dp[i - 1][j] - gap_open)
            I2[i][j] = max(I2[i][j - 1] - gap_extend, dp[i][j - 1] - gap_open)

            scores = (diag_score, I1[i][j], I2[i][j])
            max_direction = np.argmax(scores)

            dp[i][j] = scores[max_direction]
            traceback[i][j] = max_direction

    align_prof1 = [[] for i in range(n1)]
    align_prof2 = [[] for i in range(n2)]
    i, j = l1, l2

    while i > 0 and j > 0:
        if traceback[i][j] == 0:
            add_value_to_profile(align_prof1, seq_list1[:, i - 1])
            add_value_to_profile(align_prof2, seq_list2[:, j - 1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 1:
            add_value_to_profile(align_prof1, seq_list1[:, i - 1])
            add_gap_to_profile(align_prof2)
            i -= 1
        else:
            add_gap_to_profile(align_prof1)
            add_value_to_profile(align_prof2, seq_list2[:, j - 1])
            j -= 1
    while i > 0:
        add_value_to_profile(align_prof1, seq_list1[:, i - 1])
        add_gap_to_profile(align_prof2)
        i -= 1
    while j > 0:
        add_gap_to_profile(align_prof1)
        add_value_to_profile(align_prof2, seq_list2[:, j - 1])
        j -= 1
    profile = list(map(lambda x: ''.join(reversed(x)), align_prof1)) + list(
        map(lambda x: ''.join(reversed(x)), align_prof2))

    return Cluster(profile)


def profile_alignment(Cluster1, Cluster2, match=1, mismatch=1, gap_open=1, gap_extend = None):
    """
            Args:
                Cluster1 (Cluster): the first sequence profile
                Cluster2 (Cluster): the second sequence profile
                match (int): score for matching characters
                mismatch (int): penalty for mismatching characters
                gap_open (int): penalty for opening gaps
                gap_extend (float): penalty for extending gaps

            Returns a new aligned profile (Cluster). Uses function 'profile_alignment_affine' a if gap_extend is specified, otherwise uses function 'profile_alignment_simple'

        """
    if gap_extend is None:
        return profile_alignment_simple(Cluster1, Cluster2,match, mismatch, gap_open)
    else:
        return profile_alignment_affine(Cluster1, Cluster2, match, mismatch, gap_open, gap_extend)