from typing import Callable
import numpy as np


def needleman_wunsch(seq1, seq2, match=1, mismatch=1, gap=1):
    '''
    seq1 (str): first sequence to align
    seq2 (str): first sequence to align
    match (int): score for matching characters
    mismatch (int): penalty for mismatching characters
    gap (int): penalty for gaps

    Returns a tuple of two aligned sequences, a score matrix and a final score of the alignment
    '''

    l1, l2 = len(seq1), len(seq2)
    dp = np.zeros((l1 + 1, l2 + 1))
    dp[:,0] = np.linspace(0, -l1 * gap, l1 + 1)
    dp[0,:] = np.linspace(0, -l2 * gap, l2 + 1)
    traceback = np.zeros((l1 + 1, l2 + 1))

    for i in range(1, l1+1):
        for j in range(1, l2+1):

            if seq1[i-1] == seq2[j-1]:
                diag_score = dp[i-1][j-1] + match
            else:
                diag_score = dp[i-1][j-1] - mismatch
            up_score = dp[i-1][j] - gap
            left_score = dp[i][j-1] - gap

            scores = (diag_score, up_score, left_score)
            max_direction = np.argmax(scores)
            
            dp[i][j] = scores[max_direction]
            traceback[i][j] = max_direction
            
    align_seq1, align_seq2 = [], []
    i, j = l1, l2

    while i > 0 and j > 0:
        if traceback[i][j] == 0:
            align_seq1.append(seq1[i-1])
            align_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 1:
            align_seq1.append(seq1[i-1])
            align_seq2.append('-')
            i -= 1
        else:
            align_seq1.append('-')
            align_seq2.append(seq2[j-1])
            j -= 1

    align_seq1 = ''.join(reversed(align_seq1))
    align_seq2 = ''.join(reversed(align_seq2))
    
    return align_seq1, align_seq2, dp, dp[l1][l2]


def needleman_wunsch_affine(seq1, seq2, match=1, mismatch=1, gap_open=1, gap_extend=0.5):
    '''
    seq1 (str): first sequence to align
    seq2 (str): first sequence to align
    match (int): score for matching characters
    mismatch (int): penalty for mismatching characters
    gap_open (int): penalty for opening gaps
    gap_extend (float): penalty for extending gaps

    Returns a tuple of two aligned sequences, a score matrix and a final score of the alignment
    '''

    l1, l2 = len(seq1), len(seq2)
    dp = np.zeros((l1 + 1, l2 + 1))
    dp[1:,0] = np.linspace(-gap_open, -(l1-1) * gap_extend, l1)
    dp[0,1:] = np.linspace(-gap_open, -(l2-1) * gap_extend, l2)

    I1, I2 = np.zeros((l1 + 1, l2 + 1)), np.zeros((l1 + 1, l2 + 1))
    I1[1:,0] = np.linspace(-gap_open, -(l1-1) * gap_extend, l1)
    I2[0,1:] = np.linspace(-gap_open, -(l2-1) * gap_extend, l2)

    traceback = np.zeros((l1 + 1, l2 + 1))

    for i in range(1, l1+1):
        for j in range(1, l2+1):
            if seq1[i-1] == seq2[j-1]:
                diag_score = dp[i-1][j-1] + match
            else:
                diag_score = dp[i-1][j-1] - mismatch
            I1[i][j] = max(I1[i-1][j] - gap_extend, dp[i-1][j] - gap_open)
            I2[i][j] = max(I2[i][j-1] - gap_extend, dp[i][j-1] - gap_open)

            scores = (diag_score, I1[i][j], I2[i][j])
            max_direction = np.argmax(scores)
            
            dp[i][j] = scores[max_direction]
            traceback[i][j] = max_direction
            
    align_seq1, align_seq2 = [], []
    i, j = l1, l2

    while i > 0 and j > 0:
        if traceback[i][j] == 0:
            align_seq1.append(seq1[i-1])
            align_seq2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 1:
            align_seq1.append(seq1[i-1])
            align_seq2.append('-')
            i -= 1
        else:
            align_seq1.append('-')
            align_seq2.append(seq2[j-1])
            j -= 1
    
    while i > 0:
        align_seq1.append(seq1[i - 1]) 
        align_seq2.append('-')
        i -= 1
    while j > 0:
        align_seq1.append('-')
        align_seq2.append(seq2[j - 1])
        j -= 1

    align_seq1 = ''.join(reversed(align_seq1))
    align_seq2 = ''.join(reversed(align_seq2))
    
    return align_seq1, align_seq2, dp, dp[l1][l2]


def base_needleman_wunsch_affine(
    seq1_len: int,
    seq2_len: int,
    weight_gen: Callable[[int, int], float],
    gap_open: float = 1.0,
    gap_extend: float = 0.5,
):
    '''
    seq1_len: (int) length of first sequence
    seq2_len: (int) first sequence to align
    weight_gen: (callable) universal weight generator
    gap_open (int): penalty for opening gaps
    gap_extend (float): penalty for extending gaps

    Returns a tuple of two aligned sequences, a score matrix and a final score of the alignment
    '''

    l1, l2 = seq1_len, seq2_len
    dp = np.zeros((l1 + 1, l2 + 1))
    dp[1:,0] = np.linspace(-gap_open, -(l1-1) * gap_extend, l1)
    dp[0,1:] = np.linspace(-gap_open, -(l2-1) * gap_extend, l2)

    I1, I2 = np.zeros((l1 + 1, l2 + 1)), np.zeros((l1 + 1, l2 + 1))
    I1[1:,0] = np.linspace(-gap_open, -(l1-1) * gap_extend, l1)
    I2[0,1:] = np.linspace(-gap_open, -(l2-1) * gap_extend, l2)

    traceback = np.zeros((l1 + 1, l2 + 1))

    for i in range(1, l1 +1 ):
        for j in range(1, l2 + 1):
            diag_score = dp[i - 1][j - 1] + weight_gen(i - 1, j - 1)
            I1[i][j] = max(I1[i - 1][j] - gap_extend, dp[i - 1][j] - gap_open)
            I2[i][j] = max(I2[i][j - 1] - gap_extend, dp[i][j - 1] - gap_open)

            scores = (diag_score, I1[i][j], I2[i][j])
            max_direction = np.argmax(scores)
            
            dp[i][j] = scores[max_direction]
            traceback[i][j] = max_direction
            
    align_seq1, align_seq2 = [], []
    i, j = l1, l2

    while i > 0 and j > 0:
        if traceback[i][j] == 0:
            align_seq1.append(i - 1)
            align_seq2.append(j - 1)
            i -= 1
            j -= 1
        elif traceback[i][j] == 1:
            align_seq1.append(i - 1)
            align_seq2.append('-')
            i -= 1
        else:
            align_seq1.append('-')
            align_seq2.append(j - 1)
            j -= 1
    
    while i > 0:
        align_seq1.append(i - 1) 
        align_seq2.append('-')
        i -= 1
    while j > 0:
        align_seq1.append('-')
        align_seq2.append(j - 1)
        j -= 1

    align_seq1 = list(reversed(align_seq1))
    align_seq2 = list(reversed(align_seq2))

    return align_seq1, align_seq2, dp, dp[l1][l2]
