# Clustal implementation on Python

This is a study project.

Multiple sequence alignments (MSAs) are essential for figuring out evolutionary relationships between biological sequences. The most straightforward approach to computing a MSA results in a complexity of `O(L**N)` for `N` sequences of length `L`. Thus, the heuristics are required. `Clustal` utilises the progressive multiple sequence alignment algorithm. [1] 

The algorithm used in this implementation is based on articles by Desmond G. Higgins and Paul M. Sharp (1988) and Julie D.Thompson, Desmond G.Higgins and Toby J.Gibson (1994). [2, 3]

### Description
In this implementation, the clustering method `UPGMA`[4] is used to build the guide tree (dendrogram). This algorithm creates a rooted tree based on a similarity matrix, so that similar sequences are located next to each other in the tree structure. The values in the similarity matrix for each pair of sequences are calculated using the formula `$dist = 1 - \frac{score}{\max(l_1, l_2)}$`, in which `$score$` is the score of the global alignment of sequences constructed by the Needleman-Wunsch algorithm, `$l_1,l_2$` are the lengths of the sequences. The UPGMA algorithm produces an ultrametric tree in which the distances from the root to every branch tip are equal. The multiple alignment is constructed by first aligning the most similar sequence pairs, then iteratively incorporating more distant sequences through profile-profile alignment, and finally merging existing clusters along a guide tree. 

### Input options

#### Required options

  -f FILENAME, --filename FILENAME
        Input file with sequences to align

  -a {unaligned,aligned}, --alignment_mode {unaligned,aligned}
        Choose the format of sequences: {unaligned, aligned}

  -m {DNA,protein}, --molecule {DNA,protein}
        Choose the type of sequences: {DNA, protein}


To input file with sequences to align, use -- filename option. This program accepts  a sequence file with un-aligned (FASTA format) or aligned sequences (Clustal format). To select one of the modes, use option -- alignment_mode.  The program constucts alignment for both DNA sequences and protein sequences. To align  the  DNA sequence , use option -- molecule DNA. To align  the protein sequence , use option -- molecule protein.

#### Additional options

  --gap-open GAP_OPEN   
        Penalty for gap opening

  --gap-extension GAP_EXTENSION
        Penalty for gap extension


  --match MATCH        
        Bonus for match (for DNA alignment)

  --mismatch MISMATCH   
        Penalty for mismatch (for DNA alignment)

    -h, --help            
        Show this help message and exit






### Literature
[1] Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol Syst Biol. 2011 Oct 11;7:539. doi: 10.1038/msb.2011.75. PMID: 21988835; PMCID: PMC3261699.

[2] Higgins DG, Sharp PM. CLUSTAL: a package for performing multiple sequence alignment on a microcomputer. Gene. 1988 Dec 15;73(1):237-44. doi: 10.1016/0378-1119(88)90330-7. PMID: 3243435.

[3] Thompson JD, Higgins DG, Gibson TJ. CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice. Nucleic Acids Res. 1994 Nov 11;22(22):4673-80. doi: 10.1093/nar/22.22.4673. PMID: 7984417; PMCID: PMC308517.

[4] Sokal, Robert R. “The Principles and Practice of Numerical Taxonomy.” Taxon 12, no. 5 (1963): 190–99. doi: 10.2307/1217562.