# Clustal implementation on Python

This is a study project.

Multiple sequence alignments (MSAs) are essential for figuring out evolutionary relationships between biological sequences. The most straightforward approach to computing a MSA results in a complexity of `O(L**N)` for `N` sequences of length `L`. Thus, the heuristics are required. `Clustal` utilises the progressive multiple sequence alignment algorithm. [1] 

The algorithm used in this implementation is based on articles by Desmond G. Higgins and Paul M. Sharp (1988) and Julie D.Thompson, Desmond G.Higgins and Toby J.Gibson (1994). [2, 3]

### Description
?

### Input options
add parser options later

### Literature
[1] Sievers F, Wilm A, Dineen D, Gibson TJ, Karplus K, Li W, Lopez R, McWilliam H, Remmert M, Söding J, Thompson JD, Higgins DG. Fast, scalable generation of high-quality protein multiple sequence alignments using Clustal Omega. Mol Syst Biol. 2011 Oct 11;7:539. doi: 10.1038/msb.2011.75. PMID: 21988835; PMCID: PMC3261699.

[2] Higgins DG, Sharp PM. CLUSTAL: a package for performing multiple sequence alignment on a microcomputer. Gene. 1988 Dec 15;73(1):237-44. doi: 10.1016/0378-1119(88)90330-7. PMID: 3243435.

[3] Thompson JD, Higgins DG, Gibson TJ. CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice. Nucleic Acids Res. 1994 Nov 11;22(22):4673-80. doi: 10.1093/nar/22.22.4673. PMID: 7984417; PMCID: PMC308517.