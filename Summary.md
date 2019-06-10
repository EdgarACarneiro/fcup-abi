# Contents Summary

## Molecular Biology Concepts

BioinformaticsAlgorithms.C2a.pdf

---


## Basic Processing of Biological Sequences

BioIinformaticsAlgorithms.C3.pdf

---

## Finding patterns in Sequences

Examples where finding patterns might be useful:
* _mRNA_ contain in their promoter region special signal for the binding of the transcription factors that regulate gene expression.
* Proteins contain sub-sequences, also called domains, that determine the function of the protein, including the binding sites for ligands or structural domains that determine the structure of the protein.

__Sequence motifs__ are relatively short sub-sequences, shared among several (related) sequences that are presumed to have a similar biological function.

_Hamming distance_ between two sequences, with the same size, are the number of different characters between the 2 sequences.

Improved search rules, that the _Booyer-Moore_ algorithm incldues:
* __Bad-Character Rule:__ The search can advance the pattern to the next ocurrence of the symbol in the sequence at the position of the mismatch.
* __Good Suffix Rule:__ In case of a mismatch we can move forward to the next instance in the pattern of the part (suffix) that matched before of the mismatch.

___Prosite___ is a database of protein families and domains. They analyze sets of related sequences to identify regions of sequence similarity. Prosite patterns are described using a syntax similar to the RE syntax:

```python
# Convert a Prosite Pattern to a Regex Pattern
def prosite_to_regex(prosite):
    converter = {
        '-': '',
        'x': '.',
        '(': '{',
        ')': '}',
        '{': '[^',
        '}': ']'
    }
```

__Restriction enzymes__ are proteins that cut the DNA in regions that contain specific sub- sequences or motifs. They are a very useful tool in molecular biology as they allow to create restriction maps useful for cloning and sequencing techniques.

__IUB code__: way of representing ambiguous sequences, using other characters.
```python
iub_to_re = {
    "A":"A",
    "C":"C",
    "G":"G",
    "T":"T",
    "R":"[GA]",
    "Y":"[CT]",
    "M":"[AC]",
    "K":"[GT]",
    "S":"[GC]",
    "W": "[AT]",
    "B":"[CGT]",
    "D":"[AGT]",
    "H":"[ACT]",
    "V":"[ACG]",
    "N":"[ACGT]"
}
```

---

## Pairwise Sequence Alignment

Bioinformatics relies on the assumption that biological sequences with __high similarity__ share __similar functions__.

__Phylogenetics__ is the division of biology that studies evolutionary divergence and relationship between organisms, based on _Similarity_ (resemblance and differences between organisms) and _Homology_ (common ground between the organisms and if they share any ancestral characteristics).

The __edit distance__ between 2 sequences is the minimum number of editing operations to transform one sequence into the other. Operations:
* __Substitution__: Replace symbols - `T -> C`.
* __Deletion__: Delete 1 symbol - `C -> _`.
* __Insertion__: Insert 1 symbol - `_ -> C`.

__Sequence alignment__ is a way of arranging the biological sequences to identify regions of similarity. Such regions may indicate a functional, structural, or evolutionary relationships between the sequences.

A __global alignment__ aligns two sequences from beginning to end, aligning each letter in each sequence only once. An alignment is always produced. _(Needleman-Wunsch algorithm)_

A __local alignment__ maximizes the alignment of the parts of the sequences that share similarity. Finds the best aligned subsequence. An alignment may not be produced if no sufficient similarity is found. _(Smith-Waterman algorithm)_

Substitution matrices give a score for each substitution of one amino-acid by another. BLOSUM (BLOcks SUbstitution Matrix) matrix is a substitution matrix used for sequence alignment of proteins.

In __pairwise sequence alignment__ we try to arrange two sequences so that the number of matching characters is maximised.

__Dynamic Programming for Global Alignment:__ fill __score matrix__ cell by cell, using the value of adjacent cells to reach the target cell.
```python
S[i][j] = max(S[i-1][j-1] + sm(a[i], b[j]), S[i-1][j] + g, S[i][j-1] + g) for all 0 < i <= n and 0 < j <= m
```
Matrix is filled left to right and top to bottom. To calculate `S[i][j]` you need to have
calculated: `S[i-1][j], S[i][j-1], S[i-1][j-1]`.

Additionally, keep a __trace-back matrix__ (T) to keep all the possible optimal moves at each cell. From T, one can recover the optimal alignment - start from lower-right cell and trace-bacl to upper-left cell.

__Dynamic Programming for Local Alignment:__ find the best partial alignment of the sub-sequences from the 2 sequences.
```python
S[i][j] = max(S[i-1][j-1] + sm(a[i], b[j]), S[i-1][j] + g, S[i][j-1] + g, 0) for all 0 < i <= n and 0 < j <= m
```
For the optimal alignment, one now starts in the cells with highest score.

Now, the __trace-back matrix T contains 4 possible values__: three previous values and an extra value to the cases where the alignment is terminated (correspond to cells in S with 0).

---

