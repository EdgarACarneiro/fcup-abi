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


