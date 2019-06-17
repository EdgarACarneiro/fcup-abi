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

## Searching for Similar Sequences in Databases - _BLAST_

In Bioinformatics to infer the function of an unknown sequence, one has to scan large databases of sequences and compare the query sequence against all the sequences in the database, selecting those with higher degree of similarity.

Pairwise sequence alignment algorithms are not efficient enough (quadratic complexity & memory consuming because of 2 matrixes) when it is necessary to scan very large sets of sequences, so __BLAST__ _- Basic Local Alignment Search Tool -_ is used.

__BLAST__ can be used to infer functional and evolutionary relationships between sequences as well as help to identify members of gene families.

BLAST Concepts:
* __Query sequence__: sequence that will be processed with BLAST (we want to know more about).
* __Target sequence__: sequence in the database that was matched with the query sequence.
* __Dtabase (D)__: database of sequences where the search is done.

BLAST in an heuristic approach based on the idea of _K-indexing_.
_K-mer Indexing_ concepts:
* __Database Pre-processing__: Scan every word of length __K__ and keep it in a hashmap.
* __Query sequence scan__: scan every k-word in Query and get its location in the hashmap.

__BLAST Algorithm__:
1. __Seeding__: find common subwords between the query sequence and the database sequences (using the K-mer Indexing) -> _seeds_.
    * Window scanning of Q to generate K-words: L1-set
    * Find neighborhood words for each k-word until threshold T: L2-set
    * Merge L = L1 U L2
    * Look into H table where L words occur: seeds
2. __Extension__: starting from seeds, extend alignment in both directions -> _high-scoring segment paris (HSP)_.
    * BLAST 2.0. allows alignments with gaps. If two non-overlapping hits are found within distance A of one another on the same diagonal, then merge the hits into an alignment and extend the alignment in both directions until the running alignments score has dropped more than X below the maximum score yet attained.
3. __Evaluation__: assess the statistical significance of each HSP.
    * The __Expect value (E)__ is a parameter that describes the number of hits one can "expect" to see by chance when searching a database of a particular size. It decrease exponentially as the Score (S) of the match increases. The lower the E-value, or the closer it is to zero, the more "significant" the match is.

![Blast Explanation](https://i.imgur.com/nVSUG33.png)

__BLAST__ can be __summarized__ in:
1. Remove regions of low complexity (e.g. sequence repeats) from query sequence.
2. Obtain all possible “words” of size _w_ (a parameter of the algorithm), i.e. sub- sequences of length w occurring in the query sequence;
3. For each word from the previous step, compile the list of all possible words of size w that can be defined in the allowable alphabet, whose alignment score (with no gaps) is higher than a threshold _T_ (parameter of the algorithm);
4. Search in all sequences from the database, all occurrences of the words collected in the last step, which represent __matches (hits or seeds) of size w between the query and one of the database sequences__;
5. Extend all hits from the last step, in both directions, while the score follows a given criterion (typically, the criterion is dependent on the size of the extension);
6. Select the alignments in the previous step with highest scores, normalized for its size (these are named the __high-scoring pairs - HSPs__).

---

## Multiple Sequence Alignment _(MSA)_

Generalization of the Alignment Problem to multiple sequences (__pairwise alignment problem _(PSA)___ for _N > 2_ sequences).

Possible utilities of _MSA_:
* identify recurrent motifs across a family of protein sequences
* study homology relations between an evolving gene (Phylogenetic analysis)
* assess conservation on secondary and tertiary protein structures
* epidemiological studies to understand the mutation rate of a gives species strain.

MSA can be used to start phylogenetic analysis. The phylogenetic tree describes the distance between the sequences under analysis. From the alignment, each column shows if there are conservation of the residues (amino-acids), mutations or divergence from the common ancestor.

Compared to the _PSA_, we now have multiple symbols per column. The __sum of pairs _(SP)___ method considers all possible pairs of symbols within the column and sums them. Eventually, it counts the number of gaps only once.

With many sequences, dynamic programming becomes inviable and therefore __Heuristic approaches__ are used. They represent a trade-off between speed and optimality.
Some Heuristic approaches are:
* __Progressive__: start by aligning the two most similar sequences and iteratively add the other sequences to the alignment.
* __Iterative__: consider an initial alignment and then improve it by adding, removing or moving gaps.
* __Hybrid__: combine both strategies and use complementary information (e.g. protein structural information, other good alignments)

The ___CLUSTAL_ algorithm__ is a classic MSA method that is at the basis of many MSA algorithms. CLUSTAL approach __(Progressive approach)__:
1. Calculate pairwise alignment of all sequences and build a similarity matrix
2. Select the most similar sequences to form the basis of the MSA; the order of the sequences follows that of the guide tree;
3. To add more than two sequences we need to create, from the existing alignment, a summary of the content on each column.
    * __Consensus representation__: each column is represented by the most common character;
    * __Frequency representation__: each column is represented by the frequency of the characters, also called profile.
4. At each iteration align the new sequence with the profile of the current MSA.
    * The column score is the weighted average of the scores of all possible pairs.
    * Once a gap is added to the alignment it will prevail as gap throughout the next alignment steps.
    * In some situations, we may need to join to sub-alignments corresponding to two different sub-branches of the tree. Requires joining two profiles in a generalisation of the above process.

![CLUSTAL Overview](https://i.imgur.com/I3cJTMk.png)

__Cons of Heuristic__ approaches:
* If wrong decisions are made early in the alignment they will be propagated and not corrected afterwards.
* Worst results occur when sequences have low similarity.

---

## Phylogenetic Analysis

__Phylogenetics__ studies the evolutionary history and relationship among individuals or species while __Phylogenetics trees__ illustrate the relationship between these individuals.

__Phylogenetics trees__

* Leaves are sequences (typically from different species or taxonomic categories)
* Internal nodes represent common ancestors of the sequences.
* The structure of a rooted tree may be represented by clusters.
* The height of the nodes in the tree represent a measure of time (moving from the root to the leaves).
* The length of the branch represents the evolutionary distance between the ancestor and the species at the node. This is captured by the number of sequence changes between one level and the next level of the tree.
* Trees without a root illustrate the relation between the leaves without explicitly inferring the common ancestor.

__The Phylogeny Problem__

Assume we are given a set of sequences evolutionarily related, i.e. with a common ancestor.
The problem: __infer the best possible evolutionary tree__.
This is an optimisation problem, since the number of possible trees increases exponentially with the number of input sequences, and therefore requires an objective function. There are three main algorithms:
* __Distance-based algorithm__: compute a distance matrix based on the pairwise distances of the sequences. Derive trees consistent with the distances from the matrix.
* __Maximum parsimony__: search for trees that try to minimize the number of mutations (in internal nodes) to explain the variability of the sequences. Based on _MSA_ of the input sequences. Use certain columns in the alignment that are informative of the possible phylogeny.
* __Statistical/Bayesian__: probabilistic models for the occurrence of different types of mutations in the sequences. Score trees based on their probability searching the most likely trees that explain the sequences according to the assumed model.

__Distance-based methods__

Rely on measuring the consistency of the distances between the leaves in the tree (sequences) and the distances derived from sequence similarity (alignment). The structure of the tree and the length of the branches connecting the nodes reflect the pairwise distances between sequences.

Distance is the reverse of similarity (e.g. percentage of columns in the alignment with mismatches or gaps)

![Score function for Distance-based methods]()
<!-- Missing Upload of image with Imgur -->

* S: set of input sequences
* T: tree
* dij(T): distances of the leaves representing sequences i and j in the the tree.
* Dij: distance between sequences I and j in the input matrix D given from sequence * alignment.

__Ultrametric tree__: the distance between all leaves and the root is the same. The height of each leaf is 0. Thus, _duw = h(w)_ -> _duv = 2 * h(w)_.

As the number of sequences increases the solution space also increases exponentially. Heuristic methods need to be applied to obtain solutions in reasonable time. Hence, use _Unweighted Pair Group Method Using Arithmetic Averages_ - __UPGMA__, which is based on agglomerative hierarchical clustering algorithms.

__Clustering Algorithm__ (uses _Ultrametric Tree_):
* Consider each sequence (tree leaf) as its own cluster. height = 0 in the tree.
* Merge the pair of closest sequence/clusters (minimum value in the matrix D); join these sequences creating an internal node. The __height__ = half of distance between sequences. These sequences form a cluster.
* Distance of a cluster to the remaining sequences is the average of the distances. Update distance matrix D: remove cols and rows of the connected sequences. Add row and col for the new cluster.
* Iteratively: find pairs of clusters with minimum distance and repeat: join clusters, add internal node to the tree with the given height and update D.
* Stop when all sequences are within a single cluster that corresponds to the root of the tree.

![Clustering Algorithm]()
<!-- Missing Upload of image with Imgur -->

---

