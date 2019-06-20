# Contents Summary

## Molecular Biology Concepts

Evolution has three main components:
* __Inheritance__
    * Passage of characteristics from parents to offspring.
    * Determines most of the structure and functions of the organism.
    * The amount of variation passed from one generation to the next one is very small.
* __Variation__
    * Occurs with mutations, sexual recombination, random changes of genetic
material.
* __Selection__
    * Reflects the fitness of the organism to adapt to the medium.
    * This fitness capacity is expressed through reproduction, i.e. parents transmit their fitness to their offspring.

__Phylogenetics__ is the division of biology that studies evolutionary divergence and relationship between organisms, based on two important important concepts:
* __Similarity__: Measures the resemblance and differences between organisms without taking
into account any contextualization.
* __Homology__: Investigates the common ground between the organisms and if they share any
ancestral characteristics. Find the point in the evolutionary tree that they started to diverge.

> The Micro perspective: The cell

The cell is the unit of life. Each cell derives from another cell and contains all the necessary information to replicate itself. All cells have some common features, such as its composition (70% water, 23 % macro-molecules and 7% small molecules).

Organisms are categorized according to their cell type:
* __Prokaryotes__
    * No nucleus or internal membranes.
    * Lacks other membrane-bound organelles.
* __Eukaryotes__
    * Nucleus (membrane-enclosed DNA).
    * Internal membranes.
    * Organelles inside the cell that play different and specific roles.

And according to the number of cells:
* __Unicellular__
    * _Prokaryotes_: Bacteria, Archaea.
    * _Eukaryotes_: baker yeast.
* __Multicellular__
    * _Eukaryotes_: animals, plants, fungi

> What defines a cell?

* __Proteins__ perform most of the functions in the cell, they have catalytic and structural functions, from sensors and signaling to promoting chemical reactions (catalysis).
* __Enzymes__ are proteins that convert cellular molecules in other types of molecules necessary for the functions of the cell, like generating energy.

__DNA and RNA are composed of nucleic acids. Proteins are composed of amino- acids.__

Other important organic molecules:
* __Carbohydrates__, store energy (simple - immediate energy demands, complex - long term
storage of energy).
* __Lipids__, make part of the plasma membrane and also store energy and are involved in
signaling.

Other components of the cell:
* __Mitochondria__ and Chloroplasts are cellular organelles involved in the production of energy.
* __Ribosomes__ are large and complex molecules composed by a mixture of proteins and genetic
material. Their function is to assemble proteins.

Cells form tissues that themselves form organs, and eventually entire organisms.

> Information transfer in the cell: nucleic acids

Both DNA and RNA are polymers composed of four nucleic acid units, called __nucleotides__ or bases.
* Adenyne (__A__) and Guanine (__G__), belong to one group (purynes).
* Cytosine (__C__) and Timine (__T__) and Uracil (__U__), belong to another group (pyrimidine).

Timine only exists in DNA and Uracil is only found in RNA, the other three bases exist in both.

The DNA is composed of two complementary strands, in a double-helix structure, due to connections established between the bases in both strands.
* __Adenine and Timine__ (A == T), connected by two hydrogen connections
* __Guanine and Cytosine__ (G === C), connected by three hydrogen connections

Chains are __antiparallel__ because they are connected in opposite directions.

> Genetic material

* __Genome__: an organism’s genetic material (complete set of DNA).
    * human genome has 24 distinct chromosomes.
    * Each chromosome contains many genes.
    * In __Prokaryotics__ existis in the form of a circular chromosome located in the cytoplasm.
    * In __Eukaryotes__ is found in the nucleus and is tightly packaged into linear chromosomes.
* __Gene__: a discrete units of hereditary information located on the chromosomes and consisting of DNA and encode instructions on how to make proteins.
* __Genotype__: The genetic makeup of an organism.
* __Phenotype__: the physical expressed traits of an organism.

![Genetic material](https://i.imgur.com/xV6K7jh.png)

> Protein Synthesis

Cellular DNA contains instructions for building the various proteins the cell needs to survive:
1. To manufacture these proteins, specific genes within its DNA must first be transcribed into molecules of mRNA (__Transcription__);
2. these transcripts must be translated into chains of amino acids (__Translation__),
3. fold into fully functional proteins.

In most eukaryotic genes, coding regions (__exons__) are interrupted by noncoding regions (__introns__). During transcription, the entire gene is copied into a pre-mRNA, which includes exons and introns. During the process of RNA splicing, introns are removed and exons joined to form a contiguous coding sequence. This "mature" mRNA is ready for translation.

All cells in a multicellular organism contain the same set of genetic information
(Genome).

The differences in the abundance of the RNA (Transcriptome) determines the cell specificity.

> Alternative Splicing

__Alternative splicing__ (AS), the process in which the exons of the pre-RNAs are spliced in different combinations to produce distinct mRNA that lead to structurally and functionally protein variants.

Combinatorial splicing leads to the generation of of multiple isoforms from a single gene.

---


## Basic Processing of Biological Sequences

> Open Reading Frames

The translation of a protein sequence occurs for the coding region of the gene. This region start with a __start codon (ATG)__ and stops when one of the __stop codons__ is found.

A __reading frame__ is a way of dividing a DNA (or RNA) sequence into a set of consecutive non-overlapping triplets or codons. Recall that a sequence may have 6 reading frames: 3 in one strand: +1, +2, +3 starting at position 1, 2 and 3 respectively of the sequence (in python strings at index 0, 1 and 2).

An __open reading frame__ is a reading frame with the potential to be translated into protein (sarting with start codon - __M, Meteonin__ - and the stop codon - __A__).

__Putative proteins__ are the result of analysing all open reading frames of a given DNA sequence.

---

## Finding patterns in Sequences

Examples where finding patterns might be useful:
* _mRNA_ contain in their promoter region special signal for the binding of the transcription factors that regulate gene expression.
* Proteins contain sub-sequences, also called domains, that determine the function of the protein, including the binding sites for ligands or structural domains that determine the structure of the protein.

__Sequence motifs__ are relatively short sub-sequences, shared among several (related) sequences that are presumed to have a similar biological function.

__Hamming distance__ between two sequences, with the same size, are the number of different characters between the 2 sequences.

Example, where Hamming distance is 3.
```json
Q:ATTACGAT
    | | |
T:ATCAGGTT
```

Improved search rules, that the _Booyer-Moore_ algorithm includes:
* __Bad-Character Rule:__ The search can advance the pattern to the next occurrence of the symbol in the sequence at the position of the mismatch.
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

> __Dynamic Programming for Global Alignment:__

Fill the __score matrix__ cell by cell, using the value of adjacent cells to reach the target cell.
```python
S[i][j] = max(S[i-1][j-1] + sm(a[i], b[j]), S[i-1][j] + g, S[i][j-1] + g) for all 0 < i <= n and 0 < j <= m
```
Matrix is filled left to right and top to bottom. To calculate `S[i][j]` you need to have
calculated: `S[i-1][j], S[i][j-1], S[i-1][j-1]`.

* Iteration example 0:

```python
gap = -8
seq1 = 'PHSWG'
seq2 = 'HGWAG'
```

|     | gap |  H |  G  |  W  |  A  |  G  |
|:---:|:---:|:--:|:---:|:---:|:---:|:---:|
| gap |   0 | -8 | -16 | -24 | -32 | -40 |
|   P |  -8 |    |     |     |     |     |
|   H | -16 |    |     |     |     |     |
|   S | -24 |    |     |     |     |     |
|   W | -32 |    |     |     |     |     |
|   G | -40 |    |     |     |     |     |

Additionally, keep a __trace-back matrix__ (T) to keep all the possible optimal moves at each cell. From T, one can recover the optimal alignment - start from lower-right cell and trace-back to upper-left cell.

> Trace-back example

Assuming one has the trace-back matrix:

```python
DIAGONAL = 1
VERTICAL = 2
HORIZONTAL = 3
```

|     | gap | H | G | W | A | G |
|:---:|:---:|:-:|:-:|:-:|:-:|:-:|
| gap |   0 | 3 | 3 | 3 | 3 | 3 |
|   P |   2 | 1 | 1 | 3 | 1 | 3 |
|   H |   2 | 1 | 1 | 1 | 1 | 1 |
|   S |   2 | 2 | 1 | 1 | 1 | 3 |
|   W |   2 | 2 | 2 | 1 | 3 | 3 |
|   G |   2 | 2 | 1 | 2 | 1 | 1 |

We start at `TB[G][G]`. There we have a 1, so we must proceed in the diagonal to `TB[W][A]` so our sequences look like:
```
G
G
```
Remember that since we are starting from the  bottom right and going to the upper left in a trace-back style, we must also fill our alignment sequences in the opposite natural order.

Next, we find a 3 so we must move in the horizontal and reach `TB[W][W]`. Since we only decrease our column without decreasing the row, we are only moving through the second sequence. Our first sequence of the alignment will get a __gap__.
```
-G
AG
```
Next we find 1 so we reach `TB[S][G]`.
```
W-G
WAG
```
Next we have a 1 so we reach `TB[H][H]`.
```
SW-G
GWAG
```
Next we find 1 so we reach `TB[P][gap]`.
```
HW-G
HWAG
```
Next we find a 2, so we must move in the vertical and reach `TB[gap][gap]`. Since we only moved in the first sequence, the second sequence will get a __gap__.
```
PHW-G
-HWAG
```

We have reached `TB[0][0]` and so we have finished!


> __Dynamic Programming for Local Alignment:__

Find the best partial alignment of the sub-sequences from the 2 sequences.
```python
S[i][j] = max(S[i-1][j-1] + sm(a[i], b[j]), S[i-1][j] + g, S[i][j-1] + g, 0) for all 0 < i <= n and 0 < j <= m
```
For the optimal alignment, one now starts in the cells with highest score.

Now, the __trace-back matrix T contains 4 possible values__: three previous values and an extra value to the cases where the alignment is terminated (correspond to cells in S with 0).

> Trace-back example

Assuming we have the score matrix:

|     | gap | H | G |  W |  A |  G |
|:---:|:---:|:-:|:-:|:--:|:--:|:--:|
| gap |   0 | 0 | 0 |  0 |  0 |  0 |
|   P |   0 | 0 | 0 |  0 |  0 |  0 |
|   H |   0 | 8 | 0 |  0 |  0 |  0 |
|   S |   0 | 0 | 8 |  0 |  1 |  0 |
|   W |   0 | 0 | 0 | 19 | 11 |  3 |
|   G |   0 | 0 | 6 | 11 | 19 | 17 |

And we have the correspondent trace-back matrix:

|     | gap | H | G | W | A | G |
|:---:|:---:|:-:|:-:|:-:|:-:|:-:|
| gap |   0 | 0 | 0 | 0 | 0 | 0 |
|   P |   0 | 0 | 0 | 0 | 0 | 0 |
|   H |   0 | 1 | 0 | 0 | 0 | 0 |
|   S |   0 | 0 | 1 | 0 | 1 | 0 |
|   W |   0 | 0 | 0 | 1 | 3 | 3 |
|   G |   0 | 0 | 1 | 2 | 1 | 1 |

We start at the position of the cell with __maximum score__. In this case there are 2 cells with maximum score, so we will end up with 2 possible best alignments.

* For _alignment\_one_, we start at `TB[G][A]` nd we have a 1 so we move diagonally to `TB[W][W]`. There we have a 1 as well so we move to `TB[S][G]`. Once again a 1, to `TB[H][H]`. Again a 1 to `TB[P][gap]` were we have a 0, so we stop. We end up with the following alignment evolution:
```
G
A
```
```
WG
WA
```
```
SWG
GWA
```
```
HSWG
HGWA
```

* For _alignment\_two_ we start at `TB[W][W]` and then, form there, its actually the exact same alignment as _alignment\_one_, so:
```
HSW
HGW
```

---

## Searching for Similar Sequences in Databases - _BLAST_

In Bioinformatics to infer the function of an unknown sequence, one has to scan large databases of sequences and compare the query sequence against all the sequences in the database, selecting those with higher degree of similarity.

Pairwise sequence alignment algorithms are not efficient enough (quadratic complexity & memory consuming because of 2 matrixes) when it is necessary to scan very large sets of sequences, so __BLAST__ _- Basic Local Alignment Search Tool -_ is used.

__BLAST__ can be used to infer functional and evolutionary relationships between sequences as well as help to identify members of gene families.

BLAST Concepts:
* __Query sequence__: sequence that will be processed with BLAST (to know more about).
* __Target sequence__: sequence in the database that was matched with the query sequence.
* __Database (D)__: database of sequences where the search is done.

BLAST is an heuristic approach based on the idea of _K-indexing_.

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

> Consideration

Notice however, that the teacher _MyBlast_ implementation only creates a hash table out of the query sequence - __build_map()__ function. Regarding the database, he iterates over it, calling the __get_hits()__ method for each stored sequence.

> Exercises

Original __get_hits()__ function:
```python
def get_hits (seq, m, w):
    res = [] # list of tuples

    for i in range(len(seq)-w+1):
        subseq = seq[i:i+w]

        if subseq in m:
            l = m[subseq]

            for ind in l:
                res.append((ind, i))

    return res
```

__get_hits()__ allowing at most 1 mismatch, meaning returns all hits that have _w_ or _w-1_ matches:
 ```python
def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2), "Sequences can not have different lengths"

    return sum(
        [1 for i in range(len(seq1))
         if seq1[i] != seq2[i]]
    )

def get_hits (seq, m, w):
    res = [] # list of tuples

    for i in range(len(seq)-w+1):
        subseq = seq[i:i+w]

        for word in m.keys():
            if haming_distance(subseq, word) <= 1:

            for ind in m[word]:
                if (ind, i) not in res: # No duplicates
                    res.append((ind, i))

    return res
```

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

![Score function for Distance-based methods](https://i.imgur.com/WBuKTvF.png)

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

![Clustering Algorithm](https://i.imgur.com/h3fs6c6.png)

> Practical example of _UPGMA_:

Having the sequences:
```js
s1 = A-CATATC-AT-
s2 = A-GATATT-AG-
s3 = AACAGATC-T--
s4 = G-CAT--CGATT
```

Distance Matrix (using __Hamming Distances__):

* __Iteration 0__

|   | s1 | s2 | s3 | s4 |
| -:|:--:|:--:|:--:|:--:|
| s1|  0 |    |    |    |
| s2|  3 |  0 |    |    |
| s3|  4 |  6 |  0 |    |
| s4|  5 |  8 |  9 |  0 |

The most similar sequences are the __s1 & s2__, and so we start our tree, with the height being `3/2`:

```js
     ___
1.5 |   | 1.5
    |   |
   s1   s2
```

Now, lets compute the new distance matrix, remembering that:

![Imgur](https://i.imgur.com/olswI4A.png)

|A| represents the cluster size of A, e.g. |_s1, s2_| = 2 and |_s1, s2, s3, s4_| = 4.

* __Iteration 1__

* d(_s1 U s2_), _s3_ = (1 * 4 + 1 * 6) / (1 + 1) = 5
* d(_s1 U s2_), _s4_ = (1 * 5 + 1 * 8) / (1 + 1) = 6.5

|       | s1, s2 | s3 | s4 |
| -----:|:------:|:--:|:--:|
| s1, s2|      0 |    |    |
|     s3|      5 |  0 |    |
|     s4|    6.5 |  9 |  0 |

The most similar sequences are the __s1,s2 & s3__, and so we update our tree, with the height being `5/2`:

```js
        _____________
  2.5  |             |
     __|__           | 2.5
1.5 |     | 1.5      |
    |     |          |
   s1     s2         s3
```

* __Iteration 2__

* d(_(s1, s2) U s3_), _s4_ = (2 * 6,5 + 1 * 9) / (2 + 1) = 7.33

|           | s1, s2, s3 | s4 |
| ---------:|:----------:|:--:|
| s1, s2, s3|          0 |    |
|         s3|       7.33 |  0 |

The only available sequences are the __s1,s2 & s3__, and so we update our tree, with the height being `7.33/2`:

```js
                         |
               __________|________
         3.66 |                   |
              |                   |
        ______|______             |
  2.5  |             |            | 3.66
     __|__           | 2.5        |
1.5 |     | 1.5      |            |
    |     |          |            |
   s1     s2         s3           s4
```

---

## Graphs and Biological Networks

The elements of the cell form networks that are significantly different from random networks. Examples of cell networks include:
* Protein-protein interactions;
* Metabolic;
* Signalling and Metabolic networks;
* Protein phosporylation*;
* Genetic interactions;
* Co-expression networks;
* Protein-DNA interactions;

The cell is formed by the interplay of networks and emerge as the sum of the interaction of its different elements. It is a network of networks.

Graphs can be:
* __Undirected__ if edges are unordered
* __Directed__ or __digraph__ if edges have an orientation, i.e. pairs are ordered
* __Weighted__ if numerical weights are associated to edges.

__Matrices__ or __adjacency lists__ are the typical way to represent graphs.

__Matrices Representation__

All possible combination of vertices are represented.
* __Undirect graph__: the rows and cols represent the nodes. _Cell(i,j)_ represents an edge between node i and j: 1 if connected and 0 otherwise.
* __Direct graph__: rows represent the origin node and cols the destination node. _Cell(i,j)_ represents an edge between node i and j: 1 if connected and 0 otherwise.
* __Weighted graph__: Cell(i,j) represents the weight of the edge between node i and j. 0 if not connected.

__Adjacency lists Representation__

Represent only the existing edges. Each vertex is linked to a list with associated neighbor nodes.
* __Direct graph__: list for vertex v includes all destination nodes for the edges where v is the origin.
* __Undirect graph__: the edge only exist for one of the directions.

__Some Network concepts:__
* __Degree:__ how many links a node has.
* __Degree In:__ ho many links a node has entering it.
* __Degree Out:__ ho many links a node has leaving it.
* __Shortest path__: path with the smallest number of links between selected nodes. Distances are measured as the length of the path.
* __Degree distribution P(k)__: gives the probability that a selected node has exactly k links. P(k) is obtained by counting how many nodes have k = 1,2,... links and dividing by the total number of nodes N. The degree distribution will allow to discriminate between different types of networks.
* __Clustering coefficient__: indicates how connected are the the nodes on the network. If a node A is connected to B and B is connected to C, then it is highly probable that A is connected to C. The clustering coefficient of a node, gives the number of triangles that go through that node

The last two measures capture the features and characteristics of the network and therefore allow to compare and classify the various networks.

Different type of networks according to degree distribution:
* __Scale free network__: the degree distribution approximates a power-law distribution where P(k) ~ k-a
    * Most networks on the cell are scale-free.
* __Random network__: P(k) decreases exponentially, which indicates that nodes that deviate from average are extremely rare.
* __Hierarchical network__:
    * The starting point is a small number of four nodes highly linked.
    * Communication between highly clustered models is maintained by hubs.
    * Scale free topology + modular structure
    * Clustering coefficient follows C(k) ~ K-1 is the most important signature for this type of networks.

Two strategies to traverse the graph:
* __Breadth-First Search (BFS)__: starts with the source node, then visits all its successors, followed by their successors until all nodes are visited.
* __Depth-First Search (DFS)__: starts with the source node, explores its first successor, then its first successor until no further exploration is possible. Backtracks to explore further alternatives.

Not all nodes are equally significant. Identification of nodes that are most important is a central task in network analysis.

Identification of hub nodes provides information on important properties of the network, like the robustness and resistance to failure. In biological networks a hub node in a gene-gene or gene-protein interaction network may reveal an important regulator, for instance a transcription factor that regulates many other genes.

---

## High-throughput sequencing applications

> Sequencing Technologies

Presents in slides the typical applications of massive (also called next-generation or high-throughput) sequencing technologies. 

__DNA Sequencing__: process of determining the order of nucleotides in DNA.

### Exercises

> Exercise #1: __RPKM__ (Reads per KiloBase per Million):

Used for single end RNA-sequencing.

Notice that:
* Sequencing runs with more depth will have more reads mapping to each gene.
* Longer genes will have more reads mapping to them.

__RPKM__ = 10^9 * #reads / ((Σ #reads) * Len)

| gene | Len | Reads | __RPKM__ |
|:----:|:---:|:-----:|:--------:|
|   A  |   2 |    20 |          |
|   B  |   3 |    10 |          |
|   C  |   5 |    50 |          |
|   D  |  10 |   100 |          |

We can infer that gene A's __sequencing runs have more depth__ than gene B, since they have similar gene size but A #reads is _2x_ bigger.

__Σ #reads__ = 20 + 10 + 50 + 100 = 180

_Example A_ <=> _RPKM_ = 20 * 10^9/ (180 * 2) = 5 * 10^7

| gene | Len | Reads | __RPKM__ |
|:----:|:---:|:-----:|:--------:|
|   A  |   2 |    20 |   5.5    |
|   B  |   3 |    10 |   1.8    |
|   C  |   5 |    50 |   5.5    |
|   D  |  10 |   100 |   5.5    |

> Exercise #3: Candidate mutations in sequencing read pile-up

```json
CGACGACGACGACGAATGATGTATTATCGAGCGAGCGGCAGATGCTA
-----------------------------------------------
CGACGACGACGACGAATGAT    TATCGAGCGCGCGGCAGATG
 GACGACGACGACGAATGATG      CGAGCGCGCGGCAGATGCTA
 GACGACGACGACGAACGATG      CGAGCGCGCGGCAGATGCTA
         CGACGAACGATGTATTATCG
            CGAACGATGTATTATCGAGC
                   TGTATTATCGAGCGCGCGGC
```

Reads that do not match the reference genome:

```json
CGACGACGACGACGAATGATGTATTATCGAGCGAGCGGCAGATGCTA
-----------------------------------------------
CGACGACGACGACGAATGAT    TATCGAGCGCGCGGCAGATG
 GACGACGACGACGAATGATG      CGAGCGCGCGGCAGATGCTA
 GACGACGACGACGAACGATG      CGAGCGCGCGGCAGATGCTA
         CGACGAACGATGTATTATCG
            CGAACGATGTATTATCGAGC
                   TGTATTATCGAGCGCGCGGC
                ^                ^
```

