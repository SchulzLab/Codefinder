# Codefinder
_Codefinder_ is a collection of functions which permit the learning of probabilistic nucleotide pairing matrices from sets of motifs. In practice, this code was used on motifs arising from triplexRNA-sequencing and triplexDNA-sequencing in order to learn RNA-DNA pairing probabilities, which were subsequently implemented as substition matrices in the _R_ package _TriplexAligner_.

## Principle
Faced with the challenge of enhancing the prediction of RNA:DNA:DNA triple helix formation by using data obtained from a cellular context (see 10.1093/nar/gky1305), _Codefinder_ was the collection of _R_ functions which arose in the course of designing an Expectation-Maximisation algorithm in order to learn RNA-DNA nucleotide binding affinities. The input for _Codefinder_ are two lists of nucleotide motifs (e.g. output from HOMER, MEME-suite etc.), which one is interested in learning possible binding probabilities between. The output is a user-defined length list of matrices (codes) containing probabilities that each nucleotide may match a nucleotide in the corresponding list, along with metrics describing the relative error of each code given the motif pairings it describes (code objective value).

## Input
As described above, the inputs to _Codefinder_ are the outputs of a motif enrichment algorithm, implemented on nucleotide sequences which may interact. In practice, the input for _Codefinder_ was a list of enriched triplexRNA motifs and a list of enriched triplexDNA motifs, as returned by MEME-ChIP. Such motif files can be parsed with the function _readMotifs()_. We cannot guarantee that this function will work with motif lists returned from other motif enrichment algorithms, and so would reccomend prior conversion to the MEME motif format before using the function.

## Code learning
### Quadratic programming
The functionality of _Codefinder_ revolves around the computation of mapping codes (probability matrices describing the likelyhood that a nucleotide from motif x is paired with a nucleotide in motif y) from pairs of motifs. We accomplished this in _R_ using the _quadprog_ package and its _solve.QP()_ function in order to compute these probabilities between a motif pair given a number of restraints.
### Code objective values
In order to assess the accuracy of a given code for describing a given pairing of motifs, we implmented the _objective()_ function in order to compute a metric we termed the _code objective value_. In practice, this value is the mean column-wise absolute difference between a real motif pair, and an ideal motif pair given the code in question. Where the motifs in a pair are of different lengths, this metric is calculated by sliding the shorter motif partner along the length of the longer partner, and computing the objective value at each position. Subsequently, the lowest value is taken as the objective value for that code across that motif pair. Across multiple motif pairings, the objective value is computed for each pair, and then averaged across the total number of pairings described. It is this mean code objective value which is minimised across the course of the final expectation-maximisation algorithm.
### Assigning motif pairs
Given a code, motifs may be paired with their most ideal partners in a manner which minimises the mean code objective value across the number of pairings. The function _getBestPair()_ carries out this work in _Codefinder_. In practice, this is accomplished by calculating objective values for all possible pairings of motifs, and returning those with the lowest objective values. Note: multiple matching is permitted in this function, meaning some motifs may be excluded from the returned pairings if they are not involved in optimal pairings given the code in question.
### Optimizing a code/motif pairing
By iterating over the previously described steps of learning a code from a given motif pairing, and subsequently identifying the ideal motif pairing given the learned code, the Expectation-Maximization algorithm is established. In this case, it is in fact a minimization, given that the code objective value should decrease with each iteration. As soon as the objective value cannot be minimzed further, the optimized code and motif pairing is output, along with the final code objective value. The starting point for the algorithm requires a random start. This is accomplished by initializing a random pairing of motifs from the two provided lists.
### Learning multiple codes
_Codefinder_ is capable of learning multiple codes from a given set of motifs. In this case, the processes described above are carried out for each desired code. The main point of difference is that when motif pairs are assigned, they are also assigned to the code which describes them best. In practice, this is accomplished by calculating code objective values per motif pairing per code, and the selecting the pairing/code with the lowest value. This results in multiple pairings of motifs, each of which can undergo optimization to learn the ideal code given that pairing.

## Output
As previously stated, the output of _Codefinder_ contains the codes, pairings and objective values which were optimized across the course of the Expectation-Maximization algorithm. Whilst the user may define a number of codes to learn from the input data, it should be noted that if at any point during the opimization process the code has no motif pairs assigned to it, then it is excluded from the analysis. This means that in practice, the user may define the maximum number of codes to be learned, but _Codefinder_ may return fewer final results than that parameter.

## Example

Using the functions contained within _computeCode.R_, you can run an example optimization:

(1) Generate a code
```
reverseComplement = matrix(c(0,0,0,1,
                             0,0,1,0,
                             0,1,0,0,
                             1,0,0,0), 
                             nrow = 4, byrow = T, dimnames = list(c('A', 'C', 'G', 'T'), c('A', 'C', 'G', 'T')))
```

(2) Use the code to generate two motif lists, which may be paired by the code generated in (1)
```
motifLists = createRandomPWMList(alph = 4, lengthRange = 6:8, dp = 2, listLengths = 15, code = reverseComplement, lowEntropy = T)
```

(3) Learn the original code from the motif lists using the Expectation-Maximization algorithm
```

opt = optimize(listR = motifLists$listR, listD = motifLists$listD, unique_pairs = F, dp = 2, plot = F, cores = 8)

# check the results

show(opt$code)
show(opt$pairing)
show(opt$objective)

```


