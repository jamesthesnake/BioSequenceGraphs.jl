## de Bruijn Graph type

A fundamental approach for de-novo gene assembly is to use make use of de Bruijn graphs.
The graph is used to represent fragments of reads (mostly starting with kmers) as vertices and
overlaps between these fragments as edges.
`DeBruijnGraph` type is a special type of `SequenceGraph`. It is also made up of two fields:

- nodes
- links

```
struct DeBruijnGraph
    nodes::Vector{SequenceGraphNode}
    links::Vector{Vector{SequenceGraphLink}}
end
```

We initialize a DeBruijnGraph using the `deBruijn_constructor` function.
This is mainly due to the fact that arbitrary links between two vertices are not allowed in the
de Bruijn graph formalism. The constructor receives as input a list of kmers and generates the deBruijn_Graph
where each kmer is a unique vertex and each overlap  of length $k-1$ is represented with an edge.

```
kmer_vector = generate_random_kmers(DNA,4,10)
dbg = deBruijn_constructor(kmer_vector)

DeBruijnGraph(SequenceGraphNode[SequenceGraphNode{Kmer{DNA,4}}(CGCC, true), SequenceGraphNode{Kmer{DNA,4}}(TCTG, true), SequenceGraphNode{Kmer{DNA,4}}(TGTG, true), SequenceGraphNode{Kmer{DNA,4}}(GAAG, true), SequenceGraphNode{Kmer{DNA,4}}(GGCA, true), SequenceGraphNode{Kmer{DNA,4}}(ACGA, true), SequenceGraphNode{Kmer{DNA,4}}(CGCT, true), SequenceGraphNode{Kmer{DNA,4}}(TCTC, true), SequenceGraphNode{Kmer{DNA,4}}(TACG, true), SequenceGraphNode{Kmer{DNA,4}}(GCAT, true)], Array{SequenceGraphLink,1}[[], [], [], [], [SequenceGraphLink(-5, 10, 1)], [], [], [], [SequenceGraphLink(-9, 6, 1)], []])
```
