__precompile__()

module BioSequenceGraphs

export SequenceDistanceGraph,
    new_graph_from_kmerlist

using BioSequences

#include("Nodes.jl")
#include("Links.jl")
#include("SequenceGraph.jl")
#include("IO.jl")
include("SequenceDistanceGraph.jl")
include("graph_building.jl")
#include("DeBruijnGraph.jl")
#include("GFA1/GFA1.jl")
end # module BioSequenceGraphs
