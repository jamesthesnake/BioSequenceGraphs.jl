__precompile__()

module BioSequenceGraphs

using BioSequences

include("Nodes.jl")
include("Links.jl")
include("SequenceGraph.jl")
include("IO.jl")
include("seq-distance-graphs/SequenceDistanceGraph.jl")
#include("DeBruijnGraph.jl")
#include("GFA1/GFA1.jl")
end # module BioSequenceGraphs
