struct SequenceGraph
    nodes::Vector{SequenceGraphNode}
    links::Vector{Vector{SequenceGraphLink}}
end

nodes(sg::SequenceGraph) = sg.nodes
node(sg::SequenceGraph, i::NodeID) = nodes(sg)[abs(i)]
links(sg::SequenceGraph) = sg.links



"""
    links(sg::SequenceGraph, node::NodeID)

Get all of the links of a Node of a sequence graph.
"""
function links(sg::SequenceGraph, node::NodeID)
    l = links(sg)
    @inbounds return l[abs(node)]
end


"""
    add_node!(sg::SequenceGraph, n::SequenceGraphNode)

Add a SequenceGraphNode `n` to a SequenceGraph `sg`.

Returns the node ID used to access the node added from the graph.
"""
function add_node!(sg::SequenceGraph, n::SequenceGraphNode)
    newlen = length(push!(nodes(sg),n))
    resize!(links(sg), newlen)
    return newlen
end

function remove_node!(sg::SequenceGraph, n::NodeID)
    oldlinks = copy(links(sg, n))
    for oldlink in oldlinks
        remove_link!(sg, source(oldlink), destination(oldlink))
    end
end


"""
    add_link!(sg::SequenceGraph, source::NodeID, dest::NodeID, dist::Int)

Construct a link between two nodes in a sequence Graph.
"""
function add_link!(sg::SequenceGraph, source::NodeID, dest::NodeID, dist::Int)
    push!(links(sg, source), Link(source, dest, dist))
    push!(links(sg, dest), Link(dest, source, dist))
end

function remove_link!(sg::SequenceGraph, source::NodeID, dest::NodeID)
    slinks = links(sg, source)
    slinkslen = length(slinks)
    filter!(!isequal(SequenceGraphLink(source, dest, 0)), slinks)
    dlinks = links(sg, dest)
    dlinkslen = length(dlinks)
    filter!(!isequal(SequenceGraphLink(dest, source, 0)), dlinks)
    return slinkslen != length(slinks) || dlinkslen != length(dlinks)
end



"""
    find_overlaps(X::Vector{SequenceGraphNode})

Finds the overlaps of length length(Node.sequence)-1 between all Nodes in the vector X

Returns a vector of Tuples where each Tuple (i,j) corresonds to a directed edge (forward link) from i to j

For constructing the De Bruijn Graph links we need the overlapping sequence information
We first formulate the forward links as tuples with indices
Then using these tuple indices we will create the SequenceGraphLinks and finalize the Graph Construction
"""
function find_overlaps(X::Vector{SequenceGraphNode})
    overlaps = Vector{Tuple{Int64,Int64}}()
    for i in 1:size(X)[1]
        for j in 1:size(X)[1]
            if String(X[i].sequence)[2:end]==String(X[j].sequence)[1:end-1]
                #println("Overlap between $(X[i].sequence) and $(X[j].sequence)")
                push!(overlaps,(i,j))
            end
        end
    end
    return overlaps
end


"""
    deBruijn_constructor(kmer_vector::Vector{Kmer{T,K}}) where{T<:NucleicAcid,K}

Returns a BioSequenceGraph constructed by the kmers

For now lets assume that we have Kmers prepared already
So for an unknown DNA sequence we have Spectrum(s,k) where Spectrum(s,l)
is the multiset of n-l+1 l-mers in s.

Our De Bruijn Graph constructor takes as input some kmers and puts a directed link between
every two nodes n1,n2 such that n1[2:] = n2[:k-1], i.e overlap(n1,n2) == k-1
first we find the overlaps between each node by doing an computationally expensive exhaustive search
over all pairs

deBruijn_Constructor takes as input some number of kmers kmer_i, generates Sequence Graph Nodes n_i
such that each n_i.sequence = kmer_i

The constructor generates  nodes in the same order in the input kmer vector

Also the constructor adds empty  vector for each node with no forward link
So the number of vectorsin Links match with number of nodes in Nodes
"""
function deBruijn_constructor(kmer_vector::Vector{Kmer{T,K}}) where{T<:NucleicAcid,K}
    Nodes = Vector{SequenceGraphNode}()
    for kmer in kmer_vector
        node = SequenceGraphNode(kmer,true)
        push!(Nodes,node)
    end
    overlaps = find_overlaps(Nodes)
    Links = Vector{Vector{SequenceGraphLink}}()
    prev_i = 1 ## initial NodeID
    current_node_vector = Vector{SequenceGraphLink}()
    for overlap in overlaps
        if overlap[1]!=prev_i
            push!(Links,current_node_vector)
            current_node_vector = Vector{SequenceGraphLink}()
            for i in prev_i:overlap[1]-1
                push!(Links,Vector{SequenceGraphLink}())
            end
            prev_i = overlap[1]
        end
        ## links are generated from - end of outgoing node to the + end of the incoming one
        link = SequenceGraphLink(-overlap[1],overlap[2],1) ## right now dist is initialized as 1
        push!(current_node_vector,link)
    end
    push!(Links,current_node_vector)
    for i in prev_i+1:size(Nodes)[1]
        push!(Links,Vector{SequenceGraphLink}())
    end
    deBruijn_Graph = SequenceGraph(Nodes,Links)
    deBruijn_Graph
end

"""
    forward_links(sg::SequenceGraph, n::NodeID)

Get the list of links leaving `n`.
"""
function forward_links(sg::SequenceGraph, n::NodeID)
    r = Vector{Link}(0)
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_forward_link(link, n)
            push!(r, link)
        end
    end
    return r
end

backward_links(sg::SequenceGraph, n::NodeID) = forward_links(sg, -n)
