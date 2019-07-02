
# Node and Link types for SequenceDistanceGraph
# ---------------------------------------------

const NodeID = Int64

"""
The SDGNode type represents a node in a SequenceDistanceGraph.

At present it contains only two fields, first it holds an instance of a BioSequences
Sequence type.

Secondly, it tracks a flag which indicates if the node has been
deleted or not.

!!! note
    The deleted flag allows us to mark nodes in the graph as deleted, which can be of
    help in some algorithms where the graph structure is being edited (merging nodes for example).

    Actually deleting the node would shift node IDs and require redoing all links in the graph
    and so on.

    So just marking a node as deleted and not using it anymore is a lazy but sometimes
    helpful choice.
"""
struct SDGNode{S}
    seq::S
    deleted::Bool
end

function empty_node(::Type{S}) where {S <: Sequence}
    return SDGNode{S}(empty_seq(S), true)
end

# TODO: This is a hacked copy of the dna string literal macro from BioSequences,
# except it creates 2-bit based DNA sequences rather than 4 bit based ones.
# This ability to choose the bit encoding should make its way to BioSequences.jl
# in the future, but for now, it's here.
#
# I basically want this as it lets me create a single literal empty sequence, shared
# by all deleted SDG nodes. Rather than having each deleted SDG node create a new empty
# sequence.
macro dna2_str(seq, flag)
    if flag == "s"
        return BioSequence{DNAAlphabet{2}}(BioSequences.remove_newlines(seq))
    elseif flag == "d"
        return quote
            BioSequence{DNAAlphabet{2}}($(BioSequences.remove_newlines(seq)))
        end
    end
    error("Invalid DNA flag: '$(flag)'")
end

empty_seq(::Type{BioSequence{DNAAlphabet{2}}}) = dna2""s

"""
Represents a single distance between two sequences in a SequenceDistanceGraph.
"""
struct DistanceGraphLink
    source::NodeID
    destination::NodeID
    dist::Int64
end

source(l::DistanceGraphLink) = l.source
destination(l::DistanceGraphLink) = l.destination
distance(l::DistanceGraphLink) = l.dist

"Test if link `l` is a forward link leaving node `n`."
is_forward_from(l::DistanceGraphLink, n::NodeID) = source(l) == -n
is_backward_from(l::DistanceGraphLink, n::NodeID) = source(l) == n

"""
The SequenceDistanceGraph is a representation of a genome assembly.
Sequences are contained in nodes, and the distances are represented by links.

A singe node represents a sequence *and* its reverse complement.
Every node has a correlative ID starting from 1.
For every node X in the graph, the negative ID -X is mapped to the reverse
complement of X. This mapping is virtual: Only one node is stored in the graph.
This is because every node has an orientaton: Each node has a positive end (+),
and a negative end (-).
So when a node is accessed with (or traversed by entering) the positive end the
node yields the stored sequence.
Conversely, when a node is accessed with (or traversed by entering) the negative
end the node yelds the reverse complement of the stored sequence.
In this way the positive end can be thought of as the sequence start, and the
negative end can be thought of as the sequence end.

A single distance between two sequences is represented as a single link.
Every link connects two node ends and contains a distance (they take the form
`([+, -]n1, [+, -]n2, [+, -]dist)`).
A link connects two node ends, and so the order of the signed nodes in the links
does not change the link.
If the distance in a link is negative, this represents an overlap between two
sequences. These overlaps must be "perfect overlaps".
"""
struct SequenceDistanceGraph{S<:Sequence}
    nodes::Vector{SDGNode{S}}
    links::Vector{Vector{DistanceGraphLink}}
end

function SequenceDistanceGraph{S}() where {S<:Sequence}
    return SequenceDistanceGraph{S}(Vector{SDGNode{S}}(), Vector{Vector{DistanceGraphLink}}())
end

# Graph accessor functions
# ------------------------

nodes(sg::SequenceDistanceGraph) = sg.nodes
node(sg::SequenceDistanceGraph, i::NodeID) = nodes(sg)[abs(i)]
links(sg::SequenceDistanceGraph) = sg.links



"""
    links(sg::SequenceGraph, node::NodeID)

Get all of the links of a Node of a sequence graph.

!!! note
    links accepts a NodeID that can be positive or negative.
    E.g. providing either 5 or -5 both mean node 5 in a graph,
    and so you will get the links for node 5.
"""
function links(sg::SequenceDistanceGraph, node::NodeID)
    l = links(sg)
    return l[abs(node)]
end

# Graph editing operations
# ------------------------

"""
    add_node!(sg::SequenceDistanceGraph{S}, n::SDGNode{S}) where {S<:Sequence}

Add a node to a sequence distance graph.

Returns the node ID used to access the new node added in the graph.

!!! warning
    Currently, we don't enforce the sequence in the node is canonical here.
    We just trust that it is canonical.
"""
function add_node!(sg::SequenceDistanceGraph{S}, n::SDGNode{S}) where {S<:Sequence}
    newlen = length(push!(nodes(sg),n))
    push!(links(sg), Vector{DistanceGraphLink}())
    return newlen
end

"""
   add_node!(sg::SequenceDistanceGraph{S}, seq::Sequence) where {S<:Sequence}

Add a sequence to a sequence distance graph as a node.

Returns the node ID used to access the new node added in the graph.

Can accept any sequence type and will attempt to coerce the input sequence to the
type required by the graph.
"""
function add_node!(sg::SequenceDistanceGraph{S}, seq::Sequence) where {S<:Sequence}
    return add_node!(sg, SDGNode{S}(convert(S, seq), false))
end

"""
    remove_node!(sg::SequenceDistanceGraph{S}, n::NodeID) where {S<:Sequence}

Remove a node from a sequence distance graph.
"""
function remove_node!(sg::SequenceDistanceGraph{S}, n::NodeID) where {S<:Sequence}
    oldlinks = copy(links(sg, n))
    for oldlink in oldlinks
        remove_link!(sg, source(oldlink), destination(oldlink))
    end
    # TODO: This is a lazy solution to getting rid of the node.
    nodes[n] = empty_node(S)
end


"""
    add_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID, dist::Int)

Construct a link between two nodes in a sequence Graph.
"""
function add_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID, dist::Int)
    # Guard against someone adding links using ID's bigger than the current max NodeID in the graph.
    if abs(source) > length(links(sg))
        resize!(links(sg), abs(source))
    end
    if abs(dest) > length(links(sg))
        resize!(links(sg), abs(dest))
    end
    push!(links(sg, source), DistanceGraphLink(source, dest, dist))
    push!(links(sg, dest), DistanceGraphLink(dest, source, dist))
end

"""
    remove_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID)

Remove a link between two nodes in a SequenceDistanceGraph.
Returns a boolean indicating whether the removal was successful.
Reasons this function would not return `true` include that the link
didn't exist in the graph, and so could not be removed.
"""
function remove_link!(sg::SequenceDistanceGraph, source::NodeID, dest::NodeID)
    slinks = links(sg, source)
    slinkslen = length(slinks)
    filter!(!isequal(SequenceGraphLink(source, dest, 0)), slinks)
    dlinks = links(sg, dest)
    dlinkslen = length(dlinks)
    filter!(!isequal(SequenceGraphLink(dest, source, 0)), dlinks)
    return slinkslen != length(slinks) || dlinkslen != length(dlinks)
end

remove_link!(sg::SequenceDistanceGraph, lnk::DistanceGraphLink) = remove_link!(sg, source(lnk), destination(lnk))

"""
Removes all the links in the collection from and to a given nodeID.
"""
function disconnect_node!(sg::SequenceDistanceGraph, n::NodeID)
    for flink in forward_links(sg, n)
        remove_link!(sg, flink)
    end
    for rlink in backward_links(sg, n)
        remove_link!(sg, rlink)
    end
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


# Graph traversing operations
# ---------------------------

"""
    forward_links(sg::SequenceDistanceGraph, n::NodeID)

Get a vector of the links leaving `n` forward from that node.
"""
function forward_links(sg::SequenceDistanceGraph, n::NodeID)
    r = Vector{Link}(0)
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_forward_from(link, n)
            push!(r, link)
        end
    end
    return r
end

"""
    backward_links(sg::SequenceDistanceGraph, n::NodeID)

Get a vector of the links leaving `n` backwards from that node.
"""
backward_links(sg::SequenceDistanceGraph, n::NodeID) = forward_links(sg, -n)


"""
    get_next_nodes(sg::SequenceDistanceGraph, n::NodeID)

Find node IDs for forward nodes of `n`.
"""
function get_next_nodes(sg::SequenceDistanceGraph, n::NodeID)
    r = Vector{NodeID}(0)
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_forward_from(link, n)
            push!(r, destination(link))
        end
    end
    return r
end

"""
    get_previous_nodes(sg::SequenceDistanceGraph, n::NodeID)

Find node IDs for backward nodes of `n`.
"""
function get_previous_nodes(sg::SequenceDistanceGraph, n::NodeID)
    r = Vector{NodeID}(0)
    nodelinks = links(sg, n)
    sizehint!(r, length(nodelinks))
    for link in nodelinks
        if is_backwards_from(link, n)
            push!(r, -destination(link))
        end
    end
end


# TODO: We might want to move this function to the BioSequences.jl in the future.
function mer_prefix(mer::DNAKmer{K}) where {K}
    return DNAKmer{K - 1}(UInt64(mer) >> 2)
end

# TODO: We might want to move this function to the BioSequences.jl in the future.
function mer_suffix(mer::DNAKmer{K}) where {K}
    return DNAKmer{K - 1}(UInt64(mer))
end

function build_dbg(kmers::Set{DNAKmer{K}}) where {K}
    canonical_kmers = Set{DNAKmer{K}}(canonical(mer) for mer in kmers)

    kmer_ovl_bw_nodes = Vector{Tuple{DNAKmer{K-1}, Int64}}()
    kmer_ovl_fw_nodes = Vector{Tuple{DNAKmer{K-1}, Int64}}()

    sg = SequenceDistanceGraph{DNAKmer{K}}()

    for mer in canonical_kmers
        # Add the canonical kmer to the nodes of the graph, and note it's ID.
        nodeid = add_node!(sg, mer)
        # Take the prefix (k-1) of the canonical kmer.
        pre = mer_prefix(mer)
        # If the prefix is canonical, push the tuple (canonical(prefix), +nodeid) to `kmer_ovl_fw_nodes`.
        # Else push the tuple(canonical(prefix), +nodeid) to the `kmer_ovl_bw_nodes`.
        prerc = reverse_complement(pre)
        if pre < prerc
            push!(kmer_ovl_fw_nodes, (pre, nodeid))
        else
            push!(kmer_ovl_bw_nodes, (prerc, nodeid))
        end
        
        # Take the suffix (k-1) of the canonical kmer.
        suf = mer_suffix(mer)
        # If the suffix is canonical, push the tuple (canonical(suffix), -nodeid) to `kmer_ovl_bw_nodes`.
        # Else push the tuple (canonical(suffix), - nodeid) to `kmer_ovl_fw_nodes`.
        sufrc = reverse_complement(suf)
        if suf < sufrc
            push!(kmer_ovl_bw_nodes, (suf, -nodeid))
        else
            push!(kmer_ovl_fw_nodes, (sufrc, -nodeid))
        end
    end
    
    # Sort `kmer_ovl_fw_nodes` & `kmer_ovl_bw_nodes`.
    sort!(kmer_ovl_fw_nodes)
    sort!(kmer_ovl_bw_nodes)
    
    for kbn in kmer_ovl_bw_nodes
        for kfn in kmer_ovl_fw_nodes
            if first(kbn) == first(kfn)
                # Add a link to graph with source = last(kbn), destination = last(kfn) and distance = -k + 1
                add_link!(sg, last(kbn), last(kfn), -K + 1)
            end
        end
    end
    
    return sg
end
