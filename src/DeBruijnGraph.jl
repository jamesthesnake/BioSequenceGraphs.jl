"""
    DeBruijnGraph is also SequenceGraph with a special design concept

"""

struct DeBruijnGraph
    nodes::Vector{SequenceGraphNode}
    links::Vector{Vector{SequenceGraphLink}}
end


nodes(dbg::DeBruijnGraph) = dbg.nodes
node(dbg::DeBruijnGraph, i::NodeID) = nodes(dbg)[abs(i)]
links(dbg::DeBruijnGraph) = dbg.links
indegree(dbg::DeBruijnGraph,i::NodeID) = count_indegree(dbg,i)



"""
    links(sg::SequenceGraph, node::NodeID)

Get all of the links of a Node of a sequence graph.
"""
function links(dbg::DeBruijnGraph, node::NodeID)
    l = links(dbg)
    @inbounds return l[abs(node)]
end

"""

For now this is empty. We have to decide about the design concepts first!!!
"""
function add_link(dbg::DeBruijnGraph,l::SequenceGraphLink)

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
    find_overlaps(dbg::DeBruijnGraph,n::SequenceGraphNode)

Returns a list of tuples which will be used to add links to the DeBruijnGraph
finds the overlaps between a node and all nodes in a dbg
"""
function find_overlaps(dbg::DeBruijnGraph,n::SequenceGraphNode)
    nodeid = Base.length(dbg.nodes)+1
    overlaps = Vector{Tuple{Int64,Int64}}()
    for i in 1:nodeid-1
        if String(dbg.nodes[i].sequence)[2:end]==String(n.sequence)[1:end-1]
            push!(overlaps,(i,nodeid))
            println("Overlap between $(dbg.nodes[i].sequence) and $(n.sequence)")
        elseif String(dbg.nodes[i].sequence)[1:end-1]==String(n.sequence)[2:end]
            push!(overlaps,(nodeid,i))
            println("Overlap between $(dbg.nodes[i].sequence) and $(n.sequence)")
        end
    end
    return overlaps
end


"""
    add_node(dbg::DeBruijnGraph,n::SequenceGraphNode)

Before adding the node to the graph we find the overlaps between the newly added node and then include it
Maybe we can consider checking if it already exists
For now my plan is to add each corresponding link automatically when we add a node
"""
function add_node(dbg::DeBruijnGraph,n::SequenceGraphNode)
    overlaps  = find_overlaps(dbg,n)
    push!(dbg.nodes,n)
    len = Base.length(dbg.nodes)
    push!(dbg.links,Vector{SequenceGraphLink}()) ## push an empty vector to the links
    for overlap in overlaps
        link = SequenceGraphLink(-overlap[1],overlap[2])
        push!(dbg.links[overlap[1]],link)
        println(link)
    end
    nodeid = Base.length(dbg.nodes)
    nodeid
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
            for i in prev_i:overlap[1]-2
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
    deBruijn_Graph = DeBruijnGraph(Nodes,Links)
    deBruijn_Graph
end


# Query Functions
# --------

# Counters
# ------

"""
    count_indegree(dbg::DeBruijnGraph,n::NodeID)

Returns the number of incoming edges to the SequenceGraphNode with id NodeID on the DeBruijnGraph

Makes an exhaustive search on all the SequenceGraphLinks on the dbg
Now checking only the source end of the node we have to discuss about the design
Maybe we can store incoming and outgoing edges for each node for O(N) query time where N denotes number of SequenceGraphNodes
"""
function count_indegree(dbg::DeBruijnGraph,n::NodeID)
    dest = n
    in_degree = 0
    for i in 1:Base.length(nodes(dbg))
        links_ = links(dbg,i)
        for l in links_
            if destination(l)==dest  ## checking only the source end of the node, not sure!
                in_degree +=1
            end
        end
    end
    in_degree
end



"""
    count_outdegree(dbg::DeBruijnGraph,n::NodeID)

Returns the number of incoming edges to the SequenceGraphNode with id NodeID on the DeBruijnGraph

Makes an exhaustive search on all the SequenceGraphLinks on the dbg
Now checking only the sink end of the node we have to discuss about the design
O(K) query time
"""
function count_outdegree(dbg::DeBruijnGraph,n::NodeID)
    source_ = -n  ## checking the sink end of the node for outgoing edges
    out_degree = 0
    for link in links(dbg,n)
        if source_ == source(link)
            out_degree +=1
        end
    end
    out_degree
end


# Path queries
# ------

"""

    For now we assume that only one vertex exist with associated with a certain kmer
"""
function is_a_path(dbg::DeBruijnGraph,seq::BioSequence)
    


end
