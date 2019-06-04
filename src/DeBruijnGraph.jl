"""
    DeBruijnGraph is also SequenceGraph with a special design concept
    k denotes the initial kmer length that is used during constructing the dbg
"""

struct DeBruijnGraph
    nodes::Vector{SequenceGraphNode}
    links::Vector{Vector{SequenceGraphLink}}
    k ::Int
end


nodes(dbg::DeBruijnGraph) = dbg.nodes
node(dbg::DeBruijnGraph, i::NodeID) = nodes(dbg)[abs(i)]
k_value(dbg::DeBruijnGraph) = dbg.k
links(dbg::DeBruijnGraph) = dbg.links
indegree(dbg::DeBruijnGraph,i::NodeID) = count_indegree(dbg,i)
outdegree(dbg::DeBruijnGraph,i::NodeID) = count_outdegree(dbg,i)



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
function add_link!(dbg::DeBruijnGraph,l::SequenceGraphLink)

end


"""
    is_overlap(sn1::SequenceGraphNode,sn2::SequenceGraphNdoe,k::Int)

k is given as an additional information, as the length of the node labels are subject to changes during node merging
Returns true if suffix to prefix overlap of length k-1 exists
Checks overlap of k-1 long suffix of sn1 with prefix of sn2
"""
function is_overlap(sn1::SequenceGraphNode,sn2::SequenceGraphNode,k::Int)
    l1 = length(sn1)
    l2 = length(sn2)
    String(sn1.sequence)[l1-k+2:end]==String(sn2.sequence)[1:k-1]
end

"""
    find_overlaps(X::Vector{SequenceGraphNode})

Finds the overlaps of length length(Node.sequence)-1 between all Nodes in the vector X

Returns a vector of Tuples where each Tuple (i,j) corresonds to a directed edge (forward link) from i to j

For constructing the De Bruijn Graph links we need the overlapping sequence information
We first formulate the forward links as tuples with indices
Then using these tuple indices we will create the SequenceGraphLinks and finalize the Graph Construction
"""
function find_overlaps(X::Vector{SequenceGraphNode},k::Int)
    overlaps = Vector{Tuple{Int64,Int64}}()
    for i in eachindex(X)
        for j in eachindex(X)
            if is_overlap(X[i],X[j],k)
            #if String(X[i].sequence)[2:end]==String(X[j].sequence)[1:end-1]
                #println("AAAAAAAAAAAAA")
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
finds the overlaps between a node and all nodes in a dbg of length k-1
"""
function find_overlaps(dbg::DeBruijnGraph,n::SequenceGraphNode,k::Int)
    nodeid = Base.length(dbg.nodes)+1
    overlaps = Vector{Tuple{Int64,Int64}}()
    for i in 1:nodeid-1
        if is_overlap(dbg.nodes[i],n,k)
        #if String(dbg.nodes[i].sequence)[2:end]==String(n.sequence)[1:end-1]
            push!(overlaps,(i,nodeid))
            println("Overlap between $(dbg.nodes[i].sequence) and $(n.sequence)")
        elseif is_overlap(n,dbg.nodes[i],k)
        #elseif String(dbg.nodes[i].sequence)[1:end-1]==String(n.sequence)[2:end]
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
function add_node!(dbg::DeBruijnGraph,n::SequenceGraphNode)
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
    new_deBruijn_Constructor(kmer_set::Set{Kmer{T,K}})where{T,K}

Input :  a set of kmers in their canonical form

(canonical_kmer_set function under kmer.jl can be used for creating the set out of an array of kmers)

Returns a dbg consisting of a vector of SequenceGraphNodes, a vector of SequenceGraphLinks and an integer k for the kmer length

New Constructor for the DeBruijnGraph where we simultaneously consider a kmer and its reverse complement
We define the beginning of the canonical form of the kmer as the (+) end
And end of the canonical form as the (-) ends

"""
function new_deBruijn_Constructor(kmer_set::Set{Kmer{T,K}})where{T,K}
    fw_nodes = Vector{Tuple{Kmer{T,K-1},Int64}}()
    bw_nodes = Vector{Tuple{Kmer{T,K-1},Int64}}()

    nodes = Vector{SequenceGraphNode}()
    links  = Vector{Vector{SequenceGraphLink}}()

    i = 1
    flag = 0
    ## adding unique kmers to graph in their canonical form
    for kmer in kmer_set
        println(kmer)
        can_kmer = canonical(kmer)
        push!(nodes,SequenceGraphNode(can_kmer,true))
        pref = Kmer{T,K-1}(String(can_kmer)[1:end-1])
        suf = Kmer{T,K-1}(String(can_kmer)[2:end])

        if canonical(pref) == pref ## add prefix to forward nodes
            push!(fw_nodes,(pref,i))
            if pref == reverse_complement(pref)
                push!(bw_nodes,(pref,i))
            end
        else
            push!(bw_nodes,(canonical(pref),i))
        end
        a = 1
        if canonical(suf) == suf ## add suffix to backward nodes (outgoing edge)
            push!(bw_nodes,(suf,-i))
            if suf == reverse_complement(suf)
                push!(fw_nodes,(suf,-i))
            end
        else
            push!(fw_nodes,(canonical(suf),-i))
        end
        i= i + 1
    end
    prev = last(bw_nodes[1])
    links_ = Vector{SequenceGraphLink}()
    for kbn in bw_nodes
        if abs(last(kbn))!=prev
            push!(links,links_)
            for i in prev+1:abs(last(kbn))-1
                push!(links,Vector{SequenceGraphLink}())
            end
            links_ = Vector{SequenceGraphLink}()
        end
        for kfn in fw_nodes
            if first(kfn)==first(kbn)
                print()
                push!(links_,SequenceGraphLink(last(kbn),last(kfn),-K+1)) ## i did not quiet get -k+1
            end
        end
        prev = abs(last(kbn))
    end
    push!(links,links_)
    dbg = DeBruijnGraph(nodes,links,K)
end



"""
    deBruijn_constructor(kmer_vector::Vector{Kmer{T,K}}) where{T<:NucleicAcid,K}

Returns a DeBruijnGraph constructed by the kmers

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
    overlaps = find_overlaps(Nodes,K)
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
    deBruijn_Graph = DeBruijnGraph(Nodes,Links,K)
    deBruijn_Graph
end






# Query Functions
# --------

# Counters
# ------
# Counters count the indegree and outdegree separately for the (+) and  (-) nodes
# edges to and out of the  (+) end  represent the path through canonical kmer
# edges to and out of the  (-) end represent the path through non-canonical kmer

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
    source_ = n  ## checking the sink end of the node for outgoing edges
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

Returns if a given sequence is in a dbg  and also all the internal nodes have in_degree and  out_degree equal to one.


"""

function is_simple_path(dbg::DeBruijnGraph,seq::Sequence)

end




function canonical_bio(x::BioSequence{A})where{A}
    y = reverse_complement(x)
    return x < y ? x : y
end

"""
    is_a_path2(seq::Sequence,dbg::DeBruijnGraph;min_match=3)

Version two for is_a_path for the new dbg designed
Soon will replace version 1

Returns the NodeIDs of all matches and returns true if the path exists completely in dbg

If complete match is not found still returns the indexes of all the matches

min_match is the initial kmer length k used during dbg construction
After node merging we still have k-1 overlaps between nodes because we only merge the simple paths.
Start from a node and traverse its links.
Stop at the first false during traversal as there can not be another path starting from another node!


Making use of is_suffix from Nodes.jl and is_match from sequence.jl
TO-DOs :
    1 - Make sure these dependencies function properly
    2 - Check if a kmer  is contained completely in a node
"""
function is_a_path2(seq::Sequence,dbg::DeBruijnGraph;min_match=3)
    #seq = canonical_bio(seq)
    println(seq)
    index = 0
    indexes = Vector{Int64}()
    nodes_ = nodes(dbg)
    match = -1
    for i in 1:Base.length(nodes_)
        match =  is_suffix2(seq,nodes_[i],direction=1,min_match=min_match)
        match_rev = is_suffix2(seq,nodes_[i],direction=-1,min_match=min_match)
        if match == Base.length(seq) || match_rev == Base.length(seq)
            return true, i
        end

        if match!=-1
            index = -i
            println(nodes_[i])
            break
        elseif match_rev!=-1
            index = i
            match = match_rev
            seq1 = reverse_complement(sequence(nodes_[i]))
            println("Reverse complement Match: "*String(seq1))
            break
        end
    end
    if index==0
        return false

    else ## traverse on the children until remaining is smaller than or equal to min_match length
        ## since kmers overlap with k-1 min_match+1 overlap is carried to the next step
        push!(indexes,index)
        remaining_seq = sub_seq(seq,match-min_match+1)
        rem_len = Base.length(remaining_seq)
        println("Remaining: "*String(remaining_seq))
        println("Remaining length : "* string(Base.length(remaining_seq)))
        println("Index: $(index)")
        #index = index * -1
        while rem_len > 0 ## the remaining must match the next node from index 1
            found = false
            for link in links(dbg)[abs(index)]
                if source(link)==index ## if the direction  is correct
                    dest_ind  = destination(link)
                    node = nodes_[abs(dest_ind)]
                    println("Destination $(dest_ind)")
                    if dest_ind > 0
                        node_seq = sequence(node)
                    else
                        node_seq = reverse_complement(sequence(node))
                    end
                    println("Node sequence : $(node_seq)")
                    min_length = min(rem_len,Base.length(node_seq))
                    if is_match(remaining_seq,1,node_seq,1,min_length)
                        found = true
                        push!(indexes, dest_ind)
                        index = -dest_ind
                        println(nodes_[abs(index)])
                        if min_length ==rem_len
                            println("Match Found!!")
                            return true,indexes
                        end
                        remaining_seq = sub_seq(remaining_seq,min_length-min_match+1)
                        rem_len = Base.length(remaining_seq)
                        println("Remaining : $(remaining_seq)")
                        continue
                    end
                end
            end
            if found ==false
                return false,indexes
            end
        end
    end
    return true,indexes
end

"""
    is_a_path(seq::Sequence,dbg::DeBruijnGraph;min_match=3)

Returns the NodeIDs of all matches and returns true if the path exists completely in dbg

If complete match is not found still returns the indexes of all the matches

min_match is the initial kmer length k used during dbg construction
After node merging we still have k-1 overlaps between nodes because we only merge the simple paths.
Start from a node and traverse its links.
Stop at the first false during traversal as there can not be another path starting from another node!


Making use of is_suffix from Nodes.jl and is_match from sequence.jl
TO-DOs :
    1 - Make sure these dependencies function properly
    2 - Check if a kmer  is contained completely in a node
"""
function is_a_path(seq::Sequence,dbg::DeBruijnGraph;min_match=3)
    index = -1
    indexes = Vector{Int64}()
    nodes_ = nodes(dbg)
    match = -1
    for i in 1:Base.length(nodes_)
        match =  is_suffix(seq,nodes_[i],min_match =min_match)
        if match!=-1
            index = i
            println(nodes_[i])
            break
        end
    end
    if index==-1
        return false
    else ## traverse on the children until remaining is smaller than or equal to min_match length
        ## since kmers overlap with k-1 min_match+1 overlap is carried to the next step
        push!(indexes,index)
        remaining_seq = sub_seq(seq,match+1-min_match+1)
        rem_len = Base.length(remaining_seq)
        println("Remaining: "*String(remaining_seq))
        println("Remaining length : "* string(Base.length(remaining_seq)))
        while rem_len > 0 ## the remaining must match the next node from index 1
            found = false
            for link in links(dbg)[index]
                node = nodes_[destination(link)]
                node_seq = sequence(node)
                min_length = min(rem_len,Base.length(node_seq))
                println(min_length)
                if is_match(remaining_seq,1,node_seq,1,min_length)

                    found = true
                    index = destination(link)
                    push!(indexes,index)
                    println(nodes_[index])
                    if min_length ==rem_len
                        println("Match Found!!")
                        return true,indexes
                    end
                    remaining_seq = sub_seq(remaining_seq,rem_len-min_match+1)
                    rem_len = Base.length(remaining_seq)
                    println(remaining_seq)
                    continue
                end
            end
            if found ==false
                return false,indexes
            end
        end
    end
    return true,indexes
end
