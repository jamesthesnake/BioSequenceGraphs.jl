###
### Old DBG construction from set of Kmers
###

# TODO: We might want to move this function to the BioSequences.jl in the future.
function mer_prefix(mer::DNAKmer{K}) where {K}
    return DNAKmer{K - 1}(UInt64(mer) >> 2)
end

# TODO: We might want to move this function to the BioSequences.jl in the future.
function mer_suffix(mer::DNAKmer{K}) where {K}
    return DNAKmer{K - 1}(UInt64(mer))
end

function new_dbg_from_kmerset(kmers::Set{DNAKmer{K}}) where {K}
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

###
### New unitig graph construction straight from a set of kmers
###

# An improvement on the old algorithm:
# Builds a graph of unitigs straight away.
# Uses a sorted list (vector) of canonical kmers instead of a set. Checking kmer
# membership in the list requires binary lookup, which is slower than the O(1) of
# a set's hash, but typical sets of kmers for bigger genomes use waaay to much
# memory, so a list it is. But, this algorithm tries to query the list as few times
# as possible, and keeps track of all kmers used up in unitigs.
# Also, unlike the previous algorithm, when it came to making connections, an all
# v all double for loop was used, which was a bit pointless, so we avoid that
# in this version, since the two connection vectors are sorted, we can do a single
# pass of the two vectors instead. 
# With decent Kmer counting, this thing should be able to at least put together
# medium size genome like arabidopsis without much trouble. Larger genomes are
# probably fine too, but might need a big machine.
# If we use a large K value that will be allowed by the Skipmer types comming
# in BioSequences v2, the unitig graphs this will produce for medium size genomes
# should be pretty darned good in the first place.

struct Kidx{K}
    kmer::DNAKmer{K}
    idx::UInt64
end

# TODO: Update BioSequences with new neighbour iterators instead.
function kmer_fw_neighbours(mer::DNAKmer{K}) where {K}
    d = BioSequences.encoded_data(mer)
    base = d << 2
    return (DNAKmer{K}(base), DNAKmer{K}(base + 0x01), DNAKmer{K}(base + 0x02), DNAKmer{K}(base + 0x03))
end

function kmer_bw_neighbours(mer::DNAKmer{K}) where {K}
    d = BioSequences.encoded_data(mer)
    base = d >> 2
    BT = typeof(base)
    return (DNAKmer{K}(base), DNAKmer{K}(base + (BT(1) << 2(K - 1))), DNAKmer{K}(base + (BT(2) << 2(K - 1))), DNAKmer{K}(base + (BT(3) << 2(K - 1))))
end

function is_end_bw(mer::DNAKmer{K}, merlist::Vector{DNAKmer{K}}) where {K}
    next = Vector{Kidx{K}}()
    get_bw_idxs!(next, mer, merlist)
    length(next) != 1 && return true
    @inbounds p = next[1].kmer
    get_fw_idxs!(next, p, merlist)
    length(next) != 1 && return true
    return false
end

function is_end_fw(mer::DNAKmer{K}, merlist::Vector{DNAKmer{K}}) where {K}
    next = Vector{Kidx{K}}()
    get_fw_idxs!(next, mer, merlist)
    length(next) != 1 && return true
    @inbounds p = next[1].kmer
    get_bw_idxs!(next, p, merlist)
    length(next) != 1 && return true
    return false
end

function get_bw_idxs!(out::Vector{Kidx{K}}, kmer::DNAKmer{K}, kmerlist::Vector{DNAKmer{K}}) where {K}
    empty!(out)
    for n in kmer_bw_neighbours(kmer)
        cnext = canonical(n)
        cidx = min(searchsortedfirst(kmerlist, cnext), length(kmerlist))
        if @inbounds kmerlist[cidx] == cnext
            push!(out, Kidx{K}(n, cidx))
        end
    end
end

function get_fw_idxs!(out::Vector{Kidx{K}}, kmer::DNAKmer{K}, kmerlist::Vector{DNAKmer{K}}) where {K}
    empty!(out)
    for n in kmer_fw_neighbours(kmer)
        cnext = canonical(n)
        cidx = min(searchsortedfirst(kmerlist, cnext), length(kmerlist))
        if @inbounds kmerlist[cidx] == cnext
            push!(out, Kidx{K}(n, cidx))
        end
    end
end

function new_graph_from_kmerlist(kmerlist::Vector{DNAKmer{K}}) where {K}
    @info string("Start constructing Sequence Graph from ", length(kmerlist), " ", K, "-mers")
    sg = SequenceDistanceGraph{BioSequence{DNAAlphabet{2}}}()
    used_kmers = falses(length(kmerlist))
    bw = Vector{Kidx{K}}()
    fw = Vector{Kidx{K}}()
    
    @info "Building unitigs from kmers"
    for start_kmer_idx in eachindex(kmerlist)
        # Any kmer can only occur in one unitig.
        if used_kmers[start_kmer_idx] 
            continue
        end
        
        # Check if the kmer is an end/junction of a unitig.
        start_kmer = kmerlist[start_kmer_idx]
        end_bw = is_end_bw(start_kmer, kmerlist)
        end_fw = is_end_fw(start_kmer, kmerlist)
        if !end_bw && !end_fw
            continue
        end
        
        if end_bw && end_fw
            # Kmer as unitig
            s = BioSequence{DNAAlphabet{2}}(start_kmer)
            used_kmers[start_kmer_idx] = true
        else
            # A unitig starts on this kmer.
            # Make sure the unitig starts on FW.
            current_kmer = start_kmer
            used_kmers[start_kmer_idx] = true
            
            if end_fw
                current_kmer = reverse_complement(start_kmer)
                end_fw = end_bw
            end
            # Start unitig 
            s = BioSequence{DNAAlphabet{2}}(current_kmer)
            fwn = Vector{Kidx{K}}()
            while !end_fw
                # Add end nucleotide, update current kmer.
                get_fw_idxs!(fwn, current_kmer, kmerlist)
                current_kmer = first(fwn).kmer
                if used_kmers[first(fwn).idx]
                    break # Break circular contigs into lines.
                end
                used_kmers[first(fwn).idx] = true
                s.push!(last(mer))
                end_fw = is_end_fw(current_kmer, kmerlist)
            end
        end
        add_node!(sg, canonical!(s))
    end 
    
    @info "Linking unitigs by their overlaps"
    
    # Save the (k-1)mer in (rev on first k-1 / fw on last k-1) or out ( fw on first k-1 / bw on last k-1)
    in = Vector{Tuple{DNAKmer{K-1},NodeID}}()
    out = Vector{Tuple{DNAKmer{K-1},NodeID}}()
    sizehint!(in, length(nodes(sg)))
    sizehint!(out, length(nodes(sg)))
    for nid in eachindex(nodes(sg))
        nodeseq = nodes(sg)[nid]
        firstmer = DNAKmer{K-1}(nodeseq[1:K - 1])
        if iscanonical(firstmer)
            push!(in, (firstmer, nid))
        else
            push!(out, (reverse_complement(firstmer), nid))
        end
        lastmer = DNAKmer{K-1}(nodeseq[end - (K - 2):end])
        if iscanonical(lastmer)
            push!(out, (lastmer, -nid))
        else
            push!(out, (reverse_complement(lastmer), -nid))
        end
    end
    sort!(in)
    sort!(out)
    
    # Connect all out -> in for all combinations on each kmer.
    next_out_idx = 1
    for i in in
        while first(out[next_out_idx]) < first(i)
            next_out_idx += 1
        end
        
        oidx = next_out_idx
        while oidx <= length(out) && first(out[oidx]) == first(i)
            add_link!(sg, last(i), last(out[oidx]), -K+1) # No support, although we could add the DBG operation as such.
        end
    end
    @info string("Done constructing Sequence Graph from ", length(kmerlist), " ", K, "-mers")
    return sg
end
