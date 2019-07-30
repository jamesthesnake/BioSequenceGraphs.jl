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

encoded_data(mer) = reinterpret(UInt64, mer)

function iscanonical(seq)
    i = 1
    j = length(seq)
    @inbounds while i <= j
        f = seq[i]
        r = complement(seq[j])
        f < r && return true
        r < f && return false
        i += 1
        j -= 1
    end
    return true
end

function canonical!(seq::BioSequence{DNAAlphabet{2}})
    if !iscanonical(seq)
        reverse_complement!(seq)
    end
    return seq
end

# TODO: Update BioSequences with new neighbour iterators instead. For now, these
# functions will do.
function kmer_fw_neighbours(mer::DNAKmer{K}) where {K}
    d = encoded_data(mer)
    base = d << 2
    return (DNAKmer{K}(base), DNAKmer{K}(base + 0x01), DNAKmer{K}(base + 0x02), DNAKmer{K}(base + 0x03))
end

function kmer_bw_neighbours(mer::DNAKmer{K}) where {K}
    d = encoded_data(mer)
    base = d >> 2
    BT = typeof(base)
    return (DNAKmer{K}(base), DNAKmer{K}(base + (BT(1) << 2(K - 1))), DNAKmer{K}(base + (BT(2) << 2(K - 1))), DNAKmer{K}(base + (BT(3) << 2(K - 1))))
end

function is_end_bw(mer::DNAKmer{K}, merlist::Vector{DNAKmer{K}}) where {K}
    @debug "Checking if kmer is end BW" mer
    next = Vector{Kidx{K}}()
    get_bw_idxs!(next, mer, merlist)
    @debug "BW neighbours:" next
    length(next) != 1 && return true
    @inbounds p = next[1].kmer
    get_fw_idxs!(next, p, merlist)
    @debug "FW neighbours of only BW neighbour:" p next
    length(next) != 1 && return true
    return false
end

function is_end_fw(mer::DNAKmer{K}, merlist::Vector{DNAKmer{K}}) where {K}
    @debug "Checking if kmer is end FW" mer
    next = Vector{Kidx{K}}()
    get_fw_idxs!(next, mer, merlist)
    @debug "FW neighbours:" next
    length(next) != 1 && return true
    @inbounds p = next[1].kmer
    get_bw_idxs!(next, p, merlist)
    @debug "BW neighbours of only FW neighbour:" p next
    length(next) != 1 && return true
    return false
end

function get_bw_idxs2!(out::Vector{Kidx{K}}, kmer::DNAKmer{K}, kmerlist::Vector{DNAKmer{K}}) where {K}
    empty!(out)
    direction_flag = 1
    for n in kmer_bw_neighbours(kmer)
        cnext = canonical(n)
        if n!=cnext
            direction_flag = -1
        end
        cidx = min(searchsortedfirst(kmerlist, cnext), length(kmerlist))
        if @inbounds kmerlist[cidx] == cnext
            push!(out, Kidx{K}(n, cidx*direction_flag))
        end
    end
end

function get_fw_idxs2!(out::Vector{Kidx{K}}, kmer::DNAKmer{K}, kmerlist::Vector{DNAKmer{K}}) where {K}
    empty!(out)
    direction_flag = 1
    for n in kmer_fw_neighbours(kmer)
        cnext = canonical(n)
        if n!=cnext
            direction_flag = -1 ## This is used to traverse the unitig properly
        end
        cidx = min(searchsortedfirst(kmerlist, cnext), length(kmerlist))
        if @inbounds kmerlist[cidx] == cnext
            push!(out, Kidx{K}(n, cidx*direction_flag))
        end
    end
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

const GRAPH_TYPE = SequenceDistanceGraph{BioSequence{DNAAlphabet{2}}}

## Bubble removal
## Here we would like to remove one branch of a Bubble in a list of kmerlist
## Assuming that there is a single nucleotide error in a certain read
## It will generate a branch of length k which starts and ends at the kmer with the correct gene sequence
## To detect which branch corresponds to an error a common approach is to coverage for each branch
## For now we assume that we have the coverage information for each kmer and we use this to delete the kmers on the branch with low coverage
"""
    pop_bubbles!(all_paths)

    !WARNING!
    !!!all paths in all_paths must be in same direction!! otherwise the gluing operation does not work correctly!!!

    Removes the low coverage paths in all bubbles
    all_paths contain list of triples : the node_id indices, nucleotide sequence and coverage for each path
    kmer_num : number of kmers inputted for dbg creation (used for indexing)
"""
function pop_bubbles(all_paths,kmer_num)
    #assert(Base.length(all_paths)==Base.length(kmercoverage))
    print("All paths are as below : ")
    println(all_paths)
    bubbles = Vector{Tuple{Int64,Int64}}()
    for i in eachindex(all_paths)
        for j in i+1:Base.length(all_paths)
            if all_paths[i][1][1]==all_paths[j][1][1] && all_paths[i][1][end]==all_paths[j][1][end]## bubble found
                @info string("Bubble detected for unitigs : ", all_paths[i] ," and ",all_paths[j])
                push!(bubbles,(i,j))
            end
        end
    end
    for bubble in bubbles
        p1 = bubble[1]
        p2 = bubble[2]
        if all_paths[p1][3]>all_paths[p2][3]
            @info string("Deleting  unitig: ", all_paths[p1])
            deleteat!(all_paths,p1)
        else
            @info string("Deleting  unitig: ", all_paths[p2])
            deleteat!(all_paths,p2)
        end
    end

    ## after deleting the bubbles we may have newly formed longer new_contigs
    ## at this point we check end and start points of each contig and combine if they form a contig together
    sort!(all_paths)
    start_counts = zeros(kmer_num)
    end_counts = zeros(kmer_num)
    cum_counts = zeros(kmer_num)
    #cum_end_counts = zeros(Base.length(all_paths))
    for path in all_paths
        start_ind = path[1][1]
        end_ind   = path[1][end]
        start_counts[start_ind]+=1
        end_counts[end_ind]+=1
        for i in start_ind:kmer_num
            cum_counts[i]+=1
        end
    end
    merged_paths = Vector{Tuple{Int64,Int64}}()
    used_contigs = zeros(Base.length(all_paths))
    new_contigs = Vector{Vector{Int64}}()
    ## a bit inefficient for now
    ## we can only look at indices that satisfy these conditions as we have checked them already
    for i in eachindex(all_paths)
        end_ind = all_paths[i][1][end]
        if used_contigs[i] == 1
            continue
        end
        if start_counts[end_ind]==1 && end_counts[end_ind]==1
            next_contig_ind = cum_counts[end_ind]
            used_contigs[next_contig_ind ] = 1
            push!(merged_paths,(i,next_contig_ind))
            push!(new_contigs,vcat(all_paths[i][1][1:end-1],all_paths[next_contig_ind][1]))
        else
            push!(new_contigs,all_paths[i][1])
        end
        used_contigs[i]=1
    end
    ## I think there is a problem with the start index node
    """
    for pair in merged_paths
        ind1= pair[1]
        ind2= pair[2]
        push(new_contigs,vcat(all_paths[ind1][1][1:end-1],all_paths[ind2][1]))
    end
    """
    return new_contigs
end

## Tip removal

"""

    function delete_tips(kmerlist::Vector{DNAKmer{K}}) where {K}

Now we assume that the shortest tip is always the one to be removed!

"""

function delete_tips(kmerlist::Vector{DNAKmer{K}}) where {K}
    sort!(kmerlist)
    @info string("Deleting short tips for the kmerlist :  ", kmerlist)
    used= falses(length(kmerlist))
    all_tips = Dict{Int64,Vector{Vector{Int64}}}()## only store the shorest tip from each parent kmer
    for kmer_ind in  eachindex(kmerlist)
        if used[kmer_ind]
            continue
        end
        path = Vector{Int64}()
        mer = kmerlist[kmer_ind]
        next = Vector{Kidx{K}}()
        get_fw_idxs!(next, mer, kmerlist)
        next2 = Vector{Kidx{K}}() # not used
        get_bw_idxs!(next2, mer, kmerlist)
        if Base.length(next)==0 && Base.length(next2)==1 || Base.length(next)==1 && Base.length(next2)==0## Start of a tip from the current mer
            push!(path,kmer_ind)
            prev_mer = mer
            prev_ind = kmer_ind
            next = vcat(next,next2)
            next_ind = next[1].idx ## only kmer neighbor
            next_mer = kmerlist[next_ind]
            get_fw_idxs!(next, next_mer, kmerlist)
            get_bw_idxs!(next2, next_mer, kmerlist)
            while Base.length(next)==1 && Base.length(next2)==1 ## At all the kmers on the simple path
                next = vcat(next,next2)
                push!(path,next_ind)
                temp_mer = next_mer
                temp_ind = next_ind
                next_mer  = next[1].kmer==prev_mer ? next[2].kmer  : next[1].kmer ## get the OTHER neighbor of the current kmer
                next_ind =  next[1].kmer==prev_mer ? next[2].idx :  next[1].idx
                prev_mer,prev_ind = temp_mer,temp_ind
                get_fw_idxs!(next, next_mer, kmerlist)
                get_bw_idxs!(next2, next_mer, kmerlist)
            end
            ## We know that the kmer that terminated the while loop has either multiple backward or multiple forward links
            ## How about the zero case? anyway
            if next_ind in keys(all_tips)
                push!(all_tips[next_ind],path)
            else
                all_tips[next_ind] = [path]
            end
            @info string("Found a tip ",  path , " starting from kmer : " ,kmerlist[next_ind])
        end
    end
    ## only delete shortest tips that are branching (otherwise removes all final unitigs!!!)
    @info string("All possible tips found : ",all_tips)
    dead_ends = Vector{Vector{Int64}}()
    for key in keys(all_tips)
        if Base.length(all_tips[key])>1
            shortest_ind = argmin(map(Base.length,all_tips[key]))
            push!(dead_ends,all_tips[key][shortest_ind])
        end
    end
    @info string("Shortest tips to be removed : ",dead_ends)
    deads = Vector{Int64}()
    for end1 in dead_ends
        for d in end1
            push!(deads,d)
        end
    end
    new_kmer_list = Vector{DNAKmer{K}}()
    for kmer_ind in eachindex(kmerlist)
        if !(kmer_ind in deads)
            push!(new_kmer_list,kmerlist[kmer_ind])
        end
    end
    @info string("New kmer list ", new_kmer_list)
end


## right now does not know how to get coverage so for now assume it is the average of all kmers in the list
function get_coverage(path::Vector{Int64},coverage::Vector{Int64})
    norm_cover = 0
    for x in path
        norm_cover +=coverage[x]
    end
    norm_cover/Base.length(path)
end

function generate_coverage(v::Int64)
    rand(1:100,v)
end

function get_sequence(node_list,kmerlist)
    s = BioSequence{DNAAlphabet{2}}(kmerlist[node_list[1]])
    for node_ind in eachindex(node_list[2:end])
        push!(s,last(kmerlist[node_list[node_ind+1]]))
    end
    s
end
"""
    build_unitigs_from_kmerlist2

This implementation of unitig building takes as input kmer coverage in addition to a kmerlist
And remove the low-coverage branches (bubble popping) before producing the final unitig list

"""
function build_unitigs_from_kmerlist2!(sg::GRAPH_TYPE, kmerlist::Vector{DNAKmer{K}},kmercounts::Vector{Int64}) where {K}
    @info string("Constructing unitigs from ", length(kmerlist), " ", K, "-mers")
    used_kmers = falses(length(kmerlist))
    single_kmers = Vector{Int64}()
    ## each path is a triple (node indices,nucleotide sequence , coverage)
    all_paths = Vector{Tuple{Vector{Int64},BioSequence{DNAAlphabet{2}},Float64}}()
    for start_kmer_idx in eachindex(kmerlist)
        @debug "Considering new kmer" start_kmer_idx
        @info string("Considering new kmer: ", start_kmer_idx)
        # Any kmer can only occur in one unitig.
        if used_kmers[start_kmer_idx]
            @info string("Kmer has been used ", start_kmer_idx)
            continue
        end

        # Check if the kmer is an end/junction of a unitig.
        start_kmer = kmerlist[start_kmer_idx]
        end_bw = is_end_bw(start_kmer, kmerlist)
        end_fw = is_end_fw(start_kmer, kmerlist)
        @info string(start_kmer_idx, " icin on arka ", end_bw, " " , end_fw)
        if !end_bw && !end_fw
            @info string("Kmer is middle of a unitig " ,start_kmer_idx , " ", start_kmer)
            continue
        end

        if end_bw && end_fw
            @info string("Kmer is single unitig " ,start_kmer_idx , "  ",  start_kmer)
            # Kmer as unitig
            push!(single_kmers,start_kmer_idx)
            """
            next = Vector{Kidx{K}}()
            next2 = Vector{Kidx{K}}()
            get_fw_idxs!(next, start_kmer, kmerlist)## we find the previous node for bubble popping
            get_bw_idxs!(next2, start_kmer, kmerlist)## we find the previous node for bubble popping
            s = BioSequence{DNAAlphabet{2}}(start_kmer)
            if Base.length(next)==1
                path_nodes= Vector{Int64}([next[1].idx,start_kmer_idx])
                coverage = get_coverage(path_nodes,kmercounts)
                push!(all_paths,(path_nodes,s,coverage))
            elseif Base.length(next2)==1
                path_nodes= Vector{Int64}([next2[1].idx,start_kmer_idx])
                coverage = get_coverage(path_nodes,kmercounts)
                push!(all_paths,(path_nodes,s,coverage)
            """
            ## maybe add these later as separate contigs to the new contigs list
            used_kmers[start_kmer_idx] = true

        else
            # A unitig starts on this kmer.
            # Make sure the unitig starts on FW.
            @info string("Constructing unitigs starting from ", start_kmer)
            current_kmer = start_kmer
            current_kmer_idx = start_kmer_idx
            used_kmers[start_kmer_idx] = true
            if end_fw
                current_kmer = reverse_complement(start_kmer)
                end_fw = end_bw
            end
            @debug "Start of unitig" start_kmer current_kmer end_bw end_fw
            # Start unitig
            next = Vector{Kidx{K}}()
            get_bw_idxs!(next, current_kmer, kmerlist)## we find the previous node for bubble popping
            s = BioSequence{DNAAlphabet{2}}(current_kmer)
            final_node = -1
            fwn = Vector{Kidx{K}}()
            if Base.length(next)!=0
                start_node = next[1]
                path_nodes = Vector{Int64}()
                push!(path_nodes,start_node.idx)
            end
            while !end_fw
                # Add end nucleotide, update current kmer.
                get_fw_idxs!(fwn, current_kmer, kmerlist)
                @debug "Extending unitig" fwn
                push!(path_nodes,current_kmer_idx)
                current_kmer = first(fwn).kmer
                current_kmer_idx = first(fwn).idx
                if used_kmers[first(fwn).idx]
                    @debug "New kmer is already used" current_kmer
                    break # Break circular contigs into lines.
                end
                used_kmers[first(fwn).idx] = true
                push!(s, last(current_kmer))
                end_fw = is_end_fw(current_kmer, kmerlist)
            end
            push!(path_nodes,current_kmer_idx)## end node is also stored
            get_fw_idxs!(fwn, current_kmer, kmerlist)
            if Base.length(fwn)==1 ## add the final node as well for 2-1 type nodes for finding consecutive contigs
                push!(path_nodes,fwn[1].idx)
            end
            coverage = get_coverage(path_nodes,kmercounts)
            @info string("Found a contig starting from  ", path_nodes[1]," and ending at ", path_nodes[end]," with coverage " , coverage)

            push!(all_paths,(path_nodes,s,coverage))
        end
        #add_node!(sg, canonical!(s)) lets not do it now
        ## here we have a path corresponding almost to sequence s we can combine

    end

    new_contigs = pop_bubbles(all_paths,Base.length(kmerlist))
    println("New Paths after popping:")
    println(new_contigs)
    println("Initial single kmers")
    println(single_kmers)
    for contig in new_contigs
        if contig[1] in single_kmers
            @info string("Skipping ", kmerlist[contig[1]] , " as it is contained in the single kmers list")
            s = get_sequence(contig[2:end],kmerlist)
            add_node!(sg,canonical!(s))
            println(s)
        else
            s = get_sequence(contig[2:end],kmerlist)
            add_node!(sg,canonical!(s))
            println(s)
        end
    end

    for single_kmer in single_kmers
        println(kmerlist[single_kmer])
        add_node!(sg,canonical!(BioSequence{DNAAlphabet{2}}(kmerlist[single_kmer])))
    end
    # A temporary check for circle problem for now.
    if !all(used_kmers)
        @warn "Some kmers have not been incorporated into unitigs. This may be a case of the circle problem" kmerlist[(!).(used_kmers)]
    end
    @info string("Constructed ", length(nodes(sg)), " unitigs")
    return sg
end


function build_unitigs_from_kmerlist!(sg::GRAPH_TYPE, kmerlist::Vector{DNAKmer{K}}) where {K}
    @info string("Constructing unitigs from ", length(kmerlist), " ", K, "-mers")
    used_kmers = falses(length(kmerlist))

    for start_kmer_idx in eachindex(kmerlist)
        @debug "Considering new kmer" start_kmer_idx

        # Any kmer can only occur in one unitig.
        if used_kmers[start_kmer_idx]
            @debug "Kmer has been used" start_kmer_idx
            continue
        end

        # Check if the kmer is an end/junction of a unitig.
        start_kmer = kmerlist[start_kmer_idx]
        end_bw = is_end_bw(start_kmer, kmerlist)
        end_fw = is_end_fw(start_kmer, kmerlist)

        if !end_bw && !end_fw
            @debug "Kmer is middle of a unitig" start_kmer_idx start_kmer
            continue
        end

        if end_bw && end_fw
            @debug "Kmer is single unitig" start_kmer_idx start_kmer
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
            @debug "Start of unitig" start_kmer current_kmer end_bw end_fw
            # Start unitig
            s = BioSequence{DNAAlphabet{2}}(current_kmer)
            fwn = Vector{Kidx{K}}()
            while !end_fw
                # Add end nucleotide, update current kmer.
                get_fw_idxs!(fwn, current_kmer, kmerlist)
                @debug "Extending unitig" fwn
                current_kmer = first(fwn).kmer
                if used_kmers[first(fwn).idx]
                    @debug "New kmer is already used" current_kmer
                    break # Break circular contigs into lines.
                end
                used_kmers[first(fwn).idx] = true
                push!(s, last(current_kmer))
                end_fw = is_end_fw(current_kmer, kmerlist)
            end
        end
        add_node!(sg, canonical!(s))
    end
    # A temporary check for circle problem for now.
    if !all(used_kmers)
        @warn "Some kmers have not been incorporated into unitigs. This may be a case of the circle problem" kmerlist[(!).(used_kmers)]
    end
    @info string("Constructed ", length(nodes(sg)), " unitigs")
    return sg
end

function find_unitig_overlaps(sg::GRAPH_TYPE, ::Type{DNAKmer{K}}) where {K}
    @info string("Identifying the ", K - 1, "bp (K - 1) overlaps between ", length(nodes(sg)), " unitigs")
    # Save the (k-1)mer in (rev on first k-1 / fw on last k-1) or out ( fw on first k-1 / bw on last k-1)
    @debug "Sorting K - 1 overlaps as `in` or `out`"
    in = Vector{Tuple{DNAKmer{K-1},NodeID}}()
    out = Vector{Tuple{DNAKmer{K-1},NodeID}}()
    sizehint!(in, length(nodes(sg)))
    sizehint!(out, length(nodes(sg)))
    for nid in eachindex(nodes(sg))
        nodeseq = node(sg, nid).seq
        firstmer = DNAKmer{K-1}(nodeseq[1:K - 1])
        @debug string("Considering node ", nid) nodeseq
        if iscanonical(firstmer)
            @debug "Source overlap is canonical"
            push!(in, (firstmer, nid))
        else
            @debug "Source overlap is not canonical"
            push!(out, (reverse_complement(firstmer), nid))
        end
        lastmer = DNAKmer{K-1}(nodeseq[end - (K - 2):end])
        if iscanonical(lastmer)
            @debug "Sink overlap is canonical"
            push!(out, (lastmer, -nid))
        else
            @debug "Sink overlap is not canonical"
            push!(in, (reverse_complement(lastmer), -nid))
        end
    end
    sort!(in)
    sort!(out)
    return in, out
end

function connect_unitigs_by_overlaps!(sg::GRAPH_TYPE, ::Type{DNAKmer{K}}) where {K}
    in, out = find_unitig_overlaps(sg, DNAKmer{K})
    ol = length(out)
    @info string("Linking ", length(nodes(sg)), " unitigs by their ", K - 1, "bp (K - 1) overlaps")
    # Connect all out -> in for all combinations on each kmer.
    next_out_idx = 1
    for i in in
        while next_out_idx <= ol && first(out[next_out_idx]) < first(i)
            next_out_idx += 1
        end
        oidx = next_out_idx
        while oidx <= ol && first(out[oidx]) == first(i)
            add_link!(sg, last(i), last(out[oidx]), -K + 1) # No support, although we could add the DBG operation as such.
            oidx += 1
        end
    end
end
function new_graph_from_kmerlist2(kmerlist::Vector{DNAKmer{K}},kmer_counts::Vector{Int64}) where {K}
    str = string("onstructing Sequence Distance Graph from ", length(kmerlist), ' ', K, "-mers")
    @info string('C', str)
    sg = GRAPH_TYPE()
    build_unitigs_from_kmerlist2!(sg, kmerlist,kmer_counts)
    if n_nodes(sg) > 1
        connect_unitigs_by_overlaps!(sg, DNAKmer{K})
    end
    @info string("Done c", str)
    return sg
end

function new_graph_from_kmerlist(kmerlist::Vector{DNAKmer{K}}) where {K}
    str = string("onstructing Sequence Distance Graph from ", length(kmerlist), ' ', K, "-mers")
    @info string('C', str)
    sg = GRAPH_TYPE()
    build_unitigs_from_kmerlist!(sg, kmerlist)
    if n_nodes(sg) > 1
        connect_unitigs_by_overlaps!(sg, DNAKmer{K})
    end
    @info string("Done c", str)
    return sg
end
SequenceDistanceGraph(kmerlist::Vector{DNAKmer{K}}) where {K} = new_graph_from_kmerlist(kmerlist)
SequenceDistanceGraph(kmerlist::Vector{DNAKmer{K}},coverage::Vector{Int64}) where {K} = new_graph_from_kmerlist2(kmerlist,coverage)
