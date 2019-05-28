

const NodeID = Int64

"""
A SequenceGraphNode represents nodes or vertices in a sequence graph.

Conceptually, a node in a sequence graph is a node which posesses a
biological sequence, such nodes may be connected to other nodes through
'Links'.

SequenceGraphNodes are directional. They have a start (marked '+'),
and they have an end (marked '-').

Furthermore, nodes are canonical or palindromic sequences, else they are reverted.
"""
struct SequenceGraphNode{S <:Sequence}
    sequence::S
    active::Bool
end


sequence(sn::SequenceGraphNode) = sn.sequence

isactive(sn::SequenceGraphNode) = sn.active

reverse_complement!(sn::SequenceGraphNode) = reverse_complement!(sequence(sn))

length(sn::SequenceGraphNode) = length(sequence(sn))

function iscanonical(seq::BioSequence{DNAAlphabet{2}})
    i = 1
    j = length(seq)
    @inbounds while i < j
        f = seq[i]
        r = complement(seq[j])
        f < r && return true
        r < f && return false
        i += 1
        j -= 1
    end
    return true
end

iscanonical(sn::SequenceGraphNode) = iscanonical(sequence(sn))


"""
    is_suffix(seq::Sequence,sn::SequenceGraphNode,min_match = 3)

Check if any prefix of seq is a suffix of sn starting from the longest one
min_match limits the minimum matching and it is the kmer length during dbg contstruction

returns the length of the match, returns -1 otherwise

TO-DO: replace with a faster substring matcher
"""
function is_suffix(seq::Sequence,sn::SequenceGraphNode,min_match = 3)
    seqlen = length(seq)
    snlen  = Base.length(sequence(sn))
    max_match = min(seqlen,snlen)
    for i in max_match:-1:min_match ##possible match length
        if is_match(seq,1,sequence(sn),snlen-i+1,i)

            return i
        end
    end
    return -1
end

"""
    is_in_node(seq::Sequence,sn::SequenceGraphNode)


Checks if a seq is included in the label/sequence of sn as a substring

Returns true if seq is included in the label false o/w
"""
function is_in_node(seq::Sequence,sn::SequenceGraphNode)
    seqlen = length(seq)
    snlen  = Base.length(sequence(sn))
    if seqlen > snlen
        return False
    else
        diff = snlen - seqlen
        sn_str = String(sequence(sn))
        seq_str = String(seq)
        for i in 1:1+diff
            if sn_str[i:i+seqlen-1]==seq_str
                return true
            end
        end
    end
    return false
end
