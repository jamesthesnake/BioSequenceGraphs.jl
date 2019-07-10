

sub_seq(seq::Sequence,ind::Int64,ind2::Int64=Base.length(seq)) = typeof(seq)(String(seq)[ind:ind2])
sub_seq(kmer::Kmer{T,K},ind::Int64,ind2::Int64=K) where{T,K} = Kmer{T,ind2-ind+1}(String(kmer)[ind:ind2])

"""
    get_fastq_sequences(fastq_file,::Type{S})where S <: BioSequences.Sequence


As we have some dependency issues related to the biosequences package
I implement my own  fastq reader just for testing the functionalities
"""
function get_fastq_sequences(fastq_file,::Type{S})where S <: BioSequences.Sequence
    a = readlines(open(fastq_file,"r"));
    println(S)
    seqs = Set{S}()
    for i in 2:4:Base.length(a)
        x  = S(a[i])
        push!(seqs,x)
    end
    return seqs
end



"""
    get_unique_kmers(dbg)

For the contigging stage we first find the unique kmers in the unitig graph.
Then for each read we find the position of each kmer that occur in the read and see  if the path is sensible.


"""

function get_unique_kmers(dbg,k)

end

"""
    get_paired_end_reads(seqs::Set{S})where S <: BioSequences.Sequence


For testing purposes we generate some paired_end reads by masking the middle part of some reads
We represent paired-end reads as a tuple of two biosequences and an integer range representing the possible distance

"""
function get_paired_end_reads(seqs::Set{S})where S <: BioSequences.Sequence
    masked_seqs= Set{Tuple{S,S,UnitRange{Int64}}}()
    for seq in seqs
        upp = Base.length(seq)
        thres = upp/10
        rand1 = rand(1:upp)
        rand2 = rand(1:upp)
        while(abs(rand2-rand1)<thres)
            rand1 = rand(1:upp)
            rand2 = rand(1:upp)
        end
        rng1 = min(rand1,rand2):max(rand1,rand2) ## masked range
        paired_end_read = (seq[1:rng1[1]],seq[rng1[2]+1:end],rng1)
        push!(masked_seqs,paired_end_read)
    end
    return masked_seqs
end
