using FastaIO
using DelimitedFiles
using Printf
using HTTP

"""
    sequenceDownload(sequence)

Download the (protein) sequence with identifying code `sequence` in Fasta format (https://en.wikipedia.org/wiki/FASTA_format).
The site where the sequence is located is "https://www.uniprot.org/uniprotkb/sequence.fasta": for example, the protein HBB_Human has sequence code "P68871", so the site at which protein sequence is located is "https://www.uniprot.org/uniprotkb/P68871.fasta".


# Example
```julia-repl
julia> HBB_Human = sequenceDownload("P68871")
"MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH"
```
"""
function sequenceDownload(sequence)

    sequenceFile = sequence * ".fasta"

    URL = "https://www.uniprot.org/uniprotkb/" * sequenceFile

    query = HTTP.get(URL)
    fastaString=String(query.body)

    open(sequenceFile,"w") do f
        write(f,fastaString)
    end

    return FastaIO.readfasta(sequenceFile)[1][2]
end

"""
    hamming(x,y)

Compute the Hamming distance between two same length sequences `x` and `y`.
The Hamming distance between to equal length sequences simply counts the number of minima substitutions of characters needed to transform one sequence into the other.

# Example
```julia-repl
julia> x = "wheeaaa"
julia> y = "ghearpa"
julia> d = hamming(x,y)
4
```
"""
function hamming(x,y)
    # x, y sequences
    if length(x) != length(y)
        # Sequences must have equal lenghts
        println("ERROR: sequences should have equal lenghts!")
        return
    else
        d = 0 # Starting sequences distance initialized to 0
        for i = 1:eachindex(x)
            if x[i] != y[i] # If x[i] != y[i] then it occurs a character substitution, increasing the distance
                d += 1
            end
        end
        return d
    end
end

"""
    leven(x,y)

Compute the Levenshtein distance between two sequences `x` and `y` (also if they have different lengtsh) in a recursive way (https://wikimedia.org/api/rest_v1/media/math/render/svg/70962a722b0b682e398f0ee77d60c714a441c54e) throughout the use of a Dictionary data structure.
# Example
```julia-repl
julia> x = "wheeaaa"
julia> y = "ghearpa"
julia> d = lev(x,y)
4
```
"""
function leven(x,y)
    D = Dict() # Dict data struct used to implement recursive approach

    function levenshtein(x,y)
        isempty(x) && return length(y) # Base case 1
        isempty(y) && return length(x) # Base case 2
        haskey(D,(x,y)) && return D[(x,y)] # Base case 3 (residues comparison already computed, then no need to compute it again)
        D[(x,y)] = min(1 - (x[end] == y[end]) + levenshtein(x[1:end-1],y[1:end-1]), 1 + levenshtein(x[1:end-1],y), 1 + levenshtein(x,y[1:end-1]))
    end

    levenshtein(x,y)
end