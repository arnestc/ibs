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

    return FastaIO.readfasta(seqFile)[1][2]
end

"""
    hammington(x,y)

Compute the Hammington distance between two same length sequences `x` and `y`.
The Hammington distance between to equal length sequences simply counts the number of minima substitutions of characters needed to transform one sequence into the other.

# Example
```julia-repl
julia> x = "wheeaaa"
julia> y = "ghearpa"
julia> d = hammington(x,y)
4
```
"""
function hammington(x,y)
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