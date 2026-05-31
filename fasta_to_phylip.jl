#!/usr/bin/env julia
using Logging


### generic stripped down fasta to phylip useful because it does a single job and does not need flags. first argument fasta file. second argument outfile that will be in phylip format
### usage julia fasta_to_phylip.jl infile outfile

function main()
    if length(ARGS) != 2
        println(stderr, "Usage: julia fasta_to_phylip.jl input.fasta output.phy")
        exit(1)
    end

    if isfile(ARGS[1]) == false
        @error "File not found: $(ARGS[1])"
        exit(1)
    end
    
    sequences= (read_fasta(ARGS[1]))
    taxa= sort(collect(keys(sequences)))
    ntaxa= length(taxa)
    nsites= length(sequences[taxa[1]])
    open(ARGS[2], "w") do fh
        println(fh, "$ntaxa $nsites")
        for taxon in taxa
            print(fh, "$taxon  ")
            println(fh, sequences[taxon])
        end
    end
end

function read_fasta(filename::String)::Dict{String,String}
    sequences = Dict{String, String}()
    seq_id= nothing
    seq_line_chunks= String[] 
    open(filename, "r") do fh
        for line in eachline(fh)
            line= strip(line)
            if isempty(line)
                continue
            end
            
            if startswith(line, ">")   ### condition when seq_id is not empty (second sequence name onward)
                if seq_id !== nothing
                    sequences[seq_id]= join(seq_line_chunks)
                    seq_line_chunks= String[]
                end
                seq_id = replace(line, ">"=> "", r"\s+"=> "_")   #This line take care of all lines starting with > including first one because block jumped on first loop
                @info "found $seq_id" 
            else
                push!(seq_line_chunks, line)
            end
        end
        if seq_id !== nothing
            sequences[seq_id] = join(seq_line_chunks) #take care of last sequence that is added to dictionary here or will be missed
        else
            @error "file $filename is not Fasta format"
        end
    end
    return sequences
end
main()
