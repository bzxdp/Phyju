using Logging
using BioSequences

function main()
    if length(ARGS) != 2
        println(stderr, "Usage: julia entropy.jl input.phy output.txt")
        exit(1)
    end

    if isfile(ARGS[1]) == false
        @error "File not found: $(ARGS[1])"
        exit(1)
    end
    
    seqs_data= read_phylip(ARGS[1])
    entropy= sites_entropy(seqs_data)
    keffs= sites_keffs(entropy)
    sites= sort(collect(keys(entropy)))
    open(ARGS[2], "w") do fh
        println(fh, "site\tEntropy (bits)\tkeff")
        for site in sites
            println(fh, "$site\t$(entropy[site])\t$(keffs[site])")
        end
    end    
end

function sites_keffs(entropy::Dict{Int,Float64})::Dict{Int,Float64}
    sites_entropies= entropy
    keffs = Dict{Int,Float64}()
    sites = sort(collect(keys(sites_entropies)))
    for site in sites
        keffs[site]= 2^sites_entropies[site]
        @info "Keff of site $site = $(keffs[site])"
    end
    return keffs
end
    
function read_phylip(filename::String)::Dict{String,LongAA}    
    seqs = Dict{String,LongAA}()
    ntaxa = nothing
    nsites = nothing

    line_counter = 0
    open(filename, "r") do fh
        for line in eachline(fh)
            line = strip(line)
            if isempty(line)
                continue
            end
            line_counter +=1

            if line_counter ==1
                current_line = split(line, r"\s+", limit=2)
                ntaxa = parse(Int, String(current_line[1]))
                nsites = parse(Int, String(current_line[2]))
                @info "Expecting $ntaxa taxa and $nsites sites" 
            else
                current_line = split(line, r"\s+", limit=2)
                seqs[String(current_line[1])] = LongAA(String(current_line[2]))
                @info "Read sequence $(current_line[1])"
            end
        end
    end
    return seqs
end

function sites_entropy(sequences::Dict{String,LongAA})::Dict{Int,Float64}
    seqs= sequences
    entropy =Dict{Int,Float64}()
    taxa =collect(keys(seqs))
    nsites =length(seqs[taxa[1]])

    
    for site in 1:nsites
        AA_at_site=	Dict{AminoAcid,Int}()
        for taxon in taxa
            char= seqs[taxon][site]
            if isgap(char) || isambiguous(char) 
                continue
            end
            if haskey(AA_at_site, char)
                AA_at_site[char] +=1
            else
                AA_at_site[char] =1
            end
        end

        AA_list =collect(keys(AA_at_site))
        AA_freqs_at_site=Dict{AminoAcid,Float64}()
        total= 0
        for AA in AA_list
            total += AA_at_site[AA]
        end

        if total == 0
            entropy[site] = 0.0
            continue
        else
            for AA in AA_list
                AA_freq= AA_at_site[AA]/total
                AA_freqs_at_site[AA]= AA_freq
            end
        end
        
        entropy_site=0
        for AA in AA_list
            entropy_site += AA_freqs_at_site[AA]*log2(AA_freqs_at_site[AA])
        end
        entropy[site]= max(0.0, -entropy_site) # this is just to avoid in some instances -0.0 (foating point quirk).
        @info "Entropy of site $site = $(entropy[site])"
    end
    return entropy
end     
main()
