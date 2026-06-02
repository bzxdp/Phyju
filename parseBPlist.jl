using Statistics


## simple script to extract support from Phylobayes pblist files.  Takes a bplist file, number of chains, clade file with definition of clades for which support needs calculated and outfile
## e.g. julia parseBplist.jl test.bplist 2 clade_list.txt outfile.txt


infile= ARGS[1]
nchains= parse(Int, ARGS[2])
list_of_clades_to_search= ARGS[3]  ## this is a simple file, first name clade followed by definition (species in name).  E.g. CladeA taxonA taxonB taxonC
outfile= ARGS[4]

taxa_in_table = Vector{String}()
numerical_codes_of_taxa_in_table= Dict{String,Int}()
bipart_table = Vector{String}()
bipartitions_and_supports= Vector{Tuple{String,Vector{Int}}}()
clade_to_search_as_a_setof_taxa= Dict{String,Set{Int}}()  ## this has a stupid name because I changed logic. needs remaning 

## read Phylobayes output
open(infile, "r") do fh
    boul_control=false
    for line in eachline(fh)
        line= strip(line)
        if isempty(line)
            continue
        end
        if startswith(line, "Names")
            boul_control=true
            continue # this is key or Names will be parsed
        end
        if startswith(line, "End")
            boul_control=false
            continue
	end
        if boul_control == true
            push!(taxa_in_table, line)
        end
        if startswith(line, ".")
            push!(bipart_table, line)
        end
        if startswith(line, "*")
            push!(bipart_table, line)
        end
    end
end

## fill bipartitions_and_supports
for bipartition in bipart_table
    chain_supports= Int[]
    splitted_lines = String.(split(bipartition, r"\s+"))
    pattern = splitted_lines[1]
    for index in 3:2+nchains
        push!(chain_supports, parse(Int, splitted_lines[index]))
    end
    push!(bipartitions_and_supports, (pattern, chain_supports))
end

# define the number with which species are ordered in the .* string
for (i, taxon) in enumerate(taxa_in_table)
    numerical_codes_of_taxa_in_table[taxon]=i
end

## This define the clades to search as SETS of taxa represented by the numbers in numerical_codes_of_taxa_in_table
## expect a simple format for input file: a clade on each line. First element clade name second element on species names elements separated by blank spaces.
open(list_of_clades_to_search, "r") do fh2
    for line in eachline(fh2)
        line= strip(line)
        if isempty(line)
            continue
        end
        splitted_line=String.(split(line, r"\s+"))  ## remember julia is stringly typed so I need to stringify substrings.  useful to use broadcasting "." as all elements have to be stringified 
        clade_name= popfirst!(splitted_line)
        taxa = Vector{Int}()
        for taxon in splitted_line
            push!(taxa, numerical_codes_of_taxa_in_table[taxon])
        end
        taxa_as_a_set_of_integer = Set(taxa)  
        clade_to_search_as_a_setof_taxa[clade_name]= taxa_as_a_set_of_integer
    end
end

### This transform the strings of *. into splits  of taxa defined as sets we define both side of each split and assign same support.  Important because we do not know if clade of interest coded as * or .
star_defined_bipartitions_and_supports_as_sets_and_supports = Vector{Tuple{Set{Int},Vector{Int}}}()
dot_defined_bipartitions_and_supports_as_sets_and_supports = Vector{Tuple{Set{Int},Vector{Int}}}()

for (pattern, support_by_chains) in bipartitions_and_supports
    star_bipartition=Vector{Int}()
    dot_bipartition=Vector{Int}()
    for (i, char) in enumerate(pattern)
        if char =='*'
            push!(star_bipartition, i)
        else
            push!(dot_bipartition, i)
        end
    end
    star_set=Set(star_bipartition)
    dot_set=Set(dot_bipartition)
    push!(star_defined_bipartitions_and_supports_as_sets_and_supports, (star_set, support_by_chains))
    push!(dot_defined_bipartitions_and_supports_as_sets_and_supports, (dot_set, support_by_chains))        
end

all_clades_and_supports= Vector{Tuple{Set{Int},Vector{Int}}}()  ## we now put all splits (* and . defined) in the same array.
for (clade_set, supports_vector) in star_defined_bipartitions_and_supports_as_sets_and_supports
    push!(all_clades_and_supports, (clade_set, supports_vector))
end
for (clade_set, supports_vector) in dot_defined_bipartitions_and_supports_as_sets_and_supports
    push!(all_clades_and_supports, (clade_set, supports_vector))
end

clades_list= sort(collect(keys(clade_to_search_as_a_setof_taxa))) ## here we make a list of all the clades to search

open(outfile, "w") do fh3
    for clade_to_search in clades_list
        found=false
        for (clade_in_table, supports_vector) in all_clades_and_supports
            if clade_to_search_as_a_setof_taxa[clade_to_search] == clade_in_table
                average_support= mean(supports_vector)
                println(fh3, "Average Support for $clade_to_search in $nchains chains = $average_support -- individual chains support = $(supports_vector)")
                found = true
                break
            end
        end
        if !found
            println(fh3, "Clade $clade_to_search not in list: Support is 0")
        end
    end    
end


        
        
            
            
        










            
                   
        
    
