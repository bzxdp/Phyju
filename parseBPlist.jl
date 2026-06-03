#!/usr/bin/env julia
using Statistics
using ArgParse


function main()
    args= parse_arguments()
    bp_list_file= args["infile"]
    outfile= args["outfile"]
    nchains= args["number_of_chains"]
    clade_list =args["clades_to_test"]

    
    (taxa_in_tree, table_rows) = readbplist(bp_list_file)
    numerical_codes_of_taxa_in_table= enumerate_taxa(taxa_in_tree)
    star_dots_patterns_and_supports= parse_table_row(table_rows,nchains) ## This has all lines in table as the pattern representing bipartitions and support vals 
    clade_to_search_as_sets_of_taxa= define_clades_to_search_as_sets(clade_list, numerical_codes_of_taxa_in_table)
    all_clades_and_supports= single_out_bipartitions_and_their_support(star_dots_patterns_and_supports)  ## this has all bipartitions in table (i.e. every *** and ... pattern from each row -i.e. part and counterpart- as an individually listed set of taxa) with associated support
    print_results(outfile, all_clades_and_supports, clade_to_search_as_sets_of_taxa, nchains)
end

## read Phylobayes output
function readbplist(infile::String)::Tuple{Vector{String},Vector{String}}
    taxa_in_table= Vector{String}()
    bipart_table= Vector{String}()
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
    return taxa_in_table, bipart_table
end

function parse_table_row(bipart_table::Vector{String},number_of_chains::Int)::Vector{Tuple{String,Vector{Int}}}
    pattern_and_chain_supports= Vector{Tuple{String,Vector{Int}}}()
    for bipartition in bipart_table
        chain_supports= Int[]
        splitted_lines = String.(split(bipartition, r"\s+"))
        pattern = splitted_lines[1]
        for index in 3:2+number_of_chains
            push!(chain_supports, parse(Int, splitted_lines[index]))
        end
        push!(pattern_and_chain_supports, (pattern, chain_supports))
    end
    return pattern_and_chain_supports
end

# define the number with which species are ordered in the .* string
function enumerate_taxa(taxa_in_pblist::Vector{String})::Dict{String,Int}
    taxa_and_their_order= Dict{String,Int}()
    for (i, taxon) in enumerate(taxa_in_pblist)
        taxa_and_their_order[taxon]=i
    end
    return taxa_and_their_order
end

## This define the clades to search as SETS of taxa represented by the numbers in numerical_codes_of_taxa_in_table
## expect a simple format for input file: a clade on each line. First element clade name second element on species names elements separated by blank spaces.

function define_clades_to_search_as_sets(infile::String,numerical_codes::Dict{String,Int})::Dict{String,Set{Int}}
    clades_as_sets=Dict{String,Set{Int}}()
    open(infile, "r") do fh2
        for line in eachline(fh2)
            line= strip(line)
            if isempty(line)
                continue
            end
            splitted_line=String.(split(line, r"\s+"))  ## remember julia is stringly typed so I need to stringify substrings.  useful to use broadcasting "." as all elements have to be stringified 
            clade_name= popfirst!(splitted_line)
            taxa = Vector{Int}()
            for taxon in splitted_line
                push!(taxa, numerical_codes[taxon])
            end
            taxa_as_a_set_of_integer = Set(taxa)  
            clades_as_sets[clade_name]= taxa_as_a_set_of_integer
        end
    end
    return clades_as_sets
end    


### This transform the strings of *. into splits  of taxa defined as sets we define both side of each split and assign same support.  Important because we do not know if clade of interest coded as * or .
function single_out_bipartitions_and_their_support(bipartitions_and_supports_as_star_dots::Vector{Tuple{String,Vector{Int}}})::Vector{Tuple{Set{Int},Vector{Int}}}
    bipartitions_and_support_to_return=Vector{Tuple{Set{Int},Vector{Int}}}()
    star_defined_bipartitions_and_supports_as_sets_and_supports = Vector{Tuple{Set{Int},Vector{Int}}}()
    dot_defined_bipartitions_and_supports_as_sets_and_supports = Vector{Tuple{Set{Int},Vector{Int}}}()
    for (pattern, support_by_chains) in bipartitions_and_supports_as_star_dots
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
    for (clade_set, supports_vector) in star_defined_bipartitions_and_supports_as_sets_and_supports
        push!(bipartitions_and_support_to_return, (clade_set, supports_vector))
    end
    for (clade_set, supports_vector) in dot_defined_bipartitions_and_supports_as_sets_and_supports
        push!(bipartitions_and_support_to_return, (clade_set, supports_vector))
    end
    return bipartitions_and_support_to_return
end

function print_results(outfile::String,clades_and_their_support::Vector{Tuple{Set{Int},Vector{Int}}},list_of_clades_to_map::Dict{String,Set{Int}},number_of_chains::Int) 
    clades_list= sort(collect(keys(list_of_clades_to_map))) ## here we make a list of all the clades to search
    
    open(outfile, "w") do fh3
        for clade_to_search in clades_list
            found=false
            for (clade_in_table, supports_vector) in clades_and_their_support
                if list_of_clades_to_map[clade_to_search] == clade_in_table
                    average_support= mean(supports_vector)
                    println(fh3, "Average Support for $clade_to_search in $number_of_chains chains = $average_support -- individual chains support = $(supports_vector)")
                    found = true
                    break
                end
            end
            if !found
                println(fh3, "Clade $clade_to_search not found: Support is 0")
            end
        end    
    end
end


function parse_arguments()
    s =ArgParseSettings(description="reformat sequences")
    @add_arg_table s begin
        "--infile", "-i"
            help = "phylobayes output: *.bplist file"
            required = true
        "--outfile", "-o"
            help = "output file"
            required = true
        "--number_of_chains", "-n"
            help = "number of chains summarised in *bplist"
            arg_type= Int
            required =true
        "--clades_to_test", "-c"
            help = "text file with list of clades to test.  Defined as. Each clade in one row.  Start with clade name follow all species in clade. all separated by spaces.  e.g. CLADEa TaxonA TaxonB TaxonC TaxonD"
            required=true
    end
    return parse_args(s)
end

main()
        
        
             
        










            
                   
        
    
