#!/usr/bin/env julia
using ArgParse
using BioSequences
using PhyloDataIO
using FASTX

# Time complexity: O(g × t) where g = number of genes and t = number of taxa.
# For each gene, every taxon in the superalignment is visited once (O(t)),
# and taxon membership is checked in O(1) using a Set.
# This is optimal — every taxon-gene combination must be visited at least once.
# Final version of script 01/06/2026. Based on older perl-script. 
# example usage: ConcatAlignments.jl -e fas -o superalignment.phy -f p
# to visualise all options: ConcatAlignments.jl --help  


function main()
    args= parse_arguments()
    extension= args["extension"]
    outfile= args["output"]
    format= args["outformat"]

    base_outfile, extension_outfile = splitext(outfile)
    partition_file_name= "$(base_outfile)_partitions.txt"
 
    list_of_genes= filter(f->endswith(f, extension), readdir("."))
    if isempty(list_of_genes)
        println("No files found with extension $extension")
        exit(1)
    end
    superalignment, gene_names, partition_boundaries= concatenate_sequences(list_of_genes)
    
    write_PhyloData(superalignment, outfile, format)
    open(partition_file_name, "w") do fh 
        for (i, gene) in enumerate(gene_names)
            println(fh, "$(gene) = $(partition_boundaries[i]);")
        end
    end
end

function concatenate_sequences(list::Vector{String})::Tuple{Dict{String,String},Vector{String},Vector{String}}  ## stupid function to give a sensible name in main
    all_my_genes= read_many_sequence_files(list)
    superalignment, genes, partitions= concatenate(all_my_genes)
    return superalignment, genes, partitions 
end

function read_many_sequence_files(list_of_genes::Vector{String})::Dict{String, Dict{String,String}}
    all_genes= Dict{String, Dict{String,String}}()
    for gene in list_of_genes
        sequences= read_PhyloData(gene)
        all_genes[gene]=sequences
    end
    return all_genes
end

function concatenate(sequences::Dict{String, Dict{String,String}})::Tuple{Dict{String,String},Vector{String},Vector{String}}
    genes= sort(collect(keys(sequences)))
    all_taxa= Set{String}()
    super= Dict{String,String}()
    partitions= Vector{String}()
    for gene in genes  # this should create a dictionary with all the taxa.
        current_gene=sequences[gene]
        taxa=sort(collect(keys(current_gene)))
        for taxon in taxa
            push!(all_taxa, taxon)
        end
    end
    taxa_in_superalignment=sort(collect(all_taxa)) # this put all the taxa alphabetically ordered in an array.
    for taxon in taxa_in_superalignment  ## initialiase supermatrix because I need someting in there to append later on (append to nothing is ok as long as its done).                     
        super[taxon] = ""
    end
    current_position=1
    for gene in genes
        current_gene=sequences[gene]
        set_of_taxa_in_gene= Set(keys(current_gene))
        nsites_gene=length(first(values(current_gene)))
        gene_start= current_position
        gene_end = current_position + nsites_gene -1
        push!(partitions, "$gene_start - $gene_end")
        for taxon in taxa_in_superalignment
            if taxon in set_of_taxa_in_gene                
                super[taxon] *= current_gene[taxon]
            else
                super[taxon] *= repeat("?", nsites_gene)
            end
        end
        current_position= gene_end +1
    end
    return super, genes, partitions
end

function parse_arguments()
    s =ArgParseSettings(description="reformat sequences")
    @add_arg_table s begin
        "--extension", "-e"
            help = "extension of the files to merge usually .fas, .phy, .nex"
            required = true
        "--output", "-o"
            help = "output file"
            required = true
        "--outformat", "-f"
            help = "output format: n=nexus; f=fasta; p=phylip"
            required =true
    end
    return parse_args(s)
end

main()
