using ArgParse
using BioSequences
using PhyloDataIO
using FASTX

function main()
    args= parse_arguments()
    extension= args["extension"]
    outfile= args["output"]
    format= args["outformat"]
    partition_file_name= "$(outfile)_partitions.txt"
 
    list_of_genes = filter(f->endswith(f, extension), readdir("."))
    
    superalignment, gene_names, partition_boundaries= concatenate_sequences(list_of_genes)
    
    write_PhyloData(superalignment, outfile, format)
    open (partition_file_name, "w") do fh 
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
    all_taxa= Dict{String,Int}()
    super= Dict{String,String}()
    partitions= Vector{String}()
    for gene in genes  # this should create a dictionary with all the taxa.
        current_gene=sequences[gene]
        taxa=sort(collect(keys(current_gene)))
        for taxon in taxa
            all_taxa[taxon] =1
        end
    end
    taxa_in_superalignment=sort(collect(keys(all_taxa))) # this put all the taxa alphabetically ordered in an array.
    for taxon in taxa_in_superalignment  ## initialiase supermatrix because I need someting in there to append later on (append to nothing is ok as long as its done).                     
        super[taxon] = ""
    end
    current_position=1
    for gene in genes
        current_gene=sequences[gene]
        taxa_in_gene=sort(collect(keys(current_gene)))
        nsites_gene=length(current_gene[taxa_in_gene[1]])
        gene_start= current_position
        gene_end = current_position + nsites_gene
        push!(partitions, "$gene_start - $gene_end")
        for taxon in taxa_in_superalignment
            lineage_is_in_gene= false
            for lineage in taxa_in_gene
                if lineage == taxon
                    lineage_is_in_gene = true
                end
            end
            if lineage_is_in_gene == true
                super[taxon] *= current_gene[taxon]
            end
            if lineage_is_in_gene == false
                super[taxon] *= repeat ("-", nsites_gene)
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

