using ArgParse
using BioSequences
using PhyloDataIO
using FASTX

function main()
    args= parse_arguments()
    infile= args["input"]
    length_jakked_data = args["n_jack_sites"] !== nothing ? parse(Int, args["n_jack_sites"]) : nothing # remember this is a ternary operatior if then in a line if true do ? else do :
    nreps = parse(Int, (args["replicates"]))
    outfile_format=args["format"]

    basename, extension	= splitext(infile)
    sequences=  read_PhyloData(infile)
    resampling_type = nothing
    
    
    for replica in 1:nreps
        if args["n_jack_sites"] !== nothing
            taxa, data= jack_sites(sequences, length_jakked_data)
            resampling_type = "jack"
        else
            taxa, data= boot_sites(sequences)
            resampling_type = "boot"
        end
        
        resampled_seqs = matrix_to_dict(taxa, data)
        if outfile_format == "p"
            outfile_name = string(resampling_type, "_", replica, "_", basename, ".phy")
        elseif outfile_format == "n"
            outfile_name = string(resampling_type, "_", replica, "_", basename, ".nex")
        else
            outfile_name = string(resampling_type, "_", replica, "_", basename, ".fas")
        end
        write_PhyloData(resampled_seqs, outfile_name, outfile_format)
    end
end

function TaxaMatrixTuple(seqs::Dict{String,String})::Tuple{Vector{String},Matrix{Char}}
    taxa= sort(collect(keys(seqs)))
    ntaxa = length(taxa)
    nsites = length(seqs[taxa[1]])
    sequences_as_a_matrix= Matrix{Char}(undef, ntaxa, nsites)
    for (i, taxon) in enumerate(taxa)
        for j in 1:nsites
            sequences_as_a_matrix[i, j] = seqs[taxon][j]
        end
    end
    return taxa, sequences_as_a_matrix
end

function matrix_to_dict(taxa::Vector{String}, mat::Matrix{Char})::Dict{String,String}
    seqs = Dict{String,String}()
    for (i, taxon) in enumerate(taxa)
        seqs[taxon] = String(mat[i, :])
    end
    return seqs
end


function jack_sites(sequences::Dict{String,String}, nsiteJack::Int)::Tuple{Vector{String},Matrix{Char}}
    (taxa, phylo_matrix) =  TaxaMatrixTuple(sequences)
    ntaxa, nsites = size(phylo_matrix)
    jacknifed_matrix= Matrix{Char}(undef, ntaxa, nsiteJack)
    sites_in_jack_data = Int[]
    while length(sites_in_jack_data) < nsiteJack
        is_site_in_jack_data = false
        current_selected_site= rand(1:nsites)
        for site in sites_in_jack_data
            if site == current_selected_site
                is_site_in_jack_data = true
            end
        end
        if is_site_in_jack_data == true
            continue
        end
        j_index_to_populate = length(sites_in_jack_data) +1
        push! (sites_in_jack_data, current_selected_site)
        jacknifed_matrix[:,j_index_to_populate] = phylo_matrix[:, current_selected_site]
    end
    return taxa,jacknifed_matrix
end

function boot_sites(sequences::Dict{String,String})::Tuple{Vector{String},Matrix{Char}}
    (taxa, phylo_matrix) =  TaxaMatrixTuple(sequences)
    ntaxa, nsites = size(phylo_matrix)
    booted_matrix= Matrix{Char}(undef, ntaxa, nsites)
    j_index=1
    while j_index <= nsites
        booted_site= rand(1:nsites)
        booted_matrix[:,j_index] = phylo_matrix[:, booted_site]
        j_index +=  1
    end
    return taxa,booted_matrix
end

function parse_arguments()
    s =ArgParseSettings(description="reformat sequences")
    @add_arg_table s begin
        "--input", "-i"
            help = "input file"
            required = true
        "--replicates", "-r"
            help = "number of jackknife or bootstrap replicate datasets to build"
            required =true
        "--n_jack_sites", "-n"
            help = "number of sites in jackknifed datasets"
            required = false
            default = nothing
        "--format", "-f"
            help = "format of outfiles can be f=fasta, p=phylip, n=nexus"
            required =true
    end
    return parse_args(s)
end

main()
