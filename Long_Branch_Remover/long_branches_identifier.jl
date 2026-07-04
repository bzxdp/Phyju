using TreesIO
using TreeUtils
using Statistics
using TreeStats
using NewickTree
using MetaGraphsNext
using ArgParse

function main()
    args= parse_arguments()
    treefile_extension::String= args["infile__extension"]
    outfile_extension::String= args["outfile__extension"]  ## extension of list of taxa to be excluded
    clades_file::String= args["clades"] ## clade list for long branch tests
    quartet_tests_excluded::Union{String,Nothing}= args["excluded_from_quartet_tests"] ## list of clades not in local quartet tests
    safeguard_quantiles::String= args["safeguard__quantiles"] 
    tresh_triggers::String= args["threshold__triggers"]
    number_of_branches_in_long_path::Int64= args["internal_long_path__trigger"] 

    list_of_treefiles::Vector{String}= filter(f->endswith(f, treefile_extension), readdir("."))  ## array with the list of all treefiles 

    if isempty(list_of_treefiles)
        println("No files found with extension $treefile_extension")
        exit(1)
    end
    trees::Dict{String,MetaGraph}=Dict{String,MetaGraph}() ## This stores all trees as MetaGraphs

    total_files::Int64 = length(list_of_treefiles)
    for (i, file::String) in enumerate(list_of_treefiles)
        println("Processing tree $i of $total_files: $file")
        trees[file]= readtree(file, nothing, "yes")  ## here yes force the check to make sure trees have BLs.  readtree() is a function in TreeIO
    end
            
    
    ### Read file and store taxa to exclude from quartet tests
    clades_to_exclude_from_quartet_tests::Dict{String,Vector{String}}= Dict{String,Vector{String}}()
    if !isnothing(args["excluded_from_quartet_tests"])
        open(quartet_tests_excluded, "r") do fh
            for clade::String in eachline(fh)
                if isempty(strip(clade))
                    continue
                end
                clade_definition::Vector{String}= split(clade, " ")
                clade_name::String= popfirst!(clade_definition)
                clades_to_exclude_from_quartet_tests[clade_name]=clade_definition
            end
        end
    end

    ### Read file and store K values to trigger lb definitions
    triggering_values::Vector{Float64}=Vector{Int64}()
    open(tresh_triggers, "r") do fh
        for line::String in eachline(fh)
            fields::Vector{String}= split(strip(line))
            append!(triggering_values, parse.(Float64, fields))
        end
    end
    ### Read file and store quantiles values to trigger lb retentions                                                                                                                                        
    safeguarding_values::Vector{Float64}=Vector{Float64}()
    open(safeguard_quantiles, "r") do fh
	for line::String in eachline(fh)
            fields::Vector{String}= split(strip(line))
            append!(safeguarding_values, parse.(Float64, fields))
        end
    end
    
    ### This are the vales necesary to identify long branches note as we take these from the stats script they are not all useful we ignore taxon_specific_bls stem_specific_bls bls_and_their_lenght_by_tree.
    ### The values we do not need are replaced when unpacking with _ this is standard Julia to say do not care about that return value. trash it. 

    (global_bls_all_taxa::Set{Float64}, _, taxon_specific_triggering_tresholds::Dict{String,Float64})= stats_for_terminals(trees, triggering_values[1])
    (bls_of_trees, global_set_internal_blens, tree_specific_triggering_thresholds)= stats_for_internal_branches(trees, triggering_values[2])
    (global_bls_all_stems::Set{Float64}, _, stem_specific_triggering_tresholds::Dict{String,Float64})= stats_for_stems(trees, clades_file, triggering_values[3])
    (global_set_of_ratios::Set{Float64}, _)= stats_for_quartets(trees, clades_to_exclude_from_quartet_tests, triggering_values[4])


    terminal_taxa_retention_threshold::Float64 = quantile(collect(global_bls_all_taxa), safeguarding_values[1])
    internal_branches_retention_threshold::Float64 = quantile(collect(global_set_internal_blens), safeguarding_values[2])
    stem_retention_threshold::Float64 = quantile(collect(global_bls_all_stems), safeguarding_values[3])
    quartets_retention_threshold::Float64 = quantile(collect(global_set_of_ratios), safeguarding_values[4])
    
    tree_labels::Vector{String}= collect(keys(trees))

    for tree::String in tree_labels
        outfile= tree * outfile_extension ## name of outfile -- text file with list of taxa to delete
        outtree= tree * ".coloured.nex"

        (terminal_taxa_to_detele_becasue_in_a_clade_on_a_long_branch_path, outremost_branch_in_each_clade)= identify_long_branched_substrees(trees[tree], bls_of_trees[tree], tree_specific_triggering_thresholds[tree], internal_branches_retention_threshold, number_of_branches_in_long_path)
        terminal_taxa_to_delete_because_stem_globally_long= identify_clades_with_excessively_long_stems(trees[tree], tree, stem_specific_triggering_tresholds, stem_retention_threshold, clades_file)
        terminal_taxa_to_delete_because_globally_long::Vector{String} = identify_long_terminals(trees[tree], tree, taxon_specific_triggering_tresholds, terminal_taxa_retention_threshold)
        terminal_taxa_to_delete_because_of_quartet_rule= identify_local_long_branches(trees[tree], tree, clades_to_exclude_from_quartet_tests, triggering_values[4], quartets_retention_threshold)


        ##write out all the long taxa in tree following a specified order long branch subtree, clade, terminal, quartet rule in a  each block separated by a whiteline
        open(outfile, "w") do fh
            for clade::Vector{String} in terminal_taxa_to_detele_becasue_in_a_clade_on_a_long_branch_path
                for taxon::String in clade
                    println(fh, "I: $(taxon)")
                end
            end
            for clade::String in collect(keys(terminal_taxa_to_delete_because_stem_globally_long))
                for taxon::String in terminal_taxa_to_delete_because_stem_globally_long[clade]
                    println(fh, "S: $(taxon)")
                end
            end
            for taxon::String in terminal_taxa_to_delete_because_globally_long
                println(fh, "T: $(taxon)")
            end
            for taxon::String in terminal_taxa_to_delete_because_of_quartet_rule
                println(fh, "Q: $(taxon)")
            end
        end
        
        ### colour tree
        ### rule one - long subtree is ORANGE 
        for edge in outremost_branch_in_each_clade
            (node_a, node_b) = edge
            leaves_a = length([n for n in get_all_nodes_in_clade(trees[tree], node_a, node_b) if trees[tree][n].is_a_leaf])
            leaves_b = length([n for n in get_all_nodes_in_clade(trees[tree], node_b, node_a) if trees[tree][n].is_a_leaf])
            if leaves_a <= leaves_b
                lca_node = node_a
                stem_neighbour = node_b
            else
                lca_node = node_b
                stem_neighbour = node_a
            end
            all_nodes = get_all_nodes_in_clade(trees[tree], lca_node, stem_neighbour)
            for node in all_nodes
                colour_node!(trees[tree], node, "#ff8c00")
            end
        end
        ### rule two - named clade long is BLUE
        for clade::String in keys(terminal_taxa_to_delete_because_stem_globally_long)
            clade_leaves_in_tree::Vector{String} = terminal_taxa_to_delete_because_stem_globally_long[clade]
            lca::Union{String,Nothing} = find_lca(trees[tree], clade_leaves_in_tree) ## this gives us ancestor node of clade
            if isnothing(lca)
                @warn "Could not find LCA for clade $clade in tree $tree — skipping colouring"
                continue
            end
            stem_edge = identify_stem_edge(trees[tree], lca, clade_leaves_in_tree)
            if isnothing(stem_edge)
                @warn "Could not find stem edge for clade $clade in tree $tree — this should not happen. Skipping colouring. Do not trust results for this tree."
                continue
            end
            (node_a, node_b) = stem_edge


            leaves_a = length([n for n in get_all_nodes_in_clade(trees[tree], node_a, node_b) if trees[tree][n].is_a_leaf])
            leaves_b = length([n for n in get_all_nodes_in_clade(trees[tree], node_b, node_a) if trees[tree][n].is_a_leaf])
            if leaves_a <= leaves_b
                lca_node_colour = node_a
                stem_neighbour_colour = node_b
            else
                lca_node_colour = node_b
                stem_neighbour_colour = node_a
            end
            all_nodes = get_all_nodes_in_clade(trees[tree], lca_node_colour, stem_neighbour_colour)
            for node in all_nodes
                colour_node!(trees[tree], node, "#0072B2")
            end
        end
        for taxon::String in terminal_taxa_to_delete_because_globally_long
            colour_node!(trees[tree], taxon, "#ff0000")
        end
        for taxon::String in terminal_taxa_to_delete_because_of_quartet_rule
            colour_node!(trees[tree], taxon, "#D55E00")
        end
        write_nexus(trees[tree], outtree)
    end
end
        
### FUNCTIONS
### ALL FUNCTION DEVELOPED FOR THIS SCRIPT MOVED TO PHYLOSTATS.JL

     
function parse_arguments()
    s =ArgParseSettings(description="reformat sequences")
    @add_arg_table s begin
        "--infile__extension", "-e"
        help = "extension of the files is usually: tree, tre, treefile or similar"
        required = true
        "--outfile__extension", "-o"
        help = "output file"
        required = true
        "--clades", "-c"
        help = "file with clades to test: one line for clade. Name of clade first then all taxa in clade. all separated by spaces (no spaces allowed in names of clades of species: e.g. Clade_A NOT Clade A and Taxon_1 NOT Taxon 1)." 
        required = false
        "--excluded_from_quartet_tests", "-x"
        help = "file with clades expected to be naturally long branch — excluded from quartet tests"
        required = false
        "--safeguard__quantiles", "-s"
        help = "file with quantiles below which clade retained this flie has one row with three numerical values e.g. 0.95 0.95 0.95 0.95 - respectively terminal taxa retention threshold, subtreeretention threshold, names clade retention threshold. If clade_file not passed only two vales need to be in this file."
        required = true
        "--threshold__triggers", "-t"
        help = "file with trigger value (K) used to calculate when a branch is long suggested default value 3. This file has three values separated by withespaces (e.g. 3 3 3 3) - respectively terminal taxa trigger threshold, subtree trigger threshold, named clades trigger threshold. If clade_file not passed only two vales need to be in this file."
        required = true
        "--internal_long_path__trigger", "-p"
        help = "how many long branched internal branches (stacked one on the other) are necessary to define a long path and associated clade for removal - MIN = 1)"
        required = true
        arg_type = Int
        default = 1
    end
    return parse_args(s)
end

main()
