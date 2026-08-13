using TreesIO
using TreeUtils
using Statistics
using TreeStats
using NewickTree
using MetaGraphsNext
using ArgParse
using CairoMakie

### The long-branch rules themselves live next door, not in the libraries: they are
### specific to this analysis (K triggers, safeguard quantiles, clade files) and were
### cluttering TreeStats/TreeUtils, which are shared with other projects.
include("long_branch_rules.jl")

function main()
    args= parse_arguments()
    treefile_extension::String= args["infile__extension"]
    outfile_extension::String= args["outfile__extension"]  ## extension of list of taxa to be excluded
    clades_file::String= args["clades"] ## clade list for long branch tests
    quartet_tests_excluded::Union{String,Nothing}= args["excluded_from_quartet_tests"] ## list of clades not in local quartet tests
    safeguard_quantiles::String= args["safeguard__quantiles"] 
    tresh_triggers::String= args["threshold__triggers"]
    number_of_branches_in_long_path::Int64= args["internal_long_path__trigger"]
    internal_side_ratio_threshold::Float64= args["internal_side_ratio"]

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
    triggering_values::Vector{Float64}=Vector{Float64}()
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

    ### clade definition sizes — needed only to decide, for reporting, which name in
    ### a set of aliases is the "representative": the one whose FULL definition is
    ### smallest, i.e. the clade that most tightly describes the branch.  The others
    ### are travelling with it because they resolved to the same branch in that tree.
    clade_definition_sizes::Dict{String,Int64} = Dict{String,Int64}()
    open(clades_file, "r") do fh
        for line::String in eachline(fh)
            fields::Vector{String} = split(strip(line))
            if isempty(fields)
                continue
            end
            name::String = popfirst!(fields)
            clade_definition_sizes[name] = length(fields)
        end
    end

    ### deletion tallies for the clade deletion figures, accumulated across all trees.
    ### Two parallel accounts:
    ###   *_deletions          — clade lost by EITHER rule (stem or internal)
    ###   *_deletions_internal — clade lost by the INTERNAL rule alone, so a run can
    ###                          answer "would the internal rule be enough on its own,
    ###                          and could the named-clade rule be dropped?"
    clade_rep_deletions::Dict{String,Int64} = Dict{String,Int64}()
    clade_alias_deletions::Dict{String,Int64} = Dict{String,Int64}()
    clade_rep_deletions_internal::Dict{String,Int64} = Dict{String,Int64}()
    clade_alias_deletions_internal::Dict{String,Int64} = Dict{String,Int64}()

    ### fail early and clearly rather than with a BoundsError further down
    if length(triggering_values) < 4
        println("ERROR: $(tresh_triggers) must contain 4 K values (terminals, internals, stems, quartets) — found $(length(triggering_values)).")
        exit(1)
    end
    if length(safeguarding_values) < 3
        println("ERROR: $(safeguard_quantiles) must contain 3 quantiles (terminals+quartets, internals, stems) — found $(length(safeguarding_values)).")
        exit(1)
    end
    if number_of_branches_in_long_path < 1
        println("ERROR: --internal_long_path__trigger must be at least 1 — got $(number_of_branches_in_long_path).")
        exit(1)
    end
    if internal_side_ratio_threshold <= 1.0
        println("ERROR: --internal_side_ratio must be greater than 1.0 — got $(internal_side_ratio_threshold). A value of 1.0 or less accepts any difference between the two sides and disables the asymmetry test.")
        exit(1)
    end

    ### This are the vales necesary to identify long branches note as we take these from the stats script they are not all useful we ignore taxon_specific_bls stem_specific_bls bls_and_their_lenght_by_tree.
    ### The values we do not need are replaced when unpacking with _ this is standard Julia to say do not care about that return value. trash it. 

    (global_bls_all_taxa::Vector{Float64}, _, taxon_specific_triggering_tresholds::Dict{String,Float64})= stats_for_terminals(trees, triggering_values[1])
    (bls_of_trees, global_set_internal_blens, tree_specific_triggering_thresholds)= stats_for_internal_branches(trees, triggering_values[2])
    (global_bls_all_stems::Vector{Float64}, _, stem_specific_triggering_tresholds::Dict{String,Float64})= stats_for_stems(trees, clades_file, triggering_values[3])
    ## NOTE: stats_for_quartets() is deliberately NOT called here.  Its only product
    ## was the global distribution of quartet ratios, which used to supply the
    ## quartet safeguard.  That safeguard is now the terminal branch length quantile
    ## (see terminal_quartet_safeguard below), so computing the ratio distribution
    ## would mean running the quartet test over every taxon in every tree for a
    ## number nobody reads.  The ratio distribution is still produced, and plotted,
    ## by long_branches_statistical_analyses.jl.


    ## terminal_quartet_safeguard serves BOTH the terminal rule and the quartet rule:
    ## a quantile of every terminal branch length in the dataset.  One value, two
    ## consumers — as in Python, where --global_rescue rescues cross-gene terminal
    ## drops and quartet drops alike.  This is why safeguarding_values has 3 entries
    ## and not 4.
    terminal_quartet_safeguard::Float64 = quantile(global_bls_all_taxa, safeguarding_values[1])

    ## A comparator sitting on a branch the inference could not estimate carries no
    ## information about the local scale, and lets the quartet ratio explode.  The
    ## threshold scales with the dataset rather than sitting at a fixed floor.
    comparator_branch_floor::Float64 = comparator_branch_threshold(global_bls_all_taxa)
    println("\nQuartet comparators must sit on a terminal branch longer than $(comparator_branch_floor)")
    internal_branches_retention_threshold::Float64 = quantile(global_set_internal_blens, safeguarding_values[2])
    stem_retention_threshold::Float64 = quantile(global_bls_all_stems, safeguarding_values[3])
    
    tree_labels::Vector{String}= collect(keys(trees))

    ### DATASET SIZE GUARDS (Python CROSS_GENE_DISABLE_TREES / WARN_VERY_FEW_TREES /
    ### WARN_FEW_TREES).  The terminal, clade-stem and internal-branch rules all judge
    ### a branch against a median taken ACROSS trees.  With a handful of trees that
    ### median is meaningless, so those three rules are switched off entirely rather
    ### than run on noise — in doubt, retain.  The quartet rule is unaffected: it is
    ### a within-tree comparison and needs no cross-gene background.
    cross_gene_enabled::Bool = total_files > 3
    if !cross_gene_enabled
        println("\nWARNING: only $(total_files) gene tree(s) found. The terminal, clade-stem and internal-branch rules are DISABLED because the dataset is too small for cross-gene medians to mean anything. Only the quartet rule will run.")
    elseif total_files < 10
        println("\nWARNING: only $(total_files) gene trees. Cross-gene statistics are very unreliable at this sample size and may be misleading.")
    elseif total_files < 50
        println("\nWarning: cross-gene statistics are based on fewer than 50 gene trees ($(total_files)). Results may be unstable.")
    end

    for tree::String in tree_labels
        outfile= tree * outfile_extension ## name of outfile -- text file with list of taxa to delete
        outtree= tree * ".coloured.nex"

        ## the three cross-gene rules run only when there are enough trees to
        ## support a cross-gene median (see cross_gene_enabled above)
        terminal_taxa_to_detele_becasue_in_a_clade_on_a_long_branch_path::Vector{Vector{String}} = Vector{Vector{String}}()
        outremost_branch_in_each_clade::Vector{Tuple{String,String}} = Vector{Tuple{String,String}}()
        terminal_taxa_to_delete_because_stem_globally_long::Dict{String,Vector{String}} = Dict{String,Vector{String}}()
        terminal_taxa_to_delete_because_globally_long::Vector{String} = Vector{String}()

        if cross_gene_enabled
            (terminal_taxa_to_detele_becasue_in_a_clade_on_a_long_branch_path, outremost_branch_in_each_clade)= identify_long_branched_substrees(trees[tree], bls_of_trees[tree], tree_specific_triggering_thresholds[tree], internal_branches_retention_threshold, internal_side_ratio_threshold, number_of_branches_in_long_path)
            terminal_taxa_to_delete_because_stem_globally_long= identify_clades_with_excessively_long_stems(trees[tree], tree, stem_specific_triggering_tresholds, stem_retention_threshold, clades_file)
            terminal_taxa_to_delete_because_globally_long = identify_long_terminals(trees[tree], tree, taxon_specific_triggering_tresholds, terminal_quartet_safeguard)
        end
        ## quartet rule: triggering_values[4] is K, the LOCAL test.  The safeguard is
        ## terminal_quartet_safeguard — the global TERMINAL BRANCH LENGTH quantile,
        ## the same value the terminal rule uses (Python --global_rescue).  NOT a
        ## quantile of ratios: a trivial branch in a short-branched neighbourhood can
        ## score a huge ratio and is not a long-branch-attraction risk.
        terminal_taxa_to_delete_because_of_quartet_rule= identify_local_long_branches(trees[tree], tree, clades_to_exclude_from_quartet_tests, triggering_values[4], terminal_quartet_safeguard; min_comparator_branch = comparator_branch_floor)


        ##write out all the long taxa in tree following a specified order long branch subtree, clade, terminal, quartet rule in a  each block separated by a whiteline
        open(outfile, "w") do fh
            for clade::Vector{String} in terminal_taxa_to_detele_becasue_in_a_clade_on_a_long_branch_path
                for taxon::String in clade
                    println(fh, "I: $(taxon)")
                end
            end
            ## bare taxon names only — this file is consumed to prune alignments,
            ## so it must contain taxon names and nothing else.  Which clade fired is
            ## reported separately, in the clade deletion table and figure.
            for clade::String in sort(collect(keys(terminal_taxa_to_delete_because_stem_globally_long)))
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
        
        ### tally clade deletions for the figure.
        ### A named clade counts as deleted in this tree if EITHER the stem rule
        ### removed it, OR the internal-branch rule removed a subtree that contains
        ### all of its taxa — if a clade is lost to a long internal branch it is
        ### still lost.  The two sources are unioned per branch so a clade removed
        ### both ways is counted once for this tree, not twice.
        ### Within each branch the name with the LARGEST full definition counts as
        ### the clade itself — the branch stands for that clade in this tree — and the
        ### narrower names count as aliases, going along because they coincided with
        ### it here.  Note this is the opposite of the rule used to pick the TRIGGER,
        ### which comes from the tightest name because its median is the comparable
        ### one.  Naming and calibration are different questions.
        if cross_gene_enabled
            names_by_removed_branch::Dict{String,Set{String}} = Dict{String,Set{String}}()

            ## (a) branches removed by the stem rule — the key already lists every
            ##     clade name that resolved to the branch, the value is its taxa
            for clade_key::String in keys(terminal_taxa_to_delete_because_stem_globally_long)
                branch_key::String = join(sort(terminal_taxa_to_delete_because_stem_globally_long[clade_key]), "\u0001")
                union!(get!(names_by_removed_branch, branch_key, Set{String}()), String.(split(clade_key, "|")))
            end

            ## (b) clades swallowed whole by the internal-branch rule
            taxa_removed_by_internal_rule::Set{String} = Set{String}()
            for clade_taxa::Vector{String} in terminal_taxa_to_detele_becasue_in_a_clade_on_a_long_branch_path
                union!(taxa_removed_by_internal_rule, Set{String}(clade_taxa))
            end

            names_by_internal_branch::Dict{String,Set{String}} = Dict{String,Set{String}}()

            if !isempty(taxa_removed_by_internal_rule)
                (fragments_in_this_tree, _, _) = resolve_clade_fragments(trees[tree], tree, clades_file, collect(keys(stem_specific_triggering_tresholds)))
                for (clade_name, fragment_taxa) in fragments_in_this_tree
                    if isempty(fragment_taxa)
                        continue
                    end
                    if issubset(Set{String}(fragment_taxa), taxa_removed_by_internal_rule)
                        branch_key_internal::String = join(sort(fragment_taxa), "\u0001")
                        push!(get!(names_by_removed_branch, branch_key_internal, Set{String}()), clade_name)
                        push!(get!(names_by_internal_branch, branch_key_internal, Set{String}()), clade_name)
                    end
                end
            end

            tally_clade_deletions!(clade_rep_deletions, clade_alias_deletions, names_by_removed_branch, clade_definition_sizes)
            tally_clade_deletions!(clade_rep_deletions_internal, clade_alias_deletions_internal, names_by_internal_branch, clade_definition_sizes)
        end

        ### colour tree
        ### PRECEDENCE: later colour_node! calls overwrite earlier ones, so the
        ### blocks below run from weakest to strongest claim.  Named clade (BLUE)
        ### is painted FIRST and the long subtree (ORANGE) after it, so the
        ### internal-branch rule TRUMPS the clade rule: whether a named clade is
        ### 'long' is partly a judgement call, whereas a long internal branch is
        ### measured from the tree itself and is the more defensible signal.
        ### Terminal (RED) and quartet (VERMILLION) follow and mark single taxa.
        ### The colours are a visual guide only — the .out file is the record of
        ### which rules fired, and a taxon may appear under more than one.
        ### rule one - named clade long is BLUE
        for clade::String in keys(terminal_taxa_to_delete_because_stem_globally_long)
            clade_leaves_in_tree::Vector{String} = terminal_taxa_to_delete_because_stem_globally_long[clade]
            ## these taxa are a fragment and so should be monophyletic — if they are not,
            ## something is wrong, so the warning is left ON here.  Name passed so it says which clade.
            lca::Union{String,Nothing} = find_lca(trees[tree], clade_leaves_in_tree, clade) ## this gives us ancestor node of clade
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
        ### rule two - long subtree is ORANGE (overwrites BLUE where they overlap)
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
        for taxon::String in terminal_taxa_to_delete_because_globally_long
            colour_node!(trees[tree], taxon, "#ff0000")
        end
        for taxon::String in terminal_taxa_to_delete_because_of_quartet_rule
            colour_node!(trees[tree], taxon, "#D55E00")
        end
        write_nexus(trees[tree], outtree)
    end

    ### -----------------------------------------------------------------------
    ### CLADE DELETION REPORTS
    ### Two figures, deliberately:
    ###   1. every named clade lost to EITHER rule
    ###   2. named clades lost to the INTERNAL-BRANCH rule only
    ### Comparing them shows how much the named-clade (stem) rule is adding over
    ### what the topology-driven internal rule already catches — i.e. whether the
    ### stem rule could be dropped.
    ### The .out files carry taxon names for pruning; this is the clade-level view.
    ### -----------------------------------------------------------------------
    write_clade_deletion_report(clade_rep_deletions, clade_alias_deletions,
                                "clade_deletion_counts.tsv",
                                "clade_deletion_frequencies.png",
                                "Clade deletion frequency (stem rule and internal-branch rule)",
                                "any rule")

    write_clade_deletion_report(clade_rep_deletions_internal, clade_alias_deletions_internal,
                                "clade_deletion_by_internal_rule_counts.tsv",
                                "clade_deletion_by_internal_rule.png",
                                "Clade deletion frequency (internal-branch rule only)",
                                "the internal-branch rule")
end
        
### FUNCTIONS
### ALL FUNCTION DEVELOPED FOR THIS SCRIPT MOVED TO PHYLOSTATS.JL

### Count one deletion per removed branch.  Within a branch the name with the
### LARGEST full definition is the clade itself; narrower names are aliases that
### coincided with it in this tree.  (The TRIGGER uses the opposite rule — the
### tightest name — because that is a question about calibration, not naming.)
function tally_clade_deletions!(rep_counts::Dict{String,Int64}, alias_counts::Dict{String,Int64},
                                names_by_branch::Dict{String,Set{String}},
                                clade_definition_sizes::Dict{String,Int64})::Nothing
    for (_, names_on_branch) in names_by_branch
        names_sharing_branch::Vector{String} = sort(collect(names_on_branch))
        representative::String = sort(names_sharing_branch, by = n -> (-get(clade_definition_sizes, n, 0), n))[1]
        rep_counts[representative] = get(rep_counts, representative, 0) + 1
        for name::String in names_sharing_branch
            if name != representative
                alias_counts[name] = get(alias_counts, name, 0) + 1
            end
        end
    end
    return nothing
end

### Write one clade deletion table and its stacked bar figure.
function write_clade_deletion_report(rep_counts::Dict{String,Int64}, alias_counts::Dict{String,Int64},
                                     tsv_path::String, png_path::String, plot_title::String,
                                     scope_description::String)::Nothing
    clade_labels::Vector{String} = sort(collect(union(keys(rep_counts), keys(alias_counts))))

    if isempty(clade_labels)
        println("\nNo named clade was removed by $(scope_description) in any tree — $(png_path) not written.")
        return nothing
    end

    open(tsv_path, "w") do fh
        println(fh, "clade\tdeleted_as_clade\tdeleted_as_alias\ttotal")
        for label::String in clade_labels
            n_rep::Int64 = get(rep_counts, label, 0)
            n_alias::Int64 = get(alias_counts, label, 0)
            println(fh, "$(label)\t$(n_rep)\t$(n_alias)\t$(n_rep + n_alias)")
        end
    end
    println("\nWrote: $(tsv_path)")

    rep_values::Vector{Int64} = [get(rep_counts, l, 0) for l in clade_labels]
    alias_values::Vector{Int64} = [get(alias_counts, l, 0) for l in clade_labels]
    positions::Vector{Int64} = collect(1:length(clade_labels))

    fig_c = Figure(size = (max(800, 45 * length(clade_labels)), 600))
    ax_c = Axis(fig_c[1,1],
                xlabel = "Clade",
                ylabel = "Number of trees in which the clade was deleted",
                title = plot_title,
                xticklabelrotation = π/4)
    barplot!(ax_c,
             vcat(positions, positions),
             vcat(rep_values, alias_values),
             stack = vcat(fill(1, length(positions)), fill(2, length(positions))),
             color = vcat(fill(:steelblue, length(positions)), fill(:orange, length(positions))))
    ax_c.xticks = (positions, clade_labels)
    Legend(fig_c[1,2],
           [PolyElement(color = :steelblue), PolyElement(color = :orange)],
           ["Clade itself", "Alias — shared a branch\nwith another clade"])
    save(png_path, fig_c)
    println("Wrote: $(png_path)")
    return nothing
end

     
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
        help = "file with clades to test: one line per clade. Name of clade first then all taxa in clade, all separated by spaces (no spaces allowed inside names of clades or species: e.g. Clade_A NOT Clade A, and Taxon_1 NOT Taxon 1). REQUIRED: the named-clade stem rule runs on every tree, so this file must always be supplied."
        required = true
        "--excluded_from_quartet_tests", "-x"
        help = "file with clades expected to be naturally long branch — excluded from quartet tests"
        required = false
        "--safeguard__quantiles", "-s"
        help = "file with the safeguard quantiles: ONE row of THREE numbers separated by whitespace, e.g. '0.95 0.95 0.95'. In order: (1) terminal taxa AND quartet rule — one value serves both, a branch must exceed this quantile of ALL terminal branch lengths before either rule can delete the taxon; (2) internal branches / long subtrees — the branch subtending the clade to remove must exceed this quantile of ALL internal branch lengths; (3) named clades — the stem must exceed this quantile of ALL clade stem lengths. A rule deletes only when its K trigger AND its safeguard are both exceeded."
        required = true
        "--threshold__triggers", "-t"
        help = "file with the K trigger values: ONE row of FOUR numbers separated by whitespace, e.g. '3 3 3 3' (fractional values such as 3.45 are allowed). In order: (1) terminal taxa, (2) internal branches / long subtrees, (3) named clade stems, (4) quartet rule. For 1-3 the trigger is (K-1)*local_median + global_median; for 4 it is the ratio of the focal taxon's distance from the anchor to the median distance of its comparators. Suggested starting value 3 throughout — use long_branches_statistical_analyses.jl to choose."
        required = true
        "--internal_side_ratio", "-r"
        help = "asymmetry threshold for the long-subtree rule. The two sides of a long internal branch are compared by their MEDIAN distance from the branch out to their own tips; the side removed must be the internally longer one AND the smaller one, and the two medians must differ by more than this ratio. So 1.2 means one side must be at least 20 percent deeper than the other before the branch is taken to say which side is the problem. Stops a small clade of short-branched taxa being deleted just because it sits beyond a long stem - there the long branch is the stem, not the clade. No default: set it deliberately. Must be greater than 1.0 (1.0 would accept any difference at all and effectively disable the test). Lower values (1.1) are permissive, higher values (1.3, 1.5) demand a starker contrast and remove fewer clades."
        required = true
        arg_type = Float64
        "--internal_long_path__trigger", "-p"
        help = "how many long internal branches, stacked one on the other, are needed before the path counts as a long-branch path and the clade beyond it is removed. No default: this must be set deliberately. The minimum accepted is 1, but 1 means ANY single long internal branch triggers removal, which is very liberal — 2 or 3 is normally more sensible."
        required = true
        arg_type = Int
    end
    return parse_args(s)
end

main()
