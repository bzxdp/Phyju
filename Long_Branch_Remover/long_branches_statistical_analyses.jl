using TreesIO
using TreeUtils
using TreeStats
using NewickTree
using MetaGraphsNext
using ArgParse
using Statistics
using Random
using CairoMakie

### The long-branch rules themselves live next door, not in the libraries: they are
### specific to this analysis (K triggers, safeguard quantiles, clade files) and were
### cluttering TreeStats/TreeUtils, which are shared with other projects.
include("long_branch_rules.jl")


function main()
    args= parse_arguments()
    treefile_extension::String= args["extension"]
    clades_file::String= args["clades"]
    quartet_tests_excluded::Union{String,Nothing}= args["excluded_from_quartet_tests"] ## list of clades not in local quartet tests

    tresh_triggers::String= args["threshold__triggers"]
    
    list_of_treefiles::Vector{String}= filter(f->endswith(f, treefile_extension), readdir("."))
    if isempty(list_of_treefiles)
        println("No files found with extension $treefile_extension")
        exit(1)
    end

    trees::Dict{String,MetaGraph}=Dict{String,MetaGraph}() ## data structure keeping all the trees in memory as a dictionary of metagraphs
    
    total_files::Int64 = length(list_of_treefiles)
    for (i, file::String) in enumerate(list_of_treefiles)
        println("Processing tree $i of $total_files: $file")
        trees[file]= readtree(file, nothing, "yes")
    end
    tree_labels= collect(keys(trees))

    ### Read file and store taxa to exclude from quartet tests                                                                                
    masked_clades::Dict{String,Vector{String}}= Dict{String,Vector{String}}()
    if !isnothing(args["excluded_from_quartet_tests"])
        open(quartet_tests_excluded, "r") do fh
            for clade::String in eachline(fh)
                if isempty(strip(clade))
                    continue
                end
                clade_definition::Vector{String}= split(clade, " ")
                clade_name::String= popfirst!(clade_definition)
                masked_clades[clade_name]=clade_definition
            end
        end
    end

    ### Triggering thresholds
    triggering_values::Dict{String,Vector{Float64}}=Dict{String,Vector{Float64}}()
    open(tresh_triggers, "r") do fh
        for line::String in eachline(fh)
            triggers::Vector{String}= split(strip(line))
            method::String = popfirst!(triggers)
            triggering_values[method]= parse.(Float64, triggers)
        end
    end
            
    ##### stats for terminal taxa, clades, internal branches and quartets
    (global_bls_all_taxa::Vector{Float64},taxon_specific_bls::Dict{String,Vector{Float64}},taxon_specific_triggering_tresholds::Dict{String,Float64})= stats_for_terminals(trees, triggering_values["terminals"][1])
    (bls_and_their_lenght_by_tree::Dict{String,Dict{Tuple{String,String},Float64}},global_set_internal_blens::Vector{Float64},triggering_values_per_tree::Dict{String,Float64})= stats_for_internal_branches(trees, triggering_values["internals"][1])
    ## comparator floor scaled to the dataset — see comparator_branch_threshold()
    comparator_branch_floor::Float64 = comparator_branch_threshold(global_bls_all_taxa)
    println("\nQuartet comparators must sit on a terminal branch longer than $(comparator_branch_floor)")
    (global_set_of_ratios::Vector{Float64},all_ratios_over_all_taxa::Dict{String,Vector{Float64}})= stats_for_quartets(trees,masked_clades,triggering_values["quartets"][1]; min_comparator_branch = comparator_branch_floor)
    (global_bls_all_stems::Vector{Float64},stem_specific_bls::Dict{String,Vector{Float64}},stem_specific_triggering_tresholds::Dict{String,Float64})= stats_for_stems(trees,clades_file,triggering_values["stems"][1])

    ##### Alternative triggers for plotting
    ## These MUST be computed the same way TreeStats computes the global term of the
    ## real triggers, otherwise the alternative-K lines drawn below would not
    ## correspond to what long_branches_identifier.jl would actually do.
    ## In every case: MEDIAN OF THE PER-UNIT MEDIANS, one vote per unit —
    ##   terminals -> one vote per taxon, clades -> one vote per clade,
    ##   internals -> one vote per gene tree.
    ## NOT the median of the pooled branch lengths, which would be weighted by gene
    ## occupancy (terminals, clades) or by tree size (internals).
    median_clades_bl::Float64= median([median(v) for v in values(stem_specific_bls) if !isempty(v)])
    median_terminals_bl::Float64= median([median(v) for v in values(taxon_specific_bls) if !isempty(v)])
    median_internals_bl::Float64= median([median(collect(values(d))) for d in values(bls_and_their_lenght_by_tree) if !isempty(d)])

    #### NOTE for quartets the triggers do not need to be defined they are just the raw values plotted. So I am not defining a clade of alternative triggers as done below for stem terminals and internals.
    
    ### ALTERNTIVE STEM K FOR PLOTTING
    stem_specific_alternative_triggering_tresholds::Dict{String,Vector{Float64}}= Dict{String,Vector{Float64}}()
    clades_names::Vector{String}= collect(keys(stem_specific_bls))

    for clade::String in clades_names
        if isempty(stem_specific_bls[clade])
            @warn "ATTENTION Clade $(clade) has NO BLs"
            continue
        end
        stem_specific_alternative_triggering_tresholds[clade] = Vector{Float64}()
        median_current_clade_bl::Float64=median(stem_specific_bls[clade])
        counter::Int64= 2
        while counter <= length(triggering_values["stems"])
            alternative_trigger::Float64 = triggering_values["stems"][counter]
            current_trigger::Float64= (alternative_trigger - 1) * median_current_clade_bl + median_clades_bl
            push!(stem_specific_alternative_triggering_tresholds[clade], current_trigger)
            counter += 1
        end
    end

    ### ALTERNATIVE K FOR TERMINALS
    terminal_taxa_specific_alternative_triggering_tresholds::Dict{String,Vector{Float64}}= Dict{String,Vector{Float64}}()
    taxon_names::Vector{String}= collect(keys(taxon_specific_bls))
    
    for	taxon::String in taxon_names
	terminal_taxa_specific_alternative_triggering_tresholds[taxon] = Vector{Float64}()
        median_current_taxon_bl::Float64=median(taxon_specific_bls[taxon])
	counter::Int64= 2
	while counter <= length(triggering_values["terminals"])
            alternative_trigger::Float64 = triggering_values["terminals"][counter]
            current_trigger::Float64= (alternative_trigger - 1) * median_current_taxon_bl + median_terminals_bl
            push!(terminal_taxa_specific_alternative_triggering_tresholds[taxon], current_trigger)
            counter += 1
        end
    end

    ### ALTERNATIVE K FOR INTERNALS
    tree_specific_alternative_triggering_tresholds::Dict{String,Vector{Float64}}= Dict{String,Vector{Float64}}()
    
    for tree::String in tree_labels
        tree_specific_alternative_triggering_tresholds[tree] = Vector{Float64}()
        current_tree_bls_vals::Vector{Float64} = collect(values(bls_and_their_lenght_by_tree[tree]))
        median_current_tree_bls::Float64=median(current_tree_bls_vals)
        counter::Int64= 2
	while counter <= length(triggering_values["internals"])
            alternative_trigger::Float64 = triggering_values["internals"][counter]
            current_trigger::Float64= (alternative_trigger - 1) * median_current_tree_bls + median_internals_bl
            push!(tree_specific_alternative_triggering_tresholds[tree], current_trigger)
            counter += 1
        end
    end

    #### compute quantiles for plotting
    quantile_vals::Vector{Float64}= [0.25, 0.5, 0.75]
    
    quantiles_taxa::Vector{Float64}=Vector{Float64}()
    quantiles_stems::Vector{Float64}=Vector{Float64}()
    quantiles_internal::Vector{Float64}=Vector{Float64}()
    quantiles_quartets::Vector{Float64}=Vector{Float64}()
    
    for quant in quantile_vals
        push!(quantiles_taxa, quantile(collect(global_bls_all_taxa), quant))
        push!(quantiles_stems, quantile(collect(global_bls_all_stems), quant))
        push!(quantiles_internal, quantile(collect(global_set_internal_blens), quant))
        push!(quantiles_quartets, quantile(collect(global_set_of_ratios), quant))
    end
    
    ## global safeguard percentiles
    p95_internal::Float64 = quantile(collect(global_set_internal_blens), 0.95)
    p99_internal::Float64 = quantile(collect(global_set_internal_blens), 0.99)
    p95_taxa::Float64    = quantile(collect(global_bls_all_taxa), 0.95)
    p99_taxa::Float64    = quantile(collect(global_bls_all_taxa), 0.99)
    p95_stems::Float64   = quantile(collect(global_bls_all_stems), 0.95)
    p99_stems::Float64   = quantile(collect(global_bls_all_stems), 0.99)
    p95_quartets::Float64   = quantile(collect(global_set_of_ratios), 0.95)
    p99_quartets::Float64   = quantile(collect(global_set_of_ratios), 0.99)
    
    k_colours = [:lightsalmon, :salmon, :red, :darkred] ### colours for plotting alternative K valules (triggers)
    
    taxa_in_trees::Vector{String} = collect(keys(all_ratios_over_all_taxa))
    clades_in_trees::Vector{String} = collect(keys(stem_specific_bls))
    
    all_q_ratios::Vector{Float64} = Vector{Float64}()
    all_q_taxa::Vector{String} = Vector{String}()
    for taxon::String in taxa_in_trees
        for qr::Float64 in all_ratios_over_all_taxa[taxon]
            push!(all_q_taxa, taxon)
            push!(all_q_ratios, qr)
        end
    end
    
    all_clade_bls::Vector{Float64} = Vector{Float64}()
    for clade::String in clades_in_trees
        for bl::Float64 in stem_specific_bls[clade]
            push!(all_clade_bls, bl)
        end
    end
    
    all_terminal_bls::Vector{Float64} = Vector{Float64}()
    all_taxa::Vector{String} = Vector{String}()
    for taxon::String in taxa_in_trees
        for bl::Float64 in taxon_specific_bls[taxon]
            push!(all_taxa, taxon)
            push!(all_terminal_bls, bl)
        end
    end
    
    #### plot distribution quartet ratios taxa
    ### this is a debug sequence.
    println("Quartet ratio range: min=$(minimum(all_q_ratios)) max=$(maximum(all_q_ratios)) n=$(length(all_q_ratios))")
    println("Number of taxa: $(length(taxa_in_trees))")
    println("Quartet ratio 99th percentile: $(quantile(all_q_ratios, 0.99))")
    println("Quartet ratio 95th percentile: $(quantile(all_q_ratios, 0.95))")
    println("Quartet ratio median: $(median(all_q_ratios))")

    fig_q = Figure(size = (max(800, 75 * length(taxa_in_trees)), 600))
    ax_q = Axis(fig_q[1,1],
            ylabel = "Quartet Ratio Values",
            xlabel = "Taxa",
            title = "Terminal taxa quartet ratio distributions",
            xticklabelrotation = π/4,
            limits = (nothing, (0, quantile(all_q_ratios, 0.99) * 1.5))
            )
    for (i, taxon) in enumerate(taxa_in_trees)
    ratios = all_ratios_over_all_taxa[taxon]
    if !isempty(ratios)
        violin!(ax_q, fill(i, length(ratios)), ratios, width = 0.8, color = :steelblue)
    end
    ## alternative K lines
    for (ki, k_val) in enumerate(triggering_values["quartets"])
        scatter!(ax_q, [i], [Float64(k_val)], color = k_colours[min(ki, length(k_colours))], markersize = 6)
    end
end
    ax_q.xticks = (1:length(taxa_in_trees), taxa_in_trees)
    for quant_val in quantiles_quartets
    hlines!(ax_q, [quant_val], color = :grey, linestyle = :dot, linewidth = 1)
end
    hlines!(ax_q, [p95_quartets], color = :green, linestyle = :dot, linewidth = 2)
    hlines!(ax_q, [p99_quartets], color = :orange, linestyle = :dot, linewidth = 2)
    save("terminal_taxa_quartet_ratio_distributions.png", fig_q)
    
    #### plot global distribution internal bls
    fig_i = Figure(size = (800, 600))
    ax_i = Axis(fig_i[1,1],
            ylabel = "Frequency",
            xlabel = "Internal branch length",
            title = "Global distribution of internal branch lengths",
            limits = (nothing, (0, maximum(collect(global_set_internal_blens)) * 1.05))
            )
    hist!(ax_i, collect(global_set_internal_blens), bins = 10, color = :steelblue)
    for value in quantiles_internal
    vlines!(ax_i, [value], color = :grey, linestyle = :dot, linewidth = 1)
end
    vlines!(ax_i, [p95_internal], color = :green, linestyle = :dot, linewidth = 1)
    vlines!(ax_i, [p99_internal], color = :orange, linestyle = :dot, linewidth = 1)
    save("internal_bls_distribution.png", fig_i)
    
    #### global violin for internal branches
    fig_iv = Figure(size = (800, 600))
    ax_iv = Axis(fig_iv[1,1],
             ylabel = "Branch length",
             title = "Internal BLs: global distribution",
             limits = (nothing, (0, maximum(collect(global_set_internal_blens)) * 1.05))
             )
    violin!(ax_iv, fill(1, length(global_set_internal_blens)), collect(global_set_internal_blens),
        color = :steelblue, width = 0.8)
    ax_iv.xticks = ([1], ["internal branches"])
    hlines!(ax_iv, [p95_internal], color = :green, linestyle = :dot, linewidth = 2)
    hlines!(ax_iv, [p99_internal], color = :orange, linestyle = :dot, linewidth = 2)
    for quant_int in quantiles_internal
    hlines!(ax_iv, [quant_int], color = :grey, linestyle = :dot, linewidth = 1)
end
    scatter!(ax_iv, [1], [median(collect(global_set_internal_blens))], marker = :diamond, color = :black, markersize = 12)
    save("internal_bls_global_violin.png", fig_iv)
    
    #### per-tree violins for internal branches
    for tree::String in tree_labels
    if !haskey(triggering_values_per_tree, tree)
        continue
    end
    tree_bls::Vector{Float64} = collect(values(bls_and_their_lenght_by_tree[tree]))
    if isempty(tree_bls)
        continue
    end
    tree_trigger::Float64 = triggering_values_per_tree[tree]
    tree_median::Float64 = median(tree_bls)
    
    fig_t = Figure(size = (800, 600))
    ax_t = Axis(fig_t[1,1],
                ylabel = "Internal branch length",
                title = "Internal BLs: $(tree)\nmedian=$(round(tree_median, digits=4)) trigger=$(round(tree_trigger, digits=4))",
                limits = (nothing, (0, maximum(tree_bls) * 1.05))
                )
    violin!(ax_t, fill(1, length(tree_bls)), tree_bls, color = :steelblue, width = 0.8)
    ax_t.xticks = ([1], ["internal branches"])
    
    ## alternative K lines
    if haskey(tree_specific_alternative_triggering_tresholds, tree)
        for (ki, trigger) in enumerate(tree_specific_alternative_triggering_tresholds[tree])
            hlines!(ax_t, [trigger], color = k_colours[min(ki, length(k_colours))], linestyle = :dash, linewidth = 1)
        end
    end
    hlines!(ax_t, [tree_trigger], color = :red, linestyle = :dot, linewidth = 2)
    hlines!(ax_t, [p95_internal], color = :green, linestyle = :dot, linewidth = 1)
    hlines!(ax_t, [p99_internal], color = :orange, linestyle = :dot, linewidth = 1)
    scatter!(ax_t, [1], [tree_median], marker = :diamond, color = :black, markersize = 12)
    scatter!(ax_t, [1], [tree_trigger], marker = :circle, color = :red, markersize = 10)
    save("violin_$(tree).png", fig_t)
end
    
    #### plot distribution bl clades
    fig_s = Figure(size = (max(800, 75 * length(clades_in_trees)), 600))
    ax_s = Axis(fig_s[1,1],
            ylabel = "Branch lengths",
            xlabel = "Clade",
            title = "Stem clades branch lengths distributions",
            xticklabelrotation = π/4,
            limits = (nothing, (0, maximum(all_clade_bls) * 1.05))
            )
    for (i, clade) in enumerate(clades_in_trees)
    bls = stem_specific_bls[clade]
    if !isempty(bls)
        violin!(ax_s, fill(i, length(bls)), bls, width = 0.8, color = :steelblue)
    end
    ## default K trigger
    if haskey(stem_specific_triggering_tresholds, clade)
        scatter!(ax_s, [i], [stem_specific_triggering_tresholds[clade]], color = :red, markersize = 8)
    end
    ## alternative K triggers
    if haskey(stem_specific_alternative_triggering_tresholds, clade)
        for (ki, trigger) in enumerate(stem_specific_alternative_triggering_tresholds[clade])
            scatter!(ax_s, [i], [trigger], color = k_colours[min(ki, length(k_colours))], markersize = 6)
        end
    end
end
    ax_s.xticks = (1:length(clades_in_trees), clades_in_trees)
    for quant_val in quantiles_stems
    hlines!(ax_s, [quant_val], color = :grey, linestyle = :dot, linewidth = 1)
end
    hlines!(ax_s, [p95_stems], color = :green, linestyle = :dot, linewidth = 2)
    hlines!(ax_s, [p99_stems], color = :orange, linestyle = :dot, linewidth = 2)
    save("clade_stems_lengths_distribution.png", fig_s)
    
    #### plot distribution bl taxa
    fig_t2 = Figure(size = (max(800, 75 * length(taxa_in_trees)), 600))
    ax_t2 = Axis(fig_t2[1,1],
             ylabel = "Branch lengths",
             xlabel = "Taxa",
             title = "Terminal taxa branch lengths distributions",
             xticklabelrotation = π/4,
             limits = (nothing, (0, maximum(all_terminal_bls) * 1.05))
             )
    for (i, taxon) in enumerate(taxa_in_trees)
    bls = taxon_specific_bls[taxon]
    if !isempty(bls)
        violin!(ax_t2, fill(i, length(bls)), bls, width = 0.8, color = :steelblue)
    end
    ## default K trigger
    if haskey(taxon_specific_triggering_tresholds, taxon)
        scatter!(ax_t2, [i], [taxon_specific_triggering_tresholds[taxon]], color = :red, markersize = 8)
    end
    ## alternative K triggers
    if haskey(terminal_taxa_specific_alternative_triggering_tresholds, taxon)
        for (ki, trigger) in enumerate(terminal_taxa_specific_alternative_triggering_tresholds[taxon])
            scatter!(ax_t2, [i], [trigger], color = k_colours[min(ki, length(k_colours))], markersize = 6)
        end
    end
end
    ax_t2.xticks = (1:length(taxa_in_trees), taxa_in_trees)
    for quant_val in quantiles_taxa
    hlines!(ax_t2, [quant_val], color = :grey, linestyle = :dot, linewidth = 1)
end
    hlines!(ax_t2, [p95_taxa], color = :green, linestyle = :dot, linewidth = 2)
    hlines!(ax_t2, [p99_taxa], color = :orange, linestyle = :dot, linewidth = 2)
    save("terminal_taxa_lengths_distribution.png", fig_t2)

end
### FUNCTIONS
### All functions developed for this script moved to PhyloStats.jl

function parse_arguments()
    s = ArgParseSettings(description="long branch statistics and detection")
    @add_arg_table s begin
        "--extension", "-e"
            help = "extension of the tree files e.g. treefile"
            required = true
        "--clades", "-c"
            help = "file with clade definitions: CladeName taxon1 taxon2 ..."
            required = true
        "--excluded_from_quartet_tests", "-x"
            help = "file with clades expected to be naturally long branch — excluded from quartet tests"
            required = false
        "--threshold__triggers", "-t"
            help = "file with trigger value (K) used to calculate when a branch is long (suggested default value 3). This file has four lines, each start with the name of one of the method followed by the default value (3 is suggested) and then a series of alternative K values to compare against the derfault (e.g. terminals 3 1 2 4 5)."  
            required = true
    end
    return parse_args(s)
end

main()
