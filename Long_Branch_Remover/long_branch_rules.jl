### =====================================================================
### LONG BRANCH REMOVER — RULES
###
### Everything here exists to serve long_branches_identifier.jl and
### long_branches_statistical_analyses.jl.  It was moved out of TreeStats and
### TreeUtils so those libraries stay general: nothing below is a plain tree
### operation, every function takes a K trigger, a safeguard quantile, a clade
### definition file or some other long-branch policy.
###
### The code and the comments are exactly as they were in the libraries — only
### the location changed.  Section banners like this one are the only additions.
###
### Both scripts pull it in with:   include("long_branch_rules.jl")
### and must already have: using TreesIO, TreeUtils, MetaGraphsNext, Statistics
### =====================================================================

## Minimum number of trees in which a taxon (or a clade stem) must be observed
## before its own median is trusted enough to test against.  Mirrors Python
## MIN_CROSS_GENE_OBS.  Below this the rule issues no trigger at all, so the
## taxon or clade is RETAINED — a median built on one or two values is noise.
const MIN_CROSS_GENE_OBS = 3

## Shortest TERMINAL BRANCH a tip may have and still serve as a comparator.
## Tree inference pins branches it cannot estimate at a minimum value (IQ-TREE uses
## 1e-6): such a branch is UNESTIMATED, not short, and a tip carrying one tells you
## nothing about the local scale.  Comparing a taxon against it is comparing against
## missing data, and the ratio explodes.
##
## The test is on the tip's OWN terminal branch, not on its distance from the anchor.
## That distinction matters: in a cluster of near-identical sequences every tip sits
## on a floor-length branch, but two or three edges out from the anchor the ACCUMULATED
## distance already exceeds 1e-6, so a distance test keeps precisely the tips that
## should be excluded.  Measured on real data, the distance form left the maximum
## quartet ratio at 20805 and made matters worse by pushing anchors outward; this form
## does not have that failure.
##
## The Python prototype has no equivalent guard at all and produces the same ratios in
## the thousands; it never notices because it never pools them.
const MIN_COMPARATOR_BRANCH = 1e-6

## The floor above is an absolute backstop only.  What actually counts as "no
## information" depends on the dataset: measured on real trees the noise branches were
## 1e-6, 2.06e-6, 2.09e-6, 2.23e-6 and 2.54e-6 while the real branches around them were
## 0.0089 and up - three orders of magnitude clear.  A fixed 1e-6 cutoff catches only
## the exactly-floored ones and lets 2.5e-6 through, which is enough to drag a local
## median to 7e-6 and produce a ratio of 1280 for a taxon whose true ratio is 0.52.
## So the working threshold is a fraction of the dataset's typical terminal branch.
const COMPARATOR_SCALE_FRACTION = 0.01

## threshold = max(absolute floor, fraction of the typical terminal branch)
function comparator_branch_threshold(global_terminal_bls::Vector{Float64})::Float64
    if isempty(global_terminal_bls)
        return MIN_COMPARATOR_BRANCH
    end
    return max(MIN_COMPARATOR_BRANCH, COMPARATOR_SCALE_FRACTION * median(global_terminal_bls))
end


### ---------------------------------------------------------------------
### moved from TreeUtils: these take a trigger, or operate on a "stack" of long
### branches, so they are rules rather than tree utilities.
### identify_clade_for_removal currently has NO caller — it was superseded by
### clade_for_removal_by_asymmetry.  Kept here rather than deleted.
### ---------------------------------------------------------------------

function find_long_branch_stacks(unrooted_tree::MetaGraph, bls::Dict{Tuple{String,String},Float64},trigger::Float64,min_stack_length::Int64)::Vector{Vector{Tuple{String,String}}}

    canonical(u, v) = u < v ? (u, v) : (v, u)

    ## find all long internal edges
    canonical_long::Set{Tuple{String,String}} = Set{Tuple{String,String}}()
    for (edge, bl) in bls
        (u, v) = edge
        if !unrooted_tree[u].is_a_leaf && !unrooted_tree[v].is_a_leaf && bl > trigger
            push!(canonical_long, canonical(u, v))
        end
    end

    ## find connected components of long internal edges
    ## two long internal edges are stacked if they share an internal node
    visited_edges::Set{Tuple{String,String}} = Set{Tuple{String,String}}()
    stacks::Vector{Vector{Tuple{String,String}}} = Vector{Vector{Tuple{String,String}}}()

    for edge in canonical_long
        if edge in visited_edges
            continue
        end
        stack::Vector{Tuple{String,String}} = Vector{Tuple{String,String}}()
        queue::Vector{Tuple{String,String}} = [edge]
        while !isempty(queue)
            current_edge = popfirst!(queue)
            if current_edge in visited_edges
                continue
            end
            push!(visited_edges, current_edge)
            push!(stack, current_edge)
            (u, v) = current_edge
            for node in [u, v]
                for nb in neighbor_labels(unrooted_tree, node)
                    if !unrooted_tree[nb].is_a_leaf
                        candidate = canonical(node, nb)
                        if candidate in canonical_long && !(candidate in visited_edges)
                            push!(queue, candidate)
                        end
                    end
                end
            end
        end
        if length(stack) >= min_stack_length
            push!(stacks, stack)
        end
    end

    return stacks
end


function get_outermost_edge(unrooted_tree::MetaGraph,stack::Vector{Tuple{String,String}})::Tuple{String,String}
    best_edge = stack[1]
    best_larger_side = 0
    n_leaves = length(get_leaves(unrooted_tree))

    for (u, v) in stack
        leaves_v::Set{String} = Set{String}()
        visited::Set{String} = Set{String}([u])
        queue::Vector{String} = [v]
        while !isempty(queue)
            current = popfirst!(queue)
            if current in visited; continue; end
            push!(visited, current)
            if unrooted_tree[current].is_a_leaf
                push!(leaves_v, current)
            end
            for nb in neighbor_labels(unrooted_tree, current)
                if !(nb in visited); push!(queue, nb); end
            end
        end
        larger_side = max(length(leaves_v), n_leaves - length(leaves_v))
        if larger_side > best_larger_side
            best_larger_side = larger_side
            best_edge = (u, v)
        end
    end
    return best_edge
end


function identify_clade_for_removal(unrooted_tree::MetaGraph,long_branch_edge::Tuple{String,String})::Vector{String}

    (u, v) = long_branch_edge

    ## BFS from v excluding u
    leaves_v::Set{String} = Set{String}()
    total_bl_v::Float64 = 0.0
    visited_v::Set{String} = Set{String}([u])
    queue_v::Vector{String} = [v]
    while !isempty(queue_v)
        current = popfirst!(queue_v)
        if current in visited_v; continue; end
        push!(visited_v, current)
        if unrooted_tree[current].is_a_leaf
            push!(leaves_v, current)
        end
        for nb in neighbor_labels(unrooted_tree, current)
            if !(nb in visited_v)
                edge_data = unrooted_tree[current, nb]
                if !ismissing(edge_data.length)
                    total_bl_v += edge_data.length
                end
                push!(queue_v, nb)
            end
        end
    end

    ## BFS from u excluding v
    leaves_u::Set{String} = Set{String}()
    total_bl_u::Float64 = 0.0
    visited_u::Set{String} = Set{String}([v])
    queue_u::Vector{String} = [u]
    while !isempty(queue_u)
        current = popfirst!(queue_u)
        if current in visited_u; continue; end
        push!(visited_u, current)
        if unrooted_tree[current].is_a_leaf
            push!(leaves_u, current)
        end
        for nb in neighbor_labels(unrooted_tree, current)
            if !(nb in visited_u)
                edge_data = unrooted_tree[current, nb]
                if !ismissing(edge_data.length)
                    total_bl_u += edge_data.length
                end
                push!(queue_u, nb)
            end
        end
    end

    ## smaller side goes — tiebreak on higher total BL
    if length(leaves_v) < length(leaves_u)
        return collect(leaves_v)
    elseif length(leaves_u) < length(leaves_v)
        return collect(leaves_u)
    else
        return total_bl_v >= total_bl_u ? collect(leaves_v) : collect(leaves_u)
    end
end

### ---------------------------------------------------------------------
### moved from TreeStats: cross-gene backgrounds, the four rules, and the
### quartet machinery.
### ---------------------------------------------------------------------

function stats_for_quartets(trees::Dict{String,MetaGraph},masked_clades::Dict{String,Vector{String}},ratio_threshold::Float64; min_comparator_branch::Float64 = MIN_COMPARATOR_BRANCH)::Tuple{Vector{Float64},Dict{String,Vector{Float64}}}

    tree_labels::Vector{String} = collect(keys(trees))
    total_set_of_taxa::Set{String} = Set{String}()
    
        
    ## returned values
    all_ratios_over_all_taxa::Dict{String,Vector{Float64}} = Dict{String,Vector{Float64}}()
    ## pooled list of EVERY quartet ratio observed — duplicates kept, so medians
    ## and quantiles describe the empirical distribution (Python keeps a list too)
    global_set_of_ratios::Vector{Float64} = Vector{Float64}()

    for tree::String in tree_labels
        leaves_current_tree::Set{String} = get_leaves(trees[tree])
        total_set_of_taxa = union!(total_set_of_taxa, leaves_current_tree)
    end

    total_taxa::Int64 = length(total_set_of_taxa)

    ## mask fragments are tree-specific — derive them ONCE PER TREE up front,
    ## otherwise the fragment search would be repeated for every taxon.
    mask_maps::Dict{String,Dict{String,String}} = Dict{String,Dict{String,String}}()
    for tree::String in tree_labels
        mask_maps[tree] = build_quartet_mask_map(trees[tree], masked_clades)
    end

    for (i, taxon::String) in enumerate(total_set_of_taxa)
        println("Computing quartet stats: taxon $i of $total_taxa")
        all_ratios_of_taxon::Vector{Float64} = Vector{Float64}()
        for tree::String in tree_labels
            if !haskey(trees[tree], taxon)
                continue
            end
            if !trees[tree][taxon].is_a_leaf
                continue
            end
            ## skip taxa with near-zero branch lengths
            ## A FOCAL taxon on a branch the inference could not estimate has no ratio
            ## worth recording either: the numerator is noise rather than the
            ## denominator.  Same scaled threshold as the comparators, so both ends of
            ## the distribution are cleaned on the same criterion.
            ## This affects the plotted distribution ONLY.  The identifier does not skip
            ## such taxa - there the absolute branch-length safeguard already makes them
            ## unremovable, which is the better mechanism because it is data-driven.
            neighbour::String = first(neighbor_labels(trees[tree], taxon))
            edge::EdgeData = trees[tree][taxon, neighbour]
            if ismissing(edge.length) || edge.length <= min_comparator_branch
                continue
            end
            (is_long, local_ratio) = quartet_test_for_taxon(trees[tree], taxon, mask_maps[tree], ratio_threshold; min_comparator_branch = min_comparator_branch)
            if !isnan(local_ratio)
                push!(all_ratios_of_taxon, local_ratio)
            end
        end
        all_ratios_over_all_taxa[taxon] = all_ratios_of_taxon
    end
    for taxon::String in total_set_of_taxa
        append!(global_set_of_ratios, all_ratios_over_all_taxa[taxon])
    end
    return global_set_of_ratios, all_ratios_over_all_taxa
end
    
    

function stats_for_terminals(trees::Dict{String,MetaGraph}, threshold::Float64)::Tuple{Vector{Float64},Dict{String,Vector{Float64}},Dict{String,Float64}}
    tree_labels::Vector{String}= collect(keys(trees))
    total_set_of_taxa::Set{String}=Set{String}()
    median_bls_for_all_individual_taxa::Dict{String,Float64}= Dict{String,Float64}()
    global_median_bl::Float64=0.0
    all_bl_over_all_taxa::Dict{String,Vector{Float64}}= Dict{String,Vector{Float64}}()
    table_of_triggering_tresholds::Dict{String,Float64}= Dict{String,Float64}()
    ## pooled list of EVERY terminal branch length — duplicates kept
    global_set_of_bls::Vector{Float64}=Vector{Float64}()
    
    for tree::String in tree_labels
        leaves_current_tree::Set{String}= get_leaves(trees[tree])
        total_set_of_taxa= union!(total_set_of_taxa, leaves_current_tree)
    end

    for taxon::String in total_set_of_taxa
        all_bls_of_taxon::Vector{Float64}=Vector{Float64}()
        for tree::String in tree_labels
            if !haskey(trees[tree], taxon)
                continue
            end
            if !trees[tree][taxon].is_a_leaf
                continue
            else
                neighbour_of_leaf::String= first(neighbor_labels(trees[tree], taxon))
                edge::EdgeData= trees[tree][taxon, neighbour_of_leaf]
                bl::Float64 = edge.length
                push!(all_bls_of_taxon, bl)
            end
        end
        all_bl_over_all_taxa[taxon]=all_bls_of_taxon
    end

    for taxon::String in total_set_of_taxa
        if isempty(all_bl_over_all_taxa[taxon])
            continue
        end
        median_bl_taxon::Float64=median(all_bl_over_all_taxa[taxon])
        median_bls_for_all_individual_taxa[taxon]=median_bl_taxon
    end

    for taxon::String in total_set_of_taxa
        append!(global_set_of_bls, all_bl_over_all_taxa[taxon])
    end

    ## GLOBAL TERM = MEDIAN OF THE PER-TAXON MEDIANS, not the median of the pooled
    ## branch lengths.  The trigger is (K-1)*m_taxon + M, so when a taxon is exactly
    ## typical (m_taxon == M) the trigger reduces to K*M: the calibration only holds
    ## if M is the same KIND of quantity as m_taxon, i.e. a per-taxon median.
    ## Pooling would also weight the constant by gene occupancy, letting the
    ## best-sampled taxa set it.  Every taxon with at least one observation votes
    ## once (the MIN_CROSS_GENE_OBS gate governs testing, not this background).
    ## NOTE: global_set_of_bls is still returned POOLED — the safeguard quantile is
    ## a percentile of real branch lengths and must stay that way.
    global_median_bl= median(collect(values(median_bls_for_all_individual_taxa)))

    ## MINIMUM OBSERVATIONS GATE (Python MIN_CROSS_GENE_OBS).
    ## A taxon-specific median computed from one or two trees is not an estimate,
    ## it is noise, so no trigger is issued for such taxa and they are therefore
    ## never removed by this rule.  Note their branch lengths STILL enter
    ## global_set_of_bls: they inform the global background, they just cannot be
    ## judged against their own median.  In doubt, retain.
    insufficient_taxa::Vector{Tuple{String,Int64}} = Vector{Tuple{String,Int64}}()
    for taxon::String in total_set_of_taxa
        if !haskey(median_bls_for_all_individual_taxa, taxon)
            continue
        end
        n_observations::Int64 = length(all_bl_over_all_taxa[taxon])
        if n_observations < MIN_CROSS_GENE_OBS
            push!(insufficient_taxa, (taxon, n_observations))
            continue
        end
        triggering_threshold_taxon::Float64= (threshold - 1) * median_bls_for_all_individual_taxa[taxon] + global_median_bl
        table_of_triggering_tresholds[taxon]=triggering_threshold_taxon
    end

    if !isempty(insufficient_taxa)
        sort!(insufficient_taxa, by = x -> (x[2], x[1]))
        println("\nTerminal-branch rule not evaluable for $(length(insufficient_taxa)) taxa (present in fewer than $(MIN_CROSS_GENE_OBS) trees). These taxa are RETAINED by this rule:")
        for (taxon, n) in insufficient_taxa
            println("  $(taxon) ($(n) tree$(n == 1 ? "" : "s"))")
        end
    end

    return global_set_of_bls, all_bl_over_all_taxa, table_of_triggering_tresholds
end


function stats_for_stems(trees::Dict{String,MetaGraph},clade_file::String, threshold::Float64)::Tuple{Vector{Float64},Dict{String,Vector{Float64}},Dict{String,Float64}}
    tree_labels= collect(keys(trees))
    clade_definitions::Dict{String,Set{String}}= Dict{String,Set{String}}()
    median_bls_for_all_clades::Dict{String,Float64}= Dict{String,Float64}()
    global_median_stem_bl::Float64= 0.0
    table_of_triggering_tresholds_clades::Dict{String,Float64}= Dict{String,Float64}()
    ## pooled list of EVERY clade stem length — duplicates kept
    global_set_of_stem_bls::Vector{Float64}= Vector{Float64}()
    clades_and_their_stem_lenghts::Dict{String,Vector{Float64}}=Dict{String,Vector{Float64}}()

    open(clade_file, "r") do fh
        for clade::String in eachline(fh)
            clade_as_array::Vector{String}= split(strip(clade))
            if isempty(clade_as_array)
                continue
            end
            clade_name::String= popfirst!(clade_as_array)
            clade_definition::Set{String}= Set(String.(strip.(clade_as_array)))
            clade_definitions[clade_name]= clade_definition
        end
    end
    clade_names::Vector{String}= collect(keys(clade_definitions))

    for clade::String in clade_names
        clades_and_their_stem_lenghts[clade] = Vector{Float64}()
    end

    ## ALIASING.  Two clade names can resolve to the SAME physical branch in a given
    ## tree — e.g. CNID and CNID_A when only the CNID_A taxa are present, so CNID's
    ## largest monophyletic fragment IS the CNID_A fragment.  The loop is therefore
    ## tree-outer / clade-inner, so the resolved fragments can be compared within a
    ## tree and the shared branch handled once:
    ##   - every name that resolves to the branch records its length, so each clade's
    ##     own median is complete (CNID must not lose a tree just because CNID_A
    ##     described the same branch there);
    ##   - the branch enters the POOLED background only ONCE, so the global term is
    ##     not inflated by how many names happen to describe one branch.
    for tree::String in tree_labels
        leaves_current_tree::Set{String}= get_leaves(trees[tree])

        ## resolve every clade in this tree to (fragment, stem length)
        fragment_of_clade::Dict{String,Vector{String}} = Dict{String,Vector{String}}()
        stem_bl_of_clade::Dict{String,Float64} = Dict{String,Float64}()

        for clade::String in clade_names
            clade_composition_in_current_tree::Vector{String}= collect(intersect(leaves_current_tree, clade_definitions[clade]))

            if isempty(clade_composition_in_current_tree)
                continue
            end

            ## pass the clade name so any warning says WHICH clade; warn_on_failure is
            ## off because a missing LCA is the normal path into the fragment search
            ## below, and find_lca_largest_subset reports it with the name anyway
            lca_node::Union{String,Nothing}= find_lca(trees[tree], clade_composition_in_current_tree, clade; warn_on_failure=false)

            if isnothing(lca_node)
                (lca_node, largest_subset_of_members) = find_lca_largest_subset(trees[tree], clade_composition_in_current_tree, clade)
            else
                largest_subset_of_members = clade_composition_in_current_tree  ## full clade is monophyletic
            end

            if isnothing(lca_node)
                @warn "No monophyletic subset of size >= 2 found for clade $clade in tree $tree — all taxa scattered — skipping"
                continue
            end

            stem_edge::Union{Tuple{String,String},Nothing}= identify_stem_edge(trees[tree], lca_node, largest_subset_of_members)

            if isnothing(stem_edge)
                @warn "Could not identify stem edge for clade $clade in tree $tree — skipping"
                continue
            end

            (node_a, node_b)=stem_edge
            edge::EdgeData= trees[tree][node_a, node_b]
            stem_bl::Union{Float64,Missing} = edge.length

            if ismissing(stem_bl)
                continue
            end

            fragment_of_clade[clade] = sort(largest_subset_of_members)
            stem_bl_of_clade[clade]  = stem_bl
        end

        ## group the clade names that landed on the same fragment in THIS tree
        names_by_fragment::Dict{String,Vector{String}} = Dict{String,Vector{String}}()
        for clade::String in sort(collect(keys(fragment_of_clade)))
            fragment_key::String = join(fragment_of_clade[clade], "\u0001")
            push!(get!(names_by_fragment, fragment_key, Vector{String}()), clade)
        end

        for fragment_key::String in sort(collect(keys(names_by_fragment)))
            names_sharing_branch::Vector{String} = names_by_fragment[fragment_key]
            shared_bl::Float64 = stem_bl_of_clade[names_sharing_branch[1]]

            for clade::String in names_sharing_branch
                push!(clades_and_their_stem_lenghts[clade], shared_bl)   ## every alias keeps the measurement
            end
            push!(global_set_of_stem_bls, shared_bl)                     ## pooled background counts it once
        end
    end

    for clade::String in clade_names
        if isempty(clades_and_their_stem_lenghts[clade])
            continue
        end
        median_bl_clade::Float64=median(clades_and_their_stem_lenghts[clade])
        median_bls_for_all_clades[clade]=median_bl_clade
    end
    
    if isempty(global_set_of_stem_bls)
        error("No valid stem branch lengths found — cannot continue")
    end

    ## GLOBAL TERM = MEDIAN OF THE PER-CLADE MEDIANS — same reasoning as for
    ## terminals: the local unit is a named clade, identifiable across trees, so the
    ## constant it is compared against must also be a median of per-clade medians.
    ## Each clade votes once regardless of how many trees it was measurable in.
    ## global_set_of_stem_bls stays POOLED for the safeguard quantile.
    global_median_stem_bl= median(collect(values(median_bls_for_all_clades)))

    ## MINIMUM OBSERVATIONS GATE (Python MIN_CROSS_GENE_OBS) — same reasoning as
    ## for terminals: a clade whose stem was measurable in fewer than
    ## MIN_CROSS_GENE_OBS trees gets no trigger, so identify_clades_with_excessively
    ## _long_stems() never sees it (it iterates over the keys of this table) and the
    ## clade is retained.  Its stem lengths still contribute to global_set_of_stem_bls.
    insufficient_clades::Vector{Tuple{String,Int64}} = Vector{Tuple{String,Int64}}()
    for clade::String in clade_names
        if !haskey(median_bls_for_all_clades, clade)
            continue
        end
        n_observations::Int64 = length(clades_and_their_stem_lenghts[clade])
        if n_observations < MIN_CROSS_GENE_OBS
            push!(insufficient_clades, (clade, n_observations))
            continue
        end
        triggering_threshold_clade::Float64 = (threshold - 1) * median_bls_for_all_clades[clade] + global_median_stem_bl
        table_of_triggering_tresholds_clades[clade]=triggering_threshold_clade
    end

    if !isempty(insufficient_clades)
        sort!(insufficient_clades, by = x -> (x[2], x[1]))
        println("\nClade-stem rule not evaluable for $(length(insufficient_clades)) clade(s) (a stem was measurable in fewer than $(MIN_CROSS_GENE_OBS) trees). These clades are RETAINED by this rule:")
        for (clade, n) in insufficient_clades
            println("  $(clade) ($(n) tree$(n == 1 ? "" : "s"))")
        end
    end
    return global_set_of_stem_bls, clades_and_their_stem_lenghts, table_of_triggering_tresholds_clades  
end


function stats_for_internal_branches(trees::Dict{String,MetaGraph}, threshold::Float64)::Tuple{Dict{String,Dict{Tuple{String,String},Float64}},Vector{Float64},Dict{String,Float64}}

    tree_labels::Vector{String} = collect(keys(trees))
    global_median_internal_bl::Float64 = 0.0
    edges_of_trees::Dict{String,Vector{Tuple{String,String}}} = Dict{String,Vector{Tuple{String,String}}}()
    bls_of_trees::Dict{String,Dict{Tuple{String,String},Float64}} = Dict{String,Dict{Tuple{String,String},Float64}}()
    median_internal_bl_per_tree::Dict{String,Float64} = Dict{String,Float64}()
    ## pooled list of EVERY internal branch length — duplicates kept
    all_internal_branches_lenghts::Vector{Float64} = Vector{Float64}()
    triggering_value_per_tree::Dict{String,Float64} = Dict{String,Float64}()
    

    for tree::String in tree_labels
        edges::Vector{Tuple{String,String}} = get_internal_edges(trees[tree])
        edges_of_trees[tree] = edges
    end

    for tree::String in tree_labels
        bls_for_this_tree::Dict{Tuple{String,String},Float64} = Dict{Tuple{String,String},Float64}()
        for edge::Tuple{String,String} in edges_of_trees[tree]
            (node_a, node_b) = edge
            edge_data::EdgeData = trees[tree][node_a, node_b]
            if !ismissing(edge_data.length)
                bls_for_this_tree[edge] = edge_data.length
            end
        end
        bls_of_trees[tree] = bls_for_this_tree
    end

    for tree::String in tree_labels
        bls_current_tree::Vector{Float64} = Vector{Float64}()
        for edge::Tuple{String,String} in keys(bls_of_trees[tree])
            bl::Float64 = bls_of_trees[tree][edge]
            push!(bls_current_tree, bl)
            push!(all_internal_branches_lenghts, bl)
        end
        if !isempty(bls_current_tree)
            median_internal_bl_per_tree[tree] = median(bls_current_tree)
        end
    end

    ## GLOBAL TERM = MEDIAN OF THE PER-TREE MEDIANS, one vote per gene.
    ## The local unit here is the TREE, not a branch: individual internal branches
    ## have no identity across trees, but a gene tree does, and its median internal
    ## length captures that gene's overall rate.  Same calibration argument as for
    ## terminals and clades — the trigger is (K-1)*m_tree + M, so M must be a median
    ## of per-tree medians for a typical tree to reduce to K*M.  Pooling every branch
    ## instead would weight the constant by tree size, letting the best-sampled genes
    ## (most taxa, hence most internal branches) set it.
    ## NOTE: all_internal_branches_lenghts is still returned POOLED — the safeguard
    ## quantile is a percentile of real branch lengths and must stay that way.
    ## Named-clade stems are included throughout: an internal branch is an internal
    ## branch, and there is no natural way to carve them into groups.
    global_median_internal_bl = median(collect(values(median_internal_bl_per_tree)))

    ## per-tree trigger — same formula as terminals and stems
    ## trigger = (K-1) × local_median_of_tree + global_median
    for tree::String in tree_labels
        if !haskey(median_internal_bl_per_tree, tree)
            continue
        end
        triggering_value_per_tree[tree] = (threshold - 1) * median_internal_bl_per_tree[tree] + global_median_internal_bl
    end

    return bls_of_trees, all_internal_branches_lenghts, triggering_value_per_tree
end

## For one side of an internal edge, collect the leaves beyond `from_node` (with
## `blocked_node` blocking the other direction) together with the MEDIAN distance
## from `from_node` out to those leaves.  The tested edge itself is excluded, so
## the median measures how long the side is INTERNALLY — i.e. whether the taxa
## beyond the branch are themselves long-branched, or merely sit behind a long stem.
## Mirrors Python component_tip_names_and_distances().
function side_leaves_and_median_depth(tree::MetaGraph, from_node::String, blocked_node::String)::Tuple{Vector{String},Float64}
    leaves::Vector{String} = Vector{String}()
    leaf_distances::Vector{Float64} = Vector{Float64}()
    visited::Set{String} = Set{String}([blocked_node])
    queue::Vector{Tuple{String,Float64}} = [(from_node, 0.0)]
    while !isempty(queue)
        (current, d) = popfirst!(queue)
        if current in visited
            continue
        end
        push!(visited, current)
        if tree[current].is_a_leaf
            push!(leaves, current)
            push!(leaf_distances, d)
        end
        for nb::String in neighbor_labels(tree, current)
            if !(nb in visited)
                e::EdgeData = tree[current, nb]
                b::Float64 = ismissing(e.length) ? 0.0 : e.length
                push!(queue, (nb, d + b))
            end
        end
    end
    median_depth::Float64 = isempty(leaf_distances) ? NaN : median(leaf_distances)
    return leaves, median_depth
end

## Decide which side of a long internal edge — if either — should be removed.
## Mirrors Python detect_internal_branch_events_for_tree():
##   1. compute the median split-node-to-tip depth of each side
##   2. require the two sides to differ by more than side_ratio_threshold,
##      otherwise the branch is not informative about which side is the problem
##   3. delete the side with the LONGER median depth — the side that is itself
##      long-branched, not merely the side sitting behind a long stem
##   4. refuse unless that side is also STRICTLY SMALLER in taxon count
## Returns an empty vector when any condition blocks the removal.
function clade_for_removal_by_asymmetry(tree::MetaGraph, edge::Tuple{String,String}, side_ratio_threshold::Float64)::Vector{String}
    (u::String, v::String) = edge
    (leaves_u::Vector{String}, median_u::Float64) = side_leaves_and_median_depth(tree, u, v)
    (leaves_v::Vector{String}, median_v::Float64) = side_leaves_and_median_depth(tree, v, u)

    if isnan(median_u) || isnan(median_v) || min(median_u, median_v) <= 1e-12
        return Vector{String}()          ## Python: side_ratio undefined -> blocked
    end

    side_ratio::Float64 = max(median_u, median_v) / min(median_u, median_v)
    if side_ratio <= side_ratio_threshold
        return Vector{String}()          ## sides too similar to attribute the length
    end

    candidate::Vector{String} = Vector{String}()
    other::Vector{String} = Vector{String}()
    if median_u > median_v
        candidate = leaves_u; other = leaves_v
    elseif median_v > median_u
        candidate = leaves_v; other = leaves_u
    else
        return Vector{String}()          ## exact tie — no decision
    end

    if length(candidate) >= length(other)
        return Vector{String}()          ## Python: candidate_side_not_smaller
    end

    return candidate
end

function identify_long_branched_substrees(tree::MetaGraph, bls_of_tree::Dict{Tuple{String,String},Float64}, tree_specific_triggering_threshold::Float64, internal_branches_retention_threshold::Float64, internal_side_ratio_threshold::Float64, min_number_of_branches_in_stack::Int64)::Tuple{Vector{Vector{String}},Vector{Tuple{String,String}}}

    ## Stacks are built on the TRIGGER ALONE — the global safeguard is deliberately
    ## NOT applied here.  A branch that is long for this tree counts as evidence of
    ## a long-branch path even if it is unremarkable dataset-wide, so a modest link
    ## in the middle of a chain cannot split one stack into two and push both below
    ## min_number_of_branches_in_stack.  This mirrors Python, where a rescued edge
    ## still counts toward the downstream path-consistency tally.
    stack_of_long_branched_paths::Vector{Vector{Tuple{String,String}}}= find_long_branch_stacks(tree, bls_of_tree, tree_specific_triggering_threshold, min_number_of_branches_in_stack)

    list_of_outermost_path_branches_in_tree::Vector{Tuple{String,String}}=Vector{Tuple{String,String}}()
    list_of_clades_as_taxa_to_remove_from_tree::Vector{Vector{String}}=Vector{Vector{String}}()

    for stack::Vector{Tuple{String,String}} in stack_of_long_branched_paths
        outermost_branch_in_path::Tuple{String,String}= get_outermost_edge(tree, stack)

        ## GLOBAL SAFEGUARD — applied AT THE POINT OF ACTION, as Python's
        ## rescued_by_global_internal does: it blocks the prune, it does not erase
        ## the branch from the evidence.  The outermost edge of the stack is the
        ## branch subtending the clade that would be deleted, i.e. that clade's stem;
        ## it must exceed the global retention threshold (a quantile of ALL internal
        ## branch lengths, safeguarding_values[2]) or nothing is pruned.
        (node_x::String, node_y::String) = outermost_branch_in_path
        bl_outermost::Float64 = haskey(bls_of_tree, outermost_branch_in_path) ?
                                bls_of_tree[outermost_branch_in_path] :
                                get(bls_of_tree, (node_y, node_x), NaN)

        if isnan(bl_outermost) || bl_outermost <= internal_branches_retention_threshold
            continue   ## rescued — stack stands, but nothing is pruned from it
        end

        ## ASYMMETRY TEST — the side removed must be both the internally longer one
        ## and the smaller one.  A small clade of short-branched taxa sitting beyond
        ## a long stem is NOT removed: there the long branch is the stem, not the
        ## clade.  Returns empty when any of Python's conditions blocks the prune.
        taxa_in_current_long_clade::Vector{String}= clade_for_removal_by_asymmetry(tree, outermost_branch_in_path, internal_side_ratio_threshold)
        if isempty(taxa_in_current_long_clade)
            continue
        end

        push!(list_of_outermost_path_branches_in_tree, outermost_branch_in_path)
        push!(list_of_clades_as_taxa_to_remove_from_tree, taxa_in_current_long_clade)
    end

    return list_of_clades_as_taxa_to_remove_from_tree, list_of_outermost_path_branches_in_tree
end

### Resolve every named clade to the branch it occupies IN THIS TREE.
### For each clade returns the taxa of its largest monophyletic fragment and the
### length of the branch subtending that fragment.  Shared by the stem rule and by
### the driver, which needs the fragments to tell whether a clade was wholly removed
### by the internal-branch rule.
function resolve_clade_fragments(tree::MetaGraph, treename::String, clade_file_name::String, clades_list::Vector{String})::Tuple{Dict{String,Vector{String}},Dict{String,Float64},Dict{String,Int64}}

    ### definition of each clade
    clade_definitions::Dict{String,Vector{String}}= Dict{String,Vector{String}}()
    open(clade_file_name, "r") do fh
        for line::String in eachline(fh)
            definition::Vector{String}= split(line)
            if isempty(definition)
                continue
            end
            clade::String= popfirst!(definition)
            clade_definitions[clade]=definition
        end
    end
    
    ### reduce clade definitions to the taxa in this tree ## uses identify_clade_members() from TreeUtils 
    definitions_of_clades_in_tree::Dict{String,Vector{String}} = Dict{String,Vector{String}}()
    for clade::String in clades_list
        current_clade_full_definition::Vector{String}=clade_definitions[clade]
        results_of_clade_members=  identify_clade_members(tree, current_clade_full_definition)
        if isnothing(results_of_clade_members)   ## need checking because function can return "nothing" - if clade not in trees
            println("Clade not found in tree $(treename)")
            continue
        end
        (definition_of_clade_in_tree, missing_taxa)= results_of_clade_members 
        if !isempty(definition_of_clade_in_tree) && !isempty(missing_taxa)
            println("Some taxa are mising from $(clade) in tree $(treename).  These are: $(missing_taxa)")
        end
        definitions_of_clades_in_tree[clade]=definition_of_clade_in_tree 
    end
    
    ### find branch lenght of stem of clades, and the fragment each clade resolves to
    clades_branch_lengths::Dict{String,Float64}= Dict{String,Float64}()
    clade_fragments::Dict{String,Vector{String}}= Dict{String,Vector{String}}()
    for clade::String in clades_list
        if !haskey(definitions_of_clades_in_tree, clade) ### not all clades in clade list will be in current tree.  So we need checking and skipping
            continue
        end
        largest_subset::Vector{String} = definitions_of_clades_in_tree[clade]
        ## clade name passed for the message; warning suppressed because the println
        ## just below already reports this by name, with the tree
        lca_of_current_clade= find_lca(tree, definitions_of_clades_in_tree[clade], clade; warn_on_failure=false)
        if isnothing(lca_of_current_clade)
            println("$(clade) not monophyletic in tree $(treename).  We will define clade as the largest monophyletic subset of taxa of $(clade) in $(treename)")
            
            (lca_of_current_clade, largest_subset) = find_lca_largest_subset(tree, definitions_of_clades_in_tree[clade], clade)
        end
        if isnothing(lca_of_current_clade)
            println("Warning: no subset of $(clade) including 2 or more taxa found in $(treename). We will skip.")
            continue
        end
        edge_subtending_lca= identify_stem_edge(tree, lca_of_current_clade, largest_subset)
        if isnothing(edge_subtending_lca)
            println("Warning: no edge found to subtend  $(clade). Something is wrong with tree. We will skip. But this is an anomalous result.  Do not trust results of analysis")
            continue
        end
        (node_a, node_b)= edge_subtending_lca
        edge_info::EdgeData= tree[node_a, node_b]
        bl::Float64= edge_info.length
        clades_branch_lengths[clade]= bl
        clade_fragments[clade]= sort(largest_subset)
    end

    ## size of each clade's FULL definition in the file — used to decide which of
    ## several names sharing a branch describes it most tightly
    clade_definition_sizes::Dict{String,Int64} = Dict{String,Int64}()
    for (name, definition) in clade_definitions
        clade_definition_sizes[name] = length(definition)
    end

    return clade_fragments, clades_branch_lengths, clade_definition_sizes
end


function identify_clades_with_excessively_long_stems(tree::MetaGraph, treename::String, stem_specific_triggering_tresholds::Dict{String,Float64}, stem_retention_threshold::Float64, clade_file_name::String)::Dict{String,Vector{String}}
    ### return value Dictionary with long branched clades definitions
    long_clades::Dict{String,Vector{String}}= Dict{String,Vector{String}}()

    ### only clades that got a trigger (i.e. passed MIN_CROSS_GENE_OBS) are testable
    clades_list::Vector{String}= collect(keys(stem_specific_triggering_tresholds))

    (clade_fragments, clades_branch_lengths, clade_definition_sizes) = resolve_clade_fragments(tree, treename, clade_file_name, clades_list)

    ### ALIASING: group the clade names that resolved to the SAME branch in this tree
    ### so one physical branch is tested once rather than once per name.
    names_by_fragment::Dict{String,Vector{String}}= Dict{String,Vector{String}}()
    for clade::String in sort(collect(keys(clade_fragments)))
        fragment_key::String = join(clade_fragments[clade], "\u0001")
        push!(get!(names_by_fragment, fragment_key, Vector{String}()), clade)
    end

    ### test one representative per branch.
    ### Which name supplies the threshold?  Each name carries its own cross-gene
    ### median, so a shared branch has several candidate thresholds.  We use the name
    ### whose FULL DEFINITION IS SMALLEST — the one that describes this branch most
    ### tightly.  Example: CNID is split across the tree and its largest monophyletic
    ### fragment is the CNID_A taxa, so CNID and CNID_A share a branch.  CNID_A's
    ### median comes from the same taxa in every tree and is therefore comparable to
    ### the branch being measured, whereas CNID's median mixes stems of the full clade
    ### with stems of whatever fragment survived in each tree.  Trigger on the clean
    ### estimate even though, for reporting, the branch is named for the wider clade.
    candidate_groups::Vector{Tuple{Vector{String},Vector{String},Float64}} = Vector{Tuple{Vector{String},Vector{String},Float64}}()
    for fragment_key::String in sort(collect(keys(names_by_fragment)))
        names_sharing_branch::Vector{String} = sort(names_by_fragment[fragment_key])
        shared_bl::Float64 = clades_branch_lengths[names_sharing_branch[1]]

        tightest_name::String = sort(names_sharing_branch, by = n -> (get(clade_definition_sizes, n, typemax(Int64)), n))[1]
        trigger_for_branch::Float64 = stem_specific_triggering_tresholds[tightest_name]

        if shared_bl > trigger_for_branch && shared_bl > stem_retention_threshold
            fragment_taxa::Vector{String} = clade_fragments[names_sharing_branch[1]]
            push!(candidate_groups, (names_sharing_branch, fragment_taxa, shared_bl))
        end
    end

    ### NESTED CLADE SUPPRESSION: if a branch is dropped, any dropped branch whose
    ### taxa are a strict subset of it is redundant — those taxa are already going.
    ### Largest first, so the outermost branch wins and the nested one is dropped
    ### from the report rather than listing the same taxa twice.
    sort!(candidate_groups, by = g -> -length(g[2]))
    accepted::Vector{Tuple{Vector{String},Vector{String},Float64}} = Vector{Tuple{Vector{String},Vector{String},Float64}}()
    for group in candidate_groups
        taxa_set::Set{String} = Set{String}(group[2])
        nested_in_accepted::Bool = any(issubset(taxa_set, Set{String}(a[2])) for a in accepted)
        if nested_in_accepted
            println("Clade(s) $(join(group[1], "|")) in tree $(treename): stem is long but its taxa are already covered by a larger dropped clade — not reported separately.")
            continue
        end
        push!(accepted, group)
    end

    for group in accepted
        ## key names every clade that resolved to this branch, so the report does not
        ## silently pick one name over its aliases
        long_clades[join(group[1], "|")] = group[2]
    end
    return long_clades
end

function identify_long_terminals(tree::MetaGraph, treename::String, taxon_specific_triggering_tresholds::Dict{String,Float64}, terminal_taxa_retention_threshold::Float64)::Vector{String}
    println("Processing tree:$(treename)")
    terminal_taxa_in_tree::Vector{String}= collect(get_leaves(tree))
    list_long_branched_taxa_in_tree::Vector{String}=Vector{String}()
    
    for taxon::String in terminal_taxa_in_tree
        ## no trigger issued means the taxon was seen in fewer than MIN_CROSS_GENE_OBS
        ## trees: not evaluable, therefore retained
        if !haskey(taxon_specific_triggering_tresholds, taxon)
            continue
        end
        neighbour_of_leaf::String= first(neighbor_labels(tree, taxon)) ## because these are leaf there is only one neighbour in neighbour labels so this cannot go wrong.                        
        edge::EdgeData= tree[taxon, neighbour_of_leaf] ## this uses structured defined in TreeUtil (i.e. edge is of type EdgeData                                                                
        bl::Float64 = edge.length  ## this gest the value of edge "."    
        if bl > taxon_specific_triggering_tresholds[taxon] && bl > terminal_taxa_retention_threshold
            push!(list_long_branched_taxa_in_tree, taxon)
        end
    end
    return list_long_branched_taxa_in_tree
end

## -----------------------------------------------------------------------
## Quartet rule — local long branch detection
## -----------------------------------------------------------------------

## Derive the masked clades AS THEY APPEAR IN THIS TREE.
## Mirrors Python detect_mask_groups_from_graph(): for each named clade keep the
## LARGEST split-side whose leaves are all members of that clade — i.e. the
## largest monophyletic fragment of the clade in this tree — provided it holds at
## least `min_size` taxa.  A clade with fewer than `min_size` taxa present, or
## one that is too scattered to form a fragment of that size, is NOT masked at
## all in this tree.
function detect_mask_groups_for_tree(tree::MetaGraph, masked_clades::Dict{String,Vector{String}},
                                     min_size::Int64=2)::Dict{String,Set{String}}
    mask_groups::Dict{String,Set{String}} = Dict{String,Set{String}}()
    if isempty(masked_clades)
        return mask_groups
    end
    tips_in_tree::Set{String} = get_leaves(tree)
    sides::Vector{Set{String}} = all_split_leaf_sets(tree)

    for clade_name::String in sort(collect(keys(masked_clades)))
        ## strip empty fields produced by runs of spaces in the clade file
        members::Set{String} = Set{String}(filter(t -> !isempty(strip(t)), masked_clades[clade_name]))
        present_taxa::Set{String} = intersect(members, tips_in_tree)
        if length(present_taxa) < min_size
            continue
        end
        best_fragment::Set{String} = Set{String}()
        for side::Set{String} in sides
            if length(side) >= min_size && issubset(side, present_taxa) && length(side) > length(best_fragment)
                best_fragment = side
            end
        end
        if length(best_fragment) >= min_size
            mask_groups[clade_name] = best_fragment
        end
    end
    return mask_groups
end

## Map each taxon to the ONE masked clade it belongs to IN THIS TREE.
## Mirrors Python build_taxon_to_group() applied to the per-tree mask groups.
## Python raises on overlap; here we warn and keep the alphabetically first
## clade so a batch run is not killed by one bad line.
function build_quartet_mask_map(tree::MetaGraph, masked_clades::Dict{String,Vector{String}})::Dict{String,String}
    taxon_to_group::Dict{String,String} = Dict{String,String}()
    mask_groups::Dict{String,Set{String}} = detect_mask_groups_for_tree(tree, masked_clades)
    for clade_name::String in sort(collect(keys(mask_groups)))
        for taxon::String in sort(collect(mask_groups[clade_name]))
            if haskey(taxon_to_group, taxon) && taxon_to_group[taxon] != clade_name
                @warn "Taxon $(taxon) falls in more than one retained quartet-mask fragment ($(taxon_to_group[taxon]) and $(clade_name)) — keeping $(taxon_to_group[taxon])"
                continue
            end
            taxon_to_group[taxon] = clade_name
        end
    end
    return taxon_to_group
end

## Python comparator_allowed():
##   focal NOT in a masked clade  -> comparator must also be unmasked
##                                   (a normal taxon is never compared against a
##                                    known long-branch taxon)
##   focal IN masked clade G      -> comparator must be unmasked OR in the same
##                                   clade G (a ctenophore may be compared
##                                   against other ctenophores, but not against
##                                   members of a different masked clade)
function comparator_allowed(taxon_to_group::Dict{String,String}, focal::String, comparator::String)::Bool
    focal_group::Union{String,Nothing} = get(taxon_to_group, focal, nothing)
    comp_group::Union{String,Nothing}  = get(taxon_to_group, comparator, nothing)
    if isnothing(focal_group)
        return isnothing(comp_group)
    end
    return isnothing(comp_group) || comp_group == focal_group
end

## Walk outward from the focal tip collecting internal nodes, LEVEL BY LEVEL.
## Mirrors Python candidate_anchors_for_tip(): nodes are visited in breadth-first
## levels and sorted within each level by vertex code, so the anchor order is
## deterministic and runs from nearest to furthest.
function candidate_anchors_for_tip(tree::MetaGraph, focal::String)::Vector{String}
    anchors::Vector{String} = Vector{String}()
    visited::Set{String} = Set{String}([focal])
    current_level::Vector{String} = collect(neighbor_labels(tree, focal))
    while !isempty(current_level)
        sort!(current_level, by = n -> code_for(tree, n))
        next_level::Vector{String} = Vector{String}()
        for node::String in current_level
            if node in visited
                continue
            end
            push!(visited, node)
            if !tree[node].is_a_leaf  ## only internal nodes are valid anchors
                push!(anchors, node)
            end
            for nb::String in neighbor_labels(tree, node)
                if !(nb in visited)
                    push!(next_level, nb)
                end
            end
        end
        current_level = next_level
    end
    return anchors
end

## Single pass around an anchor: for every other node in the tree record
##   (a) which neighbour-branch of the anchor it lies beyond, and
##   (b) its patristic distance from the anchor.
## Replaces the repeated whole-tree BFS the previous version ran per comparator.
function anchor_branch_partition(tree::MetaGraph, anchor::String)::Tuple{Dict{String,String},Dict{String,Float64}}
    branch_of::Dict{String,String} = Dict{String,String}()
    dist_from_anchor::Dict{String,Float64} = Dict{String,Float64}()
    for nb::String in neighbor_labels(tree, anchor)
        first_edge::EdgeData = tree[anchor, nb]
        first_bl::Float64 = ismissing(first_edge.length) ? 0.0 : first_edge.length
        visited::Set{String} = Set{String}([anchor])
        queue::Vector{Tuple{String,Float64}} = [(nb, first_bl)]
        while !isempty(queue)
            (current, d) = popfirst!(queue)
            if current in visited
                continue
            end
            push!(visited, current)
            branch_of[current] = nb
            dist_from_anchor[current] = d
            for nb2::String in neighbor_labels(tree, current)
                if !(nb2 in visited)
                    e::EdgeData = tree[current, nb2]
                    b::Float64 = ismissing(e.length) ? 0.0 : e.length
                    push!(queue, (nb2, d + b))
                end
            end
        end
    end
    return branch_of, dist_from_anchor
end

## At a given anchor, drop the branch leading back to the focal tip and, for each
## remaining branch, collect the tips that are allowed comparators for this focal
## taxon.  Tips within a group are sorted by (distance from anchor, name) so that
## group[1] is the CLOSEST valid comparator — matching Python
## build_allowed_groups_at_anchor().  Returns one group per branch, branches
## ordered deterministically by vertex code.
function build_comparator_groups_at_anchor(tree::MetaGraph, focal::String, anchor::String,
                                           taxon_to_group::Dict{String,String},
                                           branch_of::Dict{String,String},
                                           dist_from_anchor::Dict{String,Float64};
                                           min_comparator_branch::Float64 = MIN_COMPARATOR_BRANCH)::Vector{Vector{String}}

    focal_branch::Union{String,Nothing} = get(branch_of, focal, nothing)
    if isnothing(focal_branch)
        return Vector{Vector{String}}()
    end

    tips_by_branch::Dict{String,Vector{String}} = Dict{String,Vector{String}}()
    for (node::String, branch::String) in branch_of
        if branch == focal_branch          ## skip the focal direction
            continue
        end
        if !tree[node].is_a_leaf
            continue
        end
        if !comparator_allowed(taxon_to_group, focal, node)
            continue
        end
        ## a tip whose own terminal branch is at the inference floor is unestimated,
        ## not short, and carries no information about the local scale
        neighbour_of_tip::String = first(neighbor_labels(tree, node))
        tip_edge::EdgeData = tree[node, neighbour_of_tip]
        if ismissing(tip_edge.length) || tip_edge.length <= min_comparator_branch
            continue
        end
        push!(get!(tips_by_branch, branch, Vector{String}()), node)
    end

    groups::Vector{Vector{String}} = Vector{Vector{String}}()
    for branch::String in sort(collect(keys(tips_by_branch)), by = n -> code_for(tree, n))
        tips::Vector{String} = tips_by_branch[branch]
        sort!(tips, by = t -> (dist_from_anchor[t], t))
        push!(groups, tips)
    end
    return groups
end

## Pick exactly THREE comparators from the groups at an anchor, as Python
## choose_quartet_from_groups() does:
##   - require at least 2 groups and at least 3 available tips in total
##   - take the closest tip from each group first (one per branch)
##   - if that gives more than 3, keep the 3 closest
##   - if it gives fewer than 3, top up from the pooled remainder, closest first
## Returns nothing if 3 comparators cannot be assembled — the caller then moves
## on to the next anchor.  This is the key difference from the previous version,
## which accepted 2 comparators and never topped up.
function choose_quartet_from_groups(groups::Vector{Vector{String}},
                                    dist_from_anchor::Dict{String,Float64})::Union{Vector{String},Nothing}
    if length(groups) < 2
        return nothing
    end
    if sum(length(g) for g in groups) < 3
        return nothing
    end

    chosen::Vector{String} = [g[1] for g in groups]
    if length(chosen) > 3
        sort!(chosen, by = t -> (dist_from_anchor[t], t))
        chosen = chosen[1:3]
    end

    remaining::Vector{String} = Vector{String}()
    for g::Vector{String} in groups
        if length(g) > 1
            append!(remaining, g[2:end])
        end
    end
    sort!(remaining, by = t -> (dist_from_anchor[t], t))

    i::Int64 = 1
    while length(chosen) < 3 && i <= length(remaining)
        push!(chosen, remaining[i])
        i += 1
    end

    return length(chosen) == 3 ? chosen : nothing
end

## Main quartet test function — for a single focal taxon in a single tree.
## Returns (is_long_branch::Bool, local_ratio::Float64)
## Returns (false, NaN) if no valid quartet could be built at any anchor.
##
## Mirrors Python find_mask_aware_quartet_for_tip() + the scoring block of
## score_taxa(): anchors are tried nearest-first and the FIRST anchor that yields
## three valid comparators is used — the search does not continue past it.
function quartet_test_for_taxon(tree::MetaGraph, focal::String,
                                taxon_to_group::Dict{String,String},
                                ratio_threshold::Float64;
                                min_comparator_branch::Float64 = MIN_COMPARATOR_BRANCH)::Tuple{Bool,Float64}

    for anchor::String in candidate_anchors_for_tip(tree, focal)
        (branch_of, dist_from_anchor) = anchor_branch_partition(tree, anchor)

        if !haskey(dist_from_anchor, focal)
            continue
        end

        groups::Vector{Vector{String}} = build_comparator_groups_at_anchor(tree, focal, anchor,
                                                                          taxon_to_group,
                                                                          branch_of, dist_from_anchor;
                                                                          min_comparator_branch = min_comparator_branch)
        comparators::Union{Vector{String},Nothing} = choose_quartet_from_groups(groups, dist_from_anchor)
        if isnothing(comparators)
            continue
        end

        d_focal::Float64 = dist_from_anchor[focal]
        comp_distances::Vector{Float64} = [dist_from_anchor[c] for c in comparators]
        local_median::Float64 = median(comp_distances)

        ## Even with degenerate comparators excluded, refuse to divide by a local
        ## reference that carries no scale.  Not evaluable, therefore retained.
        if local_median <= min_comparator_branch
            return (false, NaN)
        end

        local_ratio::Float64 = d_focal / local_median
        ## strict > : a taxon sitting exactly ON the threshold is RETAINED
        return (local_ratio > ratio_threshold, local_ratio)
    end

    ## no valid quartet found for this taxon
    return (false, NaN)
end

## Convenience method keeping the original call signature — derives this tree's
## mask fragments on the fly.  Prefer the method above inside loops: build the
## map ONCE PER TREE with build_quartet_mask_map() and pass it in, because the
## fragment search traverses every split in the tree.
function quartet_test_for_taxon(tree::MetaGraph, focal::String,
                                masked_clades::Dict{String,Vector{String}},
                                ratio_threshold::Float64;
                                min_comparator_branch::Float64 = MIN_COMPARATOR_BRANCH)::Tuple{Bool,Float64}
    return quartet_test_for_taxon(tree, focal, build_quartet_mask_map(tree, masked_clades), ratio_threshold;
                                  min_comparator_branch = min_comparator_branch)
end

## -----------------------------------------------------------------------
## DIAGNOSTIC ONLY — not part of any rule.
## Re-runs the anchor search for one taxon and reports the components of the
## ratio, so an extreme value can be attributed to its cause: a collapsed
## denominator (comparators sitting on top of the anchor) or an inflated
## numerator (the anchor having walked far from the focal taxon).
## Delete this together with the debug block in the statistics script.
## -----------------------------------------------------------------------
function quartet_diagnostics(tree::MetaGraph, focal::String,
                             taxon_to_group::Dict{String,String},
                             ratio_threshold::Float64;
                             min_comparator_branch::Float64 = MIN_COMPARATOR_BRANCH)::String
    for anchor::String in candidate_anchors_for_tip(tree, focal)
        (branch_of, dist_from_anchor) = anchor_branch_partition(tree, anchor)
        if !haskey(dist_from_anchor, focal)
            continue
        end
        groups = build_comparator_groups_at_anchor(tree, focal, anchor, taxon_to_group,
                                                   branch_of, dist_from_anchor;
                                                   min_comparator_branch = min_comparator_branch)
        comparators = choose_quartet_from_groups(groups, dist_from_anchor)
        if isnothing(comparators)
            continue
        end

        d_focal::Float64 = dist_from_anchor[focal]
        comp_distances::Vector{Float64} = [dist_from_anchor[c] for c in comparators]
        local_median::Float64 = median(comp_distances)
        edges_away = unweighted_path_distance(tree, focal, anchor)

        focal_nb::String = first(neighbor_labels(tree, focal))
        focal_edge::EdgeData = tree[focal, focal_nb]
        focal_bl = ismissing(focal_edge.length) ? NaN : focal_edge.length

        parts::Vector{String} = Vector{String}()
        for c::String in comparators
            nb::String = first(neighbor_labels(tree, c))
            e::EdgeData = tree[c, nb]
            own_bl = ismissing(e.length) ? NaN : e.length
            push!(parts, string(c,
                                " d_from_anchor=", round(dist_from_anchor[c], sigdigits=4),
                                " own_bl=", round(own_bl, sigdigits=4)))
        end

        return string("anchor=", anchor,
                      "  edges_focal_to_anchor=", isnothing(edges_away) ? "NA" : edges_away,
                      "  focal_own_bl=", round(focal_bl, sigdigits=4),
                      "  d_focal=", round(d_focal, sigdigits=6),
                      "  local_median=", round(local_median, sigdigits=6),
                      "  n_groups=", length(groups),
                      "\n        comparators: ", join(parts, "  |  "))
    end
    return "no valid quartet at any anchor"
end

## Top-level function — identify locally long terminal taxa by the quartet rule.
##
## masked_clades:  clade name -> member taxa (from excluded_from_quartet_tests file)
## ratio_threshold: K value (default 3.0) — the LOCAL test
## terminal_bl_retention_threshold: the SAFEGUARD, a quantile of ALL terminal
##   branch lengths in the dataset (Python --global_rescue / safeguarding_values[1])
##
## Two orthogonal questions, as in Python:
##   K       — is this taxon extreme *relative to its local neighbourhood*?
##   safeguard — is its branch actually long enough *in absolute terms* to matter?
## Python drops on the ratio and then rescues the taxon when its own terminal
## branch falls below the global cutoff:
##     if branch_length < global_rescue_cutoff: decision = "KEEP"
## which is the same as requiring the branch to reach the cutoff before dropping.
##
## The safeguard is deliberately on BRANCH LENGTH, not on the ratio: in a very
## short-branched neighbourhood a trivial branch can produce an enormous ratio
## (comparators at 0.001, focal at 0.01 -> ratio 10) and such a branch cannot
## drive long-branch attraction, which depends on absolute amounts of change.
## Testing the ratio against a quantile of ratios would not catch that case, and
## would merely restate K as max(K, quantile).
##
## This also subsumes the old hard-coded `edge.length <= 1e-6` skip: that was the
## same test with an arbitrary cutoff instead of one taken from the data.
function identify_local_long_branches(tree::MetaGraph, treename::String, masked_clades::Dict{String,Vector{String}}, ratio_threshold::Float64, terminal_bl_retention_threshold::Float64; min_comparator_branch::Float64 = MIN_COMPARATOR_BRANCH)::Vector{String}
    long_taxa::Vector{String} = Vector{String}()
    ## mask fragments are tree-specific — derive them once for this tree
    taxon_to_group::Dict{String,String} = build_quartet_mask_map(tree, masked_clades)

    ## coverage bookkeeping: a taxon for which no valid quartet can be assembled is
    ## NOT tested, which is not the same as tested-and-passed.  It is retained either
    ## way, but the user should know how much of the tree the rule actually reached.
    not_evaluable::Vector{String} = Vector{String}()
    n_tested::Int64 = 0

    for taxon::String in collect(get_leaves(tree))
        neighbour::String = first(neighbor_labels(tree, taxon))
        edge::EdgeData = tree[taxon, neighbour]
        if ismissing(edge.length)
            push!(not_evaluable, taxon)
            continue
        end
        terminal_bl::Float64 = edge.length

        (is_long, ratio) = quartet_test_for_taxon(tree, taxon, taxon_to_group, ratio_threshold; min_comparator_branch = min_comparator_branch)

        if isnan(ratio)
            ## no anchor yielded three valid comparators after masking, or the local
            ## reference collapsed to zero
            push!(not_evaluable, taxon)
            continue
        end
        n_tested += 1

        ## locally extreme AND absolutely long enough to be worth removing
        ## strict > on the safeguard too — in doubt, retain
        if is_long && terminal_bl > terminal_bl_retention_threshold
            push!(long_taxa, taxon)
        end
    end

    n_leaves::Int64 = n_tested + length(not_evaluable)
    if !isempty(not_evaluable)
        sort!(not_evaluable)
        shown::Vector{String} = length(not_evaluable) > 20 ? not_evaluable[1:20] : not_evaluable
        extra::String = length(not_evaluable) > 20 ? " ... and $(length(not_evaluable) - 20) more" : ""
        println("Quartet rule in $(treename): tested $(n_tested)/$(n_leaves) taxa. Not evaluable (no valid quartet after masking) and therefore RETAINED: $(join(shown, ", "))$(extra)")
    else
        println("Quartet rule in $(treename): tested $(n_tested)/$(n_leaves) taxa.")
    end

    return long_taxa
end
