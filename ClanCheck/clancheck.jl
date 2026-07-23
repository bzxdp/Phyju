using MetaGraphsNext
using TreesIO
using TreeUtils
using ArgParse


function main()
    args= parse_arguments()
    treefile_extension= args["extension"]
    outfile= args["output"]
    clades_file= args["clades"]

    
    list_of_treefiles= filter(f->endswith(f, treefile_extension), readdir("."))
    if isempty(list_of_treefiles)
        println("No files found with extension $treefile_extension")
        exit(1)
    end
    trees=Dict{String,MetaGraph}()
    for file in list_of_treefiles
        trees[file]= readtree(file, nothing)  ## need to change readtree to make bl optional
    end
    
    clades_to_check::Dict{String,Vector{String}}= Dict{String,Vector{String}}()
    open(clades_file, "r") do fh
        for clade in eachline(fh)
            clade_definition= split(clade, " ")
            clade_name = popfirst!(clade_definition)
            clades_to_check[clade_name]=clade_definition
        end
    end
    tree_labels= collect(keys(trees))
    clade_labels= collect(keys(clades_to_check))
    open(outfile, "w") do fh_out
        for tree in tree_labels
            leaves::Set{String}= Set(get_leaves(trees[tree]))
            Test_clades_in_tree::Vector{Bool}= Vector{Bool}()
            for clade in clade_labels
                clade_as_a_set::Set{String}= Set(clades_to_check[clade])
                clade_taxa_in_tree= collect(intersect(clade_as_a_set,leaves))
                if isempty(clade_taxa_in_tree) || length(clade_taxa_in_tree) == 1  ## This is an important step.  If a clade is not in a tree it cannot be tested and cannot be rejected so its an automatic pass.  If it includes only one taxon, it does not have a LCA in the standard sense to find_lca will fail. But the clade is monophyletic, so we just consider these two conditions as automatic passes.
                    push!(Test_clades_in_tree, true)
                    @warn "$(tree) does not include clade $(clade). This is marked as a pass in this implementation/interpretation of  clancheck" 
                    continue
                end
                lca_current_clade= find_lca(trees[tree], clade_taxa_in_tree)
                if !isnothing(lca_current_clade)
                    push!(Test_clades_in_tree,  true)
                else
                    push!(Test_clades_in_tree, false)
                end
            end
            if all(Test_clades_in_tree)
                println(fh_out, "$(tree)")
            end
        end
    end
end
    
function parse_arguments()
    s =ArgParseSettings(description="reformat sequences")
    @add_arg_table s begin
        "--extension", "-e"
        help = "extension of the files is usually: tree, tre, treefile or similar"
        required = true
        "--output", "-o"
        help = "output file"
        required = true
        "--clades", "-c"
        help = "file with clades to test: one line for clade. Name of clade first then all taxa in clade. all separated by spaces (no spaces allowed in names of clades of species: e.g. Clade_A NOT Clade A and Taxon_1 NOT Taxon 1)." 
        required = true
    end
    return parse_args(s)
end

main()

