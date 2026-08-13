#!/usr/bin/env julia

### simple script to concatenate genes
### will be replaced when some functions
using PhyloDataIO
using Glob
using ArgParse

function main()
    args = parse_arguments()
    extension = args["extension"]
    outfile = args["outfile"]
    outfile_format = args["format"]
    genefiles= glob("*.$(extension)")


    ## store single gene alignments
    all_genes::Dict{String,Dict{String,String}}= read_many_gene_alignments(genefiles) ## call bulk read by extension from TreesIO
    ##

    ## make concatenation and list of partitions
    (concatenated_alignment, partitions)= concatenate(all_genes)                
    
    ## write output outputs
    ## write concatenated alignment
    write_PhyloData(concatenated_alignment, outfile, outfile_format)

    ## write partition boundaries
    (partition_file_basename, _)= splitext(outfile) ## splitext splits "ext" from file name . 
    partition_file= partition_file_basename * "_partitions.txt"
    open(partition_file, "w") do fh
        for (gene, (start, stop)) in partitions ### here I dereference directly the touple.  you could do two steps (gene, pointer) then (start, stop) = pointer 
            println(fh, "$(gene) = $(start) - $(stop);")
        end
    end
end

function parse_arguments()
    s =ArgParseSettings(description="concatenate sequences")
    @add_arg_table s begin
        "--extension", "-e"
        help = "extension of files to concatenate e.g. fas"
        required = true
        "--outfile", "-o"
            help = "outfile name"
            required =true
        "--format", "-f"
        help = "format of outfile: can be f=fasta, p=phylip, n=nexus"
        required =true
    end
    return parse_args(s)
end
main()



        

