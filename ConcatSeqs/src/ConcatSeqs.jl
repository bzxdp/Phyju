module ConcatSeqs

### simple script to concatenate genes - transformed into a module for binary compilation

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

## to make it into binary machine code we need main to return a Int32 (a zero for success).
## Note that is has to be 32 not 64 or else will crash
## has to adhere to C ABI (Application Binary Interface) standard, which is
## Universal across all Operating Systems.
## i.e. Op systs expect a bunary code that work correctly to return 0 as an INT32.
## Cint is alyas for Int32
## it is convetion to use it in jula to wrap main in a function called julia_main
## This is entry point and Cint signify that we are returning 0 as 32 bit binary to adhere to the
## C 32-binary standard as this is expected by Op Syst.
## it will equally work using Int32
## But this explicitly tells the reader why we are using Int32 on a system that may have 64 or more bits.
## its because of the C ABI (Application Binary Interface) convention

function julia_main()::Cint ## this is just an alias for Int32 
    main()
    return 0
end

end #module



        

