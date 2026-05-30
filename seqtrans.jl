using ArgParse
using BioSequences
using PhyloDataIO
using FASTX

function main()
    args= parse_arguments()
    infile= args["input"]
    outfile= args["output"]
    format= args["outformat"]

    sequences=  read_PhyloData(infile)
    write_PhyloData(sequences, outfile, format)
end

function parse_arguments()
    s =ArgParseSettings(description="reformat sequences")
    @add_arg_table s begin
        "--input", "-i"
            help = "input file"
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
