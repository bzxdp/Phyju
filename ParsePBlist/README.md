# ParsePBlist

Parses PhyloBayes `.bplist` output files to extract and report support values for user-defined clades across multiple MCMC chains.
Precompiling is only useful if it is to be used for simulations where you might have to summarise results of hundreds of experiments. 
For single dataset analyses use the standard julia script provided at the root directory 

Two usage options are available:
- **Script**: run `ParsePBlist.jl` directly with Julia (requires Julia and dependencies installed)
- **Binary**: compile `ParsePBlist` for faster startup and no Julia installation required for end users

## Requirements

### For script usage
- [Julia](https://julialang.org/downloads/) >= 1.11
- Dependencies installed in your Julia environment (see below)

### For binary compilation
- [Julia](https://julialang.org/downloads/) >= 1.11
- [PackageCompiler.jl](https://github.com/JuliaLang/PackageCompiler.jl)

## Dependencies

| Package | Source |
|---------|--------|
| [ArgParse.jl](https://github.com/carlobaldassi/ArgParse.jl) | Julia General Registry |
| Statistics | Julia Standard Library (bundled with Julia) |

## Option 1: Run as a script

### 1. Install dependencies

Start Julia and run:

```julia
using Pkg
Pkg.add(["ArgParse", "Statistics"])
```

### 2. Run the script

```bash
julia /path/to/Phyju/ParsePBlist/ParsePBlist.jl -i <bplist_file> -o <outfile> -n <nchains> -c <clades_file>
```

## Option 2: Compile as a binary

### 1. Clone this repository

```bash
git clone https://github.com/bzxdp/Phyju.git
cd Phyju/ParsePBlist
```

### 2. Start Julia and activate the project

```julia
using Pkg
Pkg.activate(".")
```

### 3. Instantiate dependencies

```julia
Pkg.instantiate()
```

### 4. Add PackageCompiler

```julia
Pkg.add("PackageCompiler")
```

### 5. Build the binary

```julia
using PackageCompiler
create_app(".", "ParsePBlist_binary", force=true)
```

This will take approximately 5–10 minutes. The compiled binary will be at:

```
ParsePBlist_binary/bin/ParsePBlist
```

### 6. Add the binary to your PATH

Add the following to your `~/.bashrc` (or `~/.zshrc`):

```bash
export PATH="/path/to/Phyju/ParsePBlist/ParsePBlist_binary/bin:$PATH"
```

Then reload:

```bash
source ~/.bashrc
```

> **Note:** The entire `ParsePBlist_binary/` directory must remain intact — do not copy just the `ParsePBlist` executable alone. The binary is platform-specific and must be compiled on the target architecture (Linux x86_64, macOS ARM, etc.).

## Usage

Both the script and the binary accept the same arguments:

```bash
ParsePBlist -i <bplist_file> -o <outfile> -n <nchains> -c <clades_file>
# or
julia ParsePBlist.jl -i <bplist_file> -o <outfile> -n <nchains> -c <clades_file>
```

### Arguments

| Flag | Description | Required |
|------|-------------|----------|
| `-i` / `--infile` | PhyloBayes `.bplist` output file | Yes |
| `-o` / `--outfile` | Output file name | Yes |
| `-n` / `--number_of_chains` | Number of chains summarised in the `.bplist` file | Yes |
| `-c` / `--clades_to_test` | Text file defining clades to search (see format below) | Yes |

### Clades file format

One clade per line. First element is the clade name, followed by all species in the clade, all separated by spaces:

```
CLADEa TaxonA TaxonB TaxonC TaxonD
CLADEb TaxonE TaxonF TaxonG
```

### Example

```bash
ParsePBlist -i run1run2.bplist -o clade_supports.txt -n 2 -c clades_to_test.txt
```

### Output format

```
Average Support for CLADEa in 2 chains = 95.0 -- individual chains support = [94, 96]
Clade CLADEb not found: Support is 0
```

## Notes

- Species names in the clades file must match exactly the taxon names in the `.bplist` file
- Both sides of each bipartition are searched, so it does not matter whether your clade of interest is coded as `*` or `.` in the PhyloBayes output
- `ParsePBlist_binary/` should be added to `.gitignore` — it is large and platform-specific
