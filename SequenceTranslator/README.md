# SequenceTranslator

Translates sequence alignments between different formats: FASTA, Nexus, Phylip, and plain text table.

Two usage options are available:
- **Script**: run `seqtrans.jl` directly with Julia — sufficient for occasional single-file conversions
- **Binary**: compile `SequenceTranslator` for use in loops over large numbers of files

> **When do you need the binary?** When translating hundreds or thousands of files (e.g. all gene alignments across a set of genomes), running the plain script in a `for` loop becomes very slow because Julia recompiles the code from scratch at every invocation. The precompiled binary eliminates this overhead entirely, making large-scale batch conversions practical.

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
| [PhyloDataIO](https://github.com/bzxdp/PhyloDataIO) | GitHub (unregistered) |
| [ArgParse.jl](https://github.com/carlobaldassi/ArgParse.jl) | Julia General Registry |
| [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl) | Julia General Registry |
| [FASTX.jl](https://github.com/BioJulia/FASTX.jl) | Julia General Registry |

## Option 1: Run as a script

### 1. Install dependencies

Start Julia and run:

```julia
using Pkg
Pkg.develop(url="https://github.com/bzxdp/PhyloDataIO")
Pkg.add(["ArgParse", "BioSequences", "FASTX"])
```

### 2. Run the script

```bash
julia /path/to/Phyju/SequenceTranslator/seqtrans.jl -i <infile> -o <outfile> -f <format>
```

## Option 2: Compile as a binary

### 1. Clone this repository

```bash
git clone https://github.com/bzxdp/Phyju.git
cd Phyju/SequenceTranslator
```

### 2. Start Julia and activate the project

```julia
using Pkg
Pkg.activate(".")
```

### 3. Add PhyloDataIO (unregistered dependency)

```julia
Pkg.develop(url="https://github.com/bzxdp/PhyloDataIO")
```

### 4. Instantiate remaining dependencies

```julia
Pkg.instantiate()
```

### 5. Add PackageCompiler

```julia
Pkg.add("PackageCompiler")
```

### 6. Build the binary

```julia
using PackageCompiler
create_app(".", "SequenceTranslator_binary")
```

This will take approximately 5–10 minutes. The compiled binary will be at:

```
SequenceTranslator_binary/bin/SequenceTranslator
```

### 7. Add the binary to your PATH

Add the following to your `~/.bashrc` (or `~/.zshrc`):

```bash
export PATH="/path/to/Phyju/SequenceTranslator/SequenceTranslator_binary/bin:$PATH"
```

Then reload:

```bash
source ~/.bashrc
```

> **Note:** The entire `SequenceTranslator_binary/` directory must remain intact — do not copy just the `SequenceTranslator` executable alone. The binary is platform-specific and must be compiled on the target architecture (Linux x86_64, macOS ARM, etc.).

## Usage

Both the script and the binary accept the same arguments:

```bash
SequenceTranslator -i <infile> -o <outfile> -f <format>
# or
julia seqtrans.jl -i <infile> -o <outfile> -f <format>
```

### Arguments

| Flag | Description | Required |
|------|-------------|----------|
| `-i` / `--input` | Input alignment file (any supported format) | Yes |
| `-o` / `--output` | Output file name | Yes |
| `-f` / `--outformat` | Output format: `f` = FASTA, `n` = Nexus, `p` = Phylip | Yes |

### Example — single file

```bash
SequenceTranslator -i alignment.fasta -o alignment.phy -f p
```

### Example — batch conversion of many files (where the binary shines)

```bash
for f in *.fasta; do
    SequenceTranslator -i "$f" -o "${f%.fasta}.phy" -f p
done
```

## Notes

- `SequenceTranslator_binary/` should be added to `.gitignore` — it is large and platform-specific
