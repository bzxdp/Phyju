# ConcatSeqs

Concatenates multiple gene alignments into a single alignment and writes a partition file recording the boundaries of each gene.

Two usage options are available:
- **Script**: run `ConcatAlignments.jl` directly with Julia (requires Julia and dependencies installed)
- **Binary**: compile `ConcatSeqs` for faster startup and no Julia installation required for end users

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
| [Glob.jl](https://github.com/vtjnash/Glob.jl) | Julia General Registry |

## Option 1: Run as a script

### 1. Install dependencies

Start Julia and run:

```julia
using Pkg

# Add unregistered dependency first
Pkg.develop(url="https://github.com/bzxdp/PhyloDataIO")

# Add registered dependencies
Pkg.add(["ArgParse", "Glob"])
```

### 2. Run the script

```bash
julia /path/to/Phyju/ConcatAlignments.jl -e <extension> -o <outfile> -f <format>
```

## Option 2: Compile as a binary

### 1. Clone this repository

```bash
git clone https://github.com/bzxdp/Phyju.git
cd Phyju/ConcatSeqs
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
create_app(".", "ConcatSeqs_binary")
```

This will take approximately 5–10 minutes. The compiled binary will be at:

```
ConcatSeqs_binary/bin/ConcatSeqs
```

### 7. Add the binary to your PATH

Add the following to your `~/.bashrc` (or `~/.zshrc`):

```bash
export PATH="/path/to/Phyju/ConcatSeqs/ConcatSeqs_binary/bin:$PATH"
```

Then reload:

```bash
source ~/.bashrc
```

> **Note:** The entire `ConcatSeqs_binary/` directory must remain intact — do not copy just the `ConcatSeqs` executable alone. The binary is platform-specific and must be compiled on the target architecture (Linux x86_64, macOS ARM, etc.).

## Usage

Both the script and the binary accept the same arguments:

```bash
ConcatSeqs -e <extension> -o <outfile> -f <format>
# or
julia ConcatAlignments.jl -e <extension> -o <outfile> -f <format>
```

### Arguments

| Flag | Description | Required |
|------|-------------|----------|
| `-e` / `--extension` | Extension of gene alignment files to concatenate (e.g. `fas`, `fasta`, `phy`) | Yes |
| `-o` / `--outfile` | Output file name | Yes |
| `-f` / `--format` | Output format: `f` = FASTA, `p` = Phylip, `n` = Nexus | Yes |

### Example

```bash
# Concatenate all .fas files in the current directory
ConcatSeqs -e fas -o concatenated.phy -f p
```

This will produce:
- `concatenated.phy` — the concatenated alignment in Phylip format
- `concatenated_partitions.txt` — partition boundaries for each gene

### Partition file format

```
gene1 = 1 - 350;
gene2 = 351 - 720;
gene3 = 721 - 1100;
```
