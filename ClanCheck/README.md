# ClanCheck

A Julia implementation of ClanCheck as described in:

> Siu Ting et al. (2019) *Into the Shallow: Discordance, Incomplete Lineage Sorting, and Ancient Admixture Confound Phylogenetic Systematics of the Short-tailed Opossums (Didelphidae; Monodelphis)*. Molecular Biology and Evolution, 36(6):1344–1355. https://doi.org/10.1093/molbev/msz067

Given a set of gene trees and a set of user-defined clades, ClanCheck identifies which gene trees have all specified clades monophyletic. Trees that fail any clade constraint are excluded from the output.

If you use this software please cite Siu Ting et al. (2019) and this repository.

Two usage options are available:
- **Script**: run `clancheck.jl` directly with Julia — sufficient for occasional analyses
- **Binary**: compile `ClanCheck` for batch analyses over large numbers of gene trees, where Julia's per-invocation compilation overhead becomes significant

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
| [TreesIO](https://github.com/bzxdp/TreesIO) | GitHub (unregistered) |
| [TreeUtils](https://github.com/bzxdp/TreeUtils) | GitHub (unregistered) |
| [ArgParse.jl](https://github.com/carlobaldassi/ArgParse.jl) | Julia General Registry |
| [MetaGraphsNext.jl](https://github.com/JuliaGraphs/MetaGraphsNext.jl) | Julia General Registry |

## Option 1: Run as a script

### 1. Install dependencies

Start Julia and run:

```julia
using Pkg
Pkg.develop(url="https://github.com/bzxdp/TreesIO")
Pkg.develop(url="https://github.com/bzxdp/TreeUtils")
Pkg.add(["ArgParse", "MetaGraphsNext"])
```

### 2. Run the script

```bash
julia /path/to/Phyju/ClanCheck/clancheck.jl -e <extension> -o <outfile> -c <clades_file>
```

## Option 2: Compile as a binary

### 1. Clone this repository

```bash
git clone https://github.com/bzxdp/Phyju.git
cd Phyju/ClanCheck
```

### 2. Start Julia and activate the project

```julia
using Pkg
Pkg.activate(".")
```

### 3. Add unregistered dependencies

```julia
Pkg.develop(url="https://github.com/bzxdp/TreesIO")
Pkg.develop(url="https://github.com/bzxdp/TreeUtils")
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
create_app(".", "ClanCheck_binary")
```

This will take approximately 5–10 minutes. The compiled binary will be at:

```
ClanCheck_binary/bin/ClanCheck
```

### 7. Add the binary to your PATH

Add the following to your `~/.bashrc` (or `~/.zshrc`):

```bash
export PATH="/path/to/Phyju/ClanCheck/ClanCheck_binary/bin:$PATH"
```

Then reload:

```bash
source ~/.bashrc
```

> **Note:** The entire `ClanCheck_binary/` directory must remain intact — do not copy just the `ClanCheck` executable alone. The binary is platform-specific and must be compiled on the target architecture (Linux x86_64, macOS ARM, etc.).

## Usage

Both the script and the binary accept the same arguments:

```bash
ClanCheck -e <extension> -o <outfile> -c <clades_file>
# or
julia clancheck.jl -e <extension> -o <outfile> -c <clades_file>
```

### Arguments

| Flag | Description | Required |
|------|-------------|----------|
| `-e` / `--extension` | Extension of tree files (e.g. `tree`, `tre`, `treefile`) | Yes |
| `-o` / `--output` | Output file name | Yes |
| `-c` / `--clades` | File defining clades to test (see format below) | Yes |

### Clades file format

One clade per line. The clade name comes first, followed by all taxa in the clade, all separated by spaces. **No spaces allowed in clade or taxon names** — use underscores instead:

```
Clade_A Taxon_1 Taxon_2 Taxon_3
Clade_B Taxon_4 Taxon_5
```

### Example

```bash
ClanCheck -e treefile -o passing_trees.txt -c clades.txt
```

### Output

A plain text file listing one tree filename per line — only trees where **all** specified clades are monophyletic are included. Trees failing any clade constraint are excluded.

### Notes

- ClanCheck operates on all files with the specified extension in the **current working directory**
- If a clade contains taxa not present in a given gene tree, that constraint is automatically treated as a pass for that tree (the clade cannot be falsified by absent taxa)
- If a clade contains only one taxon present in a tree, it is also treated as a pass (a single taxon is trivially monophyletic)
- `ClanCheck_binary/` should be added to `.gitignore` — it is large and platform-specific

## Citation

If you use ClanCheck please cite:

Siu Ting et al. (2019) Into the Shallow: Discordance, Incomplete Lineage Sorting, and Ancient Admixture Confound Phylogenetic Systematics of the Short-tailed Opossums (Didelphidae; Monodelphis). *Molecular Biology and Evolution*, 36(6):1344–1355. https://doi.org/10.1093/molbev/msz067

And this repository: https://github.com/bzxdp/Phyju
