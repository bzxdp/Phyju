# Long Branch Remover — how to run it

> **Preliminary.** This document was produced with Claude AI during the debugging and
> revision of the Julia implementation, which was itself carried out with Claude AI. It
> describes code that is still being tested, and should be treated as a working record
> rather than a settled specification. Anything here may change as testing continues.

Two scripts, run in this order:

1. `long_branches_statistical_analyses.jl` — plots the distributions so you can choose K.
2. `long_branches_identifier.jl` — applies the rules and writes what could be pruned.

Both read every tree in the **current working directory**, so `cd` into the directory
holding the trees and give the scripts by full path. Neither prunes anything: the
identifier writes a list of taxa per tree and you decide what to do with it.

---

## 1. Input files

### Clade definitions (`-c`, required)

One clade per line, name first, then its taxa, all whitespace-separated. No colon
after the name — that is the Python prototype's format, not this one.

```
CNID  acropora_palmata pocillopora_verrucosa ricordea_florida nematostella_vectensis
CNID_A acropora_palmata pocillopora_verrucosa ricordea_florida
```

### Clades excluded from quartet tests (`-x`, optional)

Same format. These are clades expected to be naturally long-branched. Members are
masked so the quartet rule compares them against each other rather than against the
rest of the tree.

```
CTENO bathyctena_chuni beroe_sp haeckelia_rubra mnemiopsis_leydi
HEX aphrocallistes_beatrix bolosoma_cynae oopsacas_minuta
```

### K triggers (`-t`, required) — **two different formats**

**For the identifier** — one row of **four** numbers, in this order:

```
3 3 3 3
```

| position | rule |
|---|---|
| 1 | terminal taxa |
| 2 | internal branches / long subtrees |
| 3 | named clade stems |
| 4 | quartet |

**For the statistics script** — one **named row per rule**, the first number being the
value to use and the rest alternatives to plot for comparison:

```
terminals 3 1 2 4 5
internals 3 1 2 4 5
stems     3 1 2 4 5
quartets  3 1 2 4 5
```

Order of the rows does not matter there — they are looked up by name. Fractional
values (3.45, 2.75) are accepted in both.

### Safeguard quantiles (`-s`, identifier only, required)

One row of **three** numbers — not four:

```
0.95 0.95 0.95
```

| position | applies to |
|---|---|
| 1 | terminal taxa **and** the quartet rule (one value serves both) |
| 2 | internal branches / long subtrees |
| 3 | named clade stems |

---

## 2. Choosing the parameters

```bash
cd /path/to/trees
julia /path/to/Long_Branch_Remover/long_branches_statistical_analyses.jl \
    -e treefile \
    -c known_clades.txt \
    -x known_long_branches.txt \
    -t trig.txt
```

Writes violin plots and histograms of the terminal, stem, internal and quartet-ratio
distributions, with the default and alternative K values marked, plus quantile lines at
25/50/75 (grey), 95 (green) and 99 (orange). Read your K and safeguard values off these.

**Note there is no `-s` here.** This script does not apply safeguards, it shows you
candidates for them: the quantile lines are drawn at fixed positions and you pick one to
put in `safe.txt` for the identifier. The plotted quantiles are hard-coded, so if you
want to consider a value other than 0.95 or 0.99 there is no line for it.

---

## 3. Identifying long branches

```bash
cd /path/to/trees
julia /path/to/Long_Branch_Remover/long_branches_identifier.jl \
    -e treefile \
    -o .out \
    -c known_clades.txt \
    -x known_long_branches.txt \
    -t trig.txt \
    -s safe.txt \
    -p 3 \
    -r 1.2
```

| flag | meaning |
|---|---|
| `-e` | tree file extension |
| `-o` | extension appended for the output list (`gene1.treefile` gives `gene1.treefile.out`) |
| `-c` | clade definitions — **required**, the stem rule runs on every tree |
| `-x` | clades excluded from quartet tests — optional |
| `-t` | four K triggers |
| `-s` | three safeguard quantiles |
| `-p` | minimum stacked long internal branches — **required, no default** |
| `-r` | side asymmetry ratio — **required, must be greater than 1.0** |

### `-p` and `-r`

`-p` is how many long internal branches must be stacked in a chain before the clade
beyond them is removed. **1 means any single long internal branch triggers removal** —
very liberal. 2 or 3 is normally more sensible.

`-r` compares the two sides of that branch by their median distance out to their own
tips. 1.2 means one side must be at least 20% deeper than the other before the branch
is taken to say which side is the problem. Higher values (1.3, 1.5) demand a starker
contrast and remove fewer clades. This stops a clade of short-branched taxa being
deleted just because it sits beyond a long stem.

They are independent gates — `-r` does not replace `-p`. A clade is removed only if
the chain is long enough (`-p`), the branch subtending it clears the safeguard
quantile, the two sides differ by more than `-r`, and the deeper side is also the
smaller one.

---

## 4. Outputs

**Per tree**

- `<tree><-o extension>` — taxa nominated for removal, one per line, prefixed by rule:
  - `I:` long subtree (internal branch rule)
  - `S:` named clade stem
  - `T:` terminal branch
  - `Q:` quartet
  A taxon can appear under more than one rule. Strip the prefixes before feeding the
  list to an alignment pruner.
- `<tree>.coloured.nex` — NEXUS with removed clades and taxa coloured. Orange long
  subtree, blue named clade, red terminal, vermillion quartet. Later rules overwrite
  earlier ones, so the internal rule takes precedence over the clade rule.

**Whole run**

- `clade_deletion_frequencies.png` / `.tsv` — how often each named clade was removed by
  either rule, split into the clade itself (blue) and alias deletions (orange, where it
  shared a branch with another clade).
- `clade_deletion_by_internal_rule.png` / `.tsv` — the same, counting only clades lost
  to the internal-branch rule. Compare the pair to see how much the named-clade rule
  adds over what the topology alone already catches.

---

## 5. What the console tells you

- Taxa and clades **not evaluable** because they appear in fewer than 3 trees. These are
  retained — a median from one or two values is noise.
- Per tree, how many taxa the **quartet rule could test**. A low fraction means masking
  is too aggressive for that tree.
- Clades whose stem is long but whose taxa are **already covered by a larger dropped
  clade** — reported once, not twice.
- With 3 trees or fewer, a warning that the terminal, stem and internal rules are
  **disabled** entirely; only the quartet rule runs. Warnings also below 10 and 50 trees.

---

## 6. Errors you might hit

| message | cause |
|---|---|
| `must contain 4 K values` | trigger file has the wrong number of values, or you gave the identifier the named-row format |
| `must contain 3 quantiles` | safeguard file still has the old 4 values — drop the last one |
| `--internal_side_ratio must be greater than 1.0` | `-r` set to 1.0 or below, which disables the asymmetry test |
| `--internal_long_path__trigger must be at least 1` | `-p` set to 0 or negative |
| `required argument --clades was not provided` | the stem rule runs on every tree, so `-c` is never optional |
| `File ... has a node with no branch length` | a tree lacks branch lengths; all rules need them |

---

## 7. Note on the code layout

The rules live in `long_branch_rules.jl`, which both scripts pull in with `include`.
They are kept out of TreeStats and TreeUtils deliberately: those libraries are shared
with other projects and hold only general tree operations, whereas everything in the
rules file takes a K trigger, a safeguard quantile or a clade definition file.

Full detail on the rules and why they work as they do is in
`long_branch_python_vs_julia.pdf` in this directory.
