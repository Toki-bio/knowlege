# Aligner invocations (ssearch36, minimap2, mafft) and how their output is parsed

All line numbers are for `sear` at `85421e7` unless prefixed with another file.

## 1. ssearch36 (fragment mode)

### Exact command (`sear:396`)

```sh
ssearch36 -g -3 -T "${THREADS:-$(nproc 2>/dev/null || echo 1)}" -Q -n -z 11 -E 2 -w 95 -W 70 -m 8C "$query" "$bank" 3
```

`sear2k:190` and `sear2k_bak:58` are identical except `-T 24` hard-coded.
`$query` = the single-sequence query FASTA; `$bank` = one `*.2k.part_NNN.bnk`
file (1/20th of the 2.5 kb windows), or in minus rounds `minus_<round>_<part>`.

### Flag-by-flag (FASTA36 semantics; **not run-verified here**, ssearch36 is not installed in this container)

| Token | FASTA36 meaning (my reading of the fasta36 manual) | What SINEderella's MANUAL.md:1218-1230 says | Note |
|---|---|---|---|
| `-g -3` | `-g #` is the **gap-extension penalty**; here it consumes `-3` as its argument, so gap-extend = -3 (DNA default -4). | `-g` = "Global/gapped alignment", `-3` = "Forward strand + reverse complement" | `notes`: if my reading is right, `-3` as a standalone flag (forward strand only) is never set, both strands are searched because that is the DNA default, and gap-open stays at the default (-16 for DNA, `-f`). Verify with `ssearch36 -g -3 ... ` vs `ssearch36 -g=-3 ...` on a toy before relying on either reading. Either way both strands ARE searched: the awk at `sear:399-400` sees minus-strand hits (qend < qstart) and the result tables contain `-` strands. |
| `-T N` | threads for the library scan | thread count | `THREADS` env or `nproc` |
| `-Q` | quiet (no interactive prompts) | quiet | |
| `-n` | query is DNA | nucleotide | |
| `-z 11` | statistics from shuffled library sequences (regression on shuffles) | "shuffle with composition" | makes E-values and bitscores **stochastic run to run** (SINEderella `MANUAL.md:1212`: +-1-3 %; `docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:43-44`: 99.37 % of firm/not-firm outcomes reproduce) |
| `-E 2` | report hits with E <= 2 | E-value threshold | relative to the **part file** (~1/20 of the windows), no `-Z` given; so the threshold means something different for a 100 Mb and a 3 Gb genome. SINEderella's flankscan fixes this with `-Z 20000` (`flankscan/fs3_partners.sh:95`), sear does not. |
| `-w 95` | alignment display line width | "alignment width" | no effect with `-m 8C` (tabular) |
| `-W 70` | context residues shown around an alignment | "window size" | no effect with `-m 8C` |
| `-m 8C` | BLAST tabular with comment lines (`#` lines) | tabular | parsed by `/^[^#]/` |
| trailing `3` | ktup positional; ssearch36 (full Smith-Waterman) does not use ktup | not mentioned | harmless |

No scoring matrix (`-s`), no gap-open (`-f`), no `-Z`, no `-d`/`-b` limits on
the number of reported hits, no `-S` (lowercase masking; SINEderella's
`docs/MASKING.md:36-42` explains that `-S` does nothing for DNA libraries in
the shipped build anyway).

### Parsing (`sear:397-401`)

`-m 8C` columns: 1 query, 2 subject, 3 %identity, 4 aln length, 5 mismatches,
6 gap opens, 7 qstart, 8 qend, 9 sstart, 10 send, 11 evalue, 12 bitscore.

```awk
/^[^#]/ { printf($2"\t"$9-1"\t"$10"\t"$10-$9+1"\t");           # window, 0-based start, end, hit length
          if ($8-$7>0) print $8-$7+1"\t+\t"$3"\t"$12;           # query span, strand +, %id, bits
          else         print -($8-$7+1)"\t-\t"$3"\t"$12 }       # minus strand when qend < qstart
```

Produces an 8-column window-BED: `window start end hitlen qspan strand pctid bits`.
`sort -V` -> `bedtools merge -s -c 4,5,6,7,8 -o max,max,distinct,max,max`
(`sear:405`) collapses overlapping hits on the same strand within a window and
keeps the **max** hit length, max query span, max identity, max bitscore. Note
the merged output has 8 columns: `chrom start end strand max4 max5 distinct6 max7 max8`
is NOT the layout; bedtools merge -s with -c prints `chrom start end` then the
`-c` columns in order, so columns are `1 win 2 start 3 end 4 hitlen 5 qspan 6 strand 7 pctid 8 bits`.
Filter (`sear:407`): `$7 >= Homology && $4 >= flen` where `flen = int(Slf*len)`
(`sear:245-246`) - i.e. identity >= 65 and **max hit length** (col 4) >= 0.8 x query.
`notes`: the commit message of `85421e7` says the filter was moved to "merged
span (\-\)"; at this SHA col 4 is `max` of per-hit lengths, not `$3-$2`. In
minimap2 mode the span `$3-$2` is used (`sear:337`). Treat the ssearch36
branch's length rule as "longest single alignment >= 0.8 x query" until
verified on a run.

## 2. minimap2 (query >= 1000 bp)

Command (`sear:295-298`), inside `for preset in asm10 asm20`:

```sh
minimap2 -x "$preset" -t "${THREADS:-nproc}" --secondary=yes -p 0.3 -N 10000 "$WORKING_BANK" "$query"
```

Target = whole genome (cleaned copy), query = the consensus. PAF parsed at
`sear:299-319`: strand `$5`, target `$6 $8 $9`, matches `$10`, `AS:i:` tag;
identity = `matches / min(qspan, tspan) * 100` ("to penalize gaps", `sear:309`);
output `tname tstart tend pct_id tlen strand AS` -> `bedtools sort | merge -s -c 4,5,6,7 -o max,max,distinct,max`
-> filter `$4 >= Homology && ($3-$2) >= flen` (`sear:334-338`). The cascade
stops at the first preset with **any raw alignment** (`sear:321-326`), before
filtering, so asm20 is only tried when asm10 returns zero lines.

## 3. mafft (`-b` mode, `sear:562-564`)

```sh
mafft --thread $T --threadtb $T --threadit $T --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --quiet in.fa
```

Input = query + top-N flanked hits. SINEderella later found that
`--threadit > 0` makes iterative refinement non-reproducible and switched its
own calls to `--threadit 0` (`PUBLISH_WORKFLOW.md:123`); this repo still passes
`--threadit $T`.

## 4. annotate_aligned_by_ssearch_seqkit.sh

`ssearch36 -m 8 "$f" "$db"` per single-sequence file, parallel via
`xargs -0 -n1 -P $THREADS` (`:44-52`); best hit = max of the last column.
Default statistics (no `-z 11`), no `-E`.
