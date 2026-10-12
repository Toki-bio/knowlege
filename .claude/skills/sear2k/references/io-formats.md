# Input and output formats

## Inputs accepted

| Input | Constraints | Where read |
|---|---|---|
| `query.fasta` | One SINE consensus (100-500 bp typical). Length is summed over all records (`sear:245`), so pass one record per file; SINEderella splits its bank with `seqkit split -i` first. `sear_multi` (SINEderella) is the multi-record variant. | `sear:146,245` |
| `bank.fasta` | Genome FASTA. Headers with `_` or tabs trigger a cleaned copy (`_` -> `@U@`) (`sear:180-214`). `.fai` created if missing. In SINEderella the genome has already been sanitised (`_` -> `@U@`, uppercase, 60-col) so no copy is made. | `sear:147,173-214` |
| `Slf` (positional 3) | fraction of query length a hit must reach, default 0.8 (`sear:216`); 0.9 in `sear2k:146` | |
| `Homology` (positional 4) | minimum % identity, default 65 | `sear:217` |
| `Flank` (positional 5) | flank added to `.bnk` extraction, default 50; `--slop L,R` overrides it and the best-hit flank | `sear:218,221-231` |
| `KeepSplits` (positional 6) / `-k` | keep `*.2k.*` and the clean bank for the next query | `sear:219,233-236` |
| `THREADS` (env) | threads for ssearch36/minimap2/mafft | `sear:296,396,562` |

Typical sizes (from SINEderella runs, `docs/RELEASE_CHECK_2026-10-05.md`):
genomes 46.7 Mb (hg38 chr21) to 2.1-3.0 Gb (bats, Sicista); banks of 2-58
consensuses of 105-600 bp; one `sear` call per consensus.

## Outputs produced

| File | Writer | Columns / fields |
|---|---|---|
| `<tax>-<qname>.bed` (BED7) | `sear:482-499` (ssearch36), `sear:334-338` (minimap2) | `chrom start end pctid length strand bitscore`. `pctid` = max identity of merged hits; `length` = max hit length (ssearch36) or target span (minimap2); `bitscore` = max ssearch36 bits or minimap2 `AS`. Strand-aware merged, 0-based half-open. |
| `<tax>-<qname>.bnk` | `sear:581-582` | FASTA of each BED interval +- flank, `bedtools getfasta -s` headers `chrom:start-end(strand)` |
| `<tax>-<qname>.best<N>.mafft.fa` | `sear:562-564` | MAFFT alignment of query + top-N flanked hits (names `chrom:start-end(+/-)`), only with `-b` |

`<tax>` = first three characters of the bank file name (in ssearch36 mode taken
from the part file names, `sear:387`; SINEderella always sees `gen` because its
genome is `genome.clean.fa`, hence `searches/gen-*.bed` everywhere in that
pipeline). `<qname>` = query path minus extension (`sear:250`), so for
SINEderella `gen-consensuses.clean.part_<NAME>.bed`.

## Intermediate files (cwd, removed by `cleanup_query` unless noted)

| File | Meaning |
|---|---|
| `<bank>.2k.bnk`, `<bank>.2k.bnk.fai` | all sliding windows; headers `<chrom>_sliding:<start>-<end>` (seqkit `-g` greedy last window) |
| `<bank>.2k.part_001.bnk` ... `part_020.bnk` | 20 equal-record splits; **kept** with `-k` |
| `<bank>_temp_clean.<ext>` (+ `.fai`) | header-cleaned copy; kept with `-k` |
| `_ssearch_sorted_<round>_<part>.tmp` | parsed, sorted window hits of one ssearch36 call |
| `plus_to_add_<round>_<part>.bed` | merged + filtered hits of one round, window coordinates (8 cols) |
| `<round>_<part>.bed`, `minus_<round>_<part>` | cumulative hits and the complement FASTA for the next minus round |
| `<bank>.all_hits.txt` | concatenation of all `plus_to_add_*` |
| `<bank>.all` | one offset per line of `all_hits.txt` |
| `<bank>.stripped.genome` | chrom lengths without `_sliding` suffixes |
| `<bank>.minimap_raw.txt` | parsed PAF (minimap2 mode) |
| `<tax>-<qname>.best<N>.{bed7,bed,flank50.*.bed,fa,mafft.in.fa,mafft.err}` | best-mode temporaries, deleted at `sear:614-623` |

SINEderella's copy additionally writes `<tax>-<qname>.rawhits.tsv` (every merged
hit passing the identity cut, any length, genome coordinates) for its satellite
stage (SINEderella `sear:341,387,417,540-544`); this repo's `sear` does not.

## Window header -> genome offset

Window headers look like `NC_012345.1_sliding:2001-4500`; because `_` in the
chromosome name is already `@U@`, the only `_` left is `_sliding`. The offset
awk (`sear:456-460`) replaces `_sliding:` by a space, takes field 2
(`2001-4500`), turns `-` into ` -` (so the end becomes `-4500`), and sums only
positive fields: offset = 2001, the 1-based window start. The hit's 0-based
start `$9-1` and 1-based end `$10` plus 2001 are both one too high;
`bedtools slop -l 1 -r -1` (`sear:497`) shifts them back. The `sub(/_slidin.*/)`
... `sub(/_sl.*/)` chain (`sear:466-472,487-493`) exists because ssearch36
truncates long identifiers, leaving partial suffixes. `notes`: read, not run;
test on a hand-written header before reusing this parser in a replacement tool.
