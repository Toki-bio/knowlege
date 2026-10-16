---
title: "sear2k: Performance: where the wall-clock goes"
type: reference
generated: true
---

> Part of the [sear2k](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/sear2k/references/performance.md)

{% raw %}
# Performance: where the wall-clock goes

> Type: reference for the code structure; `notes` for the inferences. This
> repo has no timing data of its own; numbers are from the SINEderella repo
> (`90d8ac4`) which runs this `sear`. See the `sinederella` skill,
> `references/performance.md`, for the full table.

## Work decomposition of one `sear` call (ssearch36 mode)

| Stage | Cost driver | Code | Parallel? |
|---|---|---|---|
| `seqkit sliding` + `split2` | one pass over the genome, written as 2.5 kb records (1.25x the genome on disk) | `sear:372-373` | no; done once per directory when `-k` is used (`step1_search_extract.sh:107-111`) |
| 20 x `ssearch36 -m 8C` | full Smith-Waterman of a ~300 bp query against every 2.5 kb window = O(querylen x genome x 1.25). Plus per-hit work: every window above the 8-bit score limit is re-aligned in 16 bits and then re-aligned for display (`docs/MASKING.md:46`), so a query with 100k copies pays for 100k alignments on top of the scan. | `sear:381-454` | ssearch36 `-T` threads inside one part; parts are serial |
| awk/sort/bedtools per part | linear in hits | `sear:397-411` | no |
| minus rounds (`-m`) | repeats the scan on the complement of found loci; rounds until none found | `sear:438-449` | serial |
| coordinate mapping + final merge | linear in hits | `sear:456-499` | no |
| `bedtools getfasta` | random access via `.fai` | `sear:582` | no |
| MAFFT (`-b`) | L-INS-i on N+1 sequences, N <= 100 | `sear:562-564` | mafft threads |

## Measured numbers that involve `sear`

| Observation | Source |
|---|---|
| Rhin-1 (short bat SINE) vs 2.3 Gb *Tadarida*: **69 s per query** unmasked, 51 s with the 136 Mb VES family masked to N; MEG queries 51-64 s. "the Smith-Waterman scan over the genome costs the same, so an unrelated query is not faster" | SINEderella `docs/MASKING.md:121-132` |
| Planted 8 Mb genome: 37-46 s per query when the 120 bp family has 1,500 copies, ~4 s when they are N-masked: "the time is the known copies being aligned", not the scan | `docs/MASKING.md:19-24,44-46` |
| hg38 chr21 (46.7 Mb), 2 consensuses: whole SINEderella full run 20 min; `--add` of one consensus 10 min | `docs/RELEASE_CHECK_2026-10-05.md:22-23` |
| rsi (2.1 Gb), 21 consensuses: step 1 (21 `sear` calls + satellite stage) 61 min of an 84 min run | `docs/RELEASE_CHECK_2026-10-05.md:75` |
| tbr (2.3 Gb): 54 min full run with a 2-consensus bank | `docs/RELEASE_CHECK_2026-10-05.md:101` |
| `sear_multi` (all queries in one ssearch36 pass per part) exists to cut the per-part overhead; `benchmark_sear.sh` measures old vs new wall time and BED identity. No result table is committed. | SINEderella `sear_multi:4-7`, `benchmark_sear.sh:15-20,279-280` |

`notes`: with ~50-70 s per consensus on a 2 Gb genome, a 20-consensus bank is
~20 min of ssearch36 and a 60-consensus bank ~1 h, scaling linearly with the
bank size because queries are searched one at a time. The per-copy
re-alignment term means hit-rich families (VES: 620k copies) dominate.

## Parallelism

- One part at a time; `ssearch36 -T` is the only knob (`sear:396`). On a
  32-core box the 20 parts are each searched with 32 threads, so the serial
  loop is not itself the bottleneck, but ssearch36's threading efficiency on a
  single 300 bp query is unknown (unverified).
- `step1_search_extract.sh:109-123` loops over consensuses serially too. The
  multi-genome wrapper `SINEderella_multi --parallel N` runs whole genomes in
  parallel (`SINEderella_multi:31-32`).
- Nothing in this repo uses GNU parallel or `xargs -P` except
  `annotate_aligned_by_ssearch_seqkit.sh:51-52`.

## Author's comments about speed and sensitivity (verbatim or near)

- `README.md:5`: minimap2 for >= 1000 bp "for faster whole-genome alignment (no fragment scanning needed)".
- SINEderella `sear:36-37`: "asm5 is too strict for diverged copies".
- `docs/MASKING.md:5,56`: "Keep ssearch36 ... Do not move the search to blastn ... blastn: default word 11 loses copies below ~85% identity; blastn-short keeps them only with a short-word seed, the class of miss full Smith-Waterman was chosen to avoid. minimap2: no SINE-scale sensitivity and no case."
- `MANUAL.md:1350-1354`: the 0.8-length / 65 %-identity rule is "structurally immune to the tandem-repeat false-positive mode" that k-mer tools need filters for.
- `docs/RELEASE_CHECK_2026-10-05.md:30`: "Recall is bounded by the length rule, not by divergence": missed RepeatMasker Alus are fragments < 240 bp, not more diverged copies.
- `flankscan/HANDOFF.md:150-153`: "ssearch36 loses hits when the library is huge and one family dominates" (622k windows, 81k hits, 91 % of VES copies without a core hit); fixed there by 20,000-record chunks with `-Z 20000`. `notes`: `sear`'s 20 parts of a 2 Gb genome are ~50k windows each and have no `-Z`; whether the same loss occurs inside `sear` is untested.
- No `TODO`/`FIXME` in this repo. Only speed-related TODO in SINEderella touching the search: `docs/SATELLITES.md:193,206` ("use raw ssearch36 hits from step 1", since done via `rawhits.tsv`).
{% endraw %}
