---
title: "bwa: Pitfalls and gotchas"
type: reference
generated: true
---

> Part of the [bwa](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/bwa/references/pitfalls.md)

{% raw %}
# Pitfalls and gotchas

- **`-c` subsamples, it does not skip.** Usage text (`fastmap.c:278`) and man page (`bwa.1:203-205`) say "skip/discard"; the code keeps `max_occ` evenly spaced SA entries (`bwamem.c:304-305`). The old `continue` is commented out at `303`. The subsample is arbitrary with respect to genomic position.
- **Seeds are maximal matches only.** Shorter k-mers shared by more copies are never emitted by pass 1 (`bwt.c:289-351`); re-seeding (`-r`, hidden `-s`) is one level deep and only for SMEMs with `<= 10` occurrences (`bwamem.c:163`). Pass 3 (`-y`) only emits seeds rarer than `-y` (`bwt.c:370`).
- **Hidden options**: `-s split_width` (parsed `fastmap.c:186`, usage line commented out `277`), `-X mask_level` (`194`), `-Q mapQ_coef_len` (`203-206`), `-F debug` (`195`), `-G max_chain_gap` (`187`), `-N max_chain_extend` (`188`) — all functional.
- **`-A` rescales other penalties** that were not given explicitly (`update_a`, `fastmap.c:125-139`), including `-T` and `-d`. `-x` presets bypass `update_a` (`359`).
- **`-D` has two jobs**: chain-weight dropping (`bwamem.c:384`) and the secondary output threshold `score >= drop_ratio * primary` (`1050`). `-a` without `-D 0` hides secondaries scoring below half the primary.
- **XA is all-or-nothing**: more than `-h` eligible hits -> no XA tag (`bwamem_extra.c:144`); `-a` disables XA (`bwamem.c:1041`).
- **MAPQ for repeats is 0 by construction**: `sub >= score -> 0` (`bwamem.c:987`), then `* (1 - frac_rep)` (`1004`). Secondaries always MAPQ 0 (`1136`), SEQ/QUAL `*` (`897-898`), no `XS` (`1056`).
- **Supplementary lines are hard-clipped** unless `-Y` (`838-849`), so their SEQ is partial.
- **`-T` only filters output** (`1048`); it does not speed anything up.
- **Default scoring is for ~0.3% divergence** (`bwamem.c:57`: `e = .75*exp(-ln4*B/A)`); at <= 80% identity expected score per base is <= 0.
- **Extension is scalar**, not SSE (`ksw_extend2`, `ksw.c:416`); SSE2 `ksw_u8/i16` serve seed-SW and mate rescue only.
- **Index sampling is fixed** by `bwa index` (`bwtindex.c:316`, `sa_intv = 32`); use `bwt2sa -i` manually. `OCC_INTERVAL` must stay 128 (`bwt.h:36`).
- **`-a is` overflows** for `2G >= 2^31` (`is.c:208` takes `int n`); auto-selection uses `l_pac > 50 Mbp` where `l_pac` is already doubled (`bwtindex.c:275`).
- **Both strands are in the BWT**; `x[2]` occurrence counts and `fastmap` counts include both strands. Reference positions `>= l_pac` are reverse strand (`bns_depos`).
- **N handling**: reference N -> random base with fixed seed (`bntseq.c:266,295`); query N -> code 4, terminates SMEM extension (`bwt.c:296,316-319`) and scores -1 (`bwa.c:142`).
- **Seeds crossing a contig or strand boundary are discarded**, not split (`bwamem.c:312-313`, TODO comment).
- **`-v 4` prints debug to stdout** and corrupts SAM (`bwa.1:368-369`).
- **Batch size depends on `-t`** (`chunk_size * n_threads`, `fastmap.c:394`); PE insert-size inference therefore depends on thread count unless `-K` (`bwa.1:316-324`). Output order is preserved by `kt_pipeline`.
- **Coordinate-sorted FASTQ** makes whole batches repetitive and slow/memory hungry (README FAQ 8). Repeat consensus libraries behave the same.
- **Man page staleness**: `bwasw` defaults `-T 37 -w 33` in the man page (`bwa.1:594-597`) vs 30 and 50 in code (`bwtsw2_aux.c:57-58`); `-c` wording as above; `bwasw` is flagged DEPRECATED (`main.c:70`); README now recommends minimap2/minibwa/bwa-mem2/bwa-mem3.
- **`samse -n`**: XA multi-hits only if total occurrences `<= n+1`, otherwise none (`bwase.c:57-58`).
- `mem_reg2aln` comment at `bwa.c:169`: the gap-free fast path is never taken ("hurts performance", not accuracy).

## Defaults that surprise

| Parameter | Default | Why it matters |
|---|---|---|
| `-k` | 19 | 19-mer exact seeds exist in only ~26% of 75%-identity 300 bp copies |
| `-c` | 500 | seeds above this are subsampled, not dropped; frac_rep -> MAPQ shrink |
| `-s` (hidden) | 10 | re-seeding never happens for SMEMs with > 10 copies |
| `-y` | 20 | third pass never emits seeds with >= 20 copies |
| `-D` | 0.5 | also hides secondaries < 50% of the primary score with `-a` |
| `-h` | 5 | > 5 alt hits -> XA omitted entirely |
| `-B` | 4 | implies 0.3% divergence; zero expected score at 80% identity |
| `-T` | 30 | with A=1 ~ 30 net matching bases |
| SA interval | 32 | locate = up to 31 LF steps per occurrence |
| chunk | 10 Mbp x threads | output identical, PE stats not, across `-t` |

## Known bugs / open upstream issues (as of pinned commit)

- TODO at `bwamem.c:313`: seeds bridging reference boundaries are dropped rather than split.
- FIXME at `bwtsw2_aux.c:587`: ambiguous bases in bwasw queries are randomised.
- README states BWA-MEM "will be replaced by minibwa"; the repository is in maintenance mode (last release 0.7.19, March 2025; NEWS says output identical to 0.7.17).
{% endraw %}
